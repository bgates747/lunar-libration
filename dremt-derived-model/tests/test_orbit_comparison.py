"""Verify the ephemeris input and Earth-centred virtual-Moon geometry."""
import json
import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from orbit_comparison import build_comparison_track, build_perigee_comparison_track, SECONDS_PER_DAY


@pytest.fixture(scope="module")
def comparison():
    return build_comparison_track()


@pytest.fixture(scope="module")
def perigee_track():
    return build_perigee_comparison_track()


def test_real_moon_uses_later_ephemeris_nodes(comparison):
    reference = json.loads((Path(__file__).resolve().parents[1] / "data" / "moon-model-reference.json").read_text())
    for index in (0, 199, 500, 672):
        day = comparison.orbit.times_seconds[index] / SECONDS_PER_DAY
        sample = comparison.sample(day)
        np.testing.assert_allclose(sample.position_km, [float(v) for v in reference["moon"]["rows"][index][2:5]], atol=1e-9)
        np.testing.assert_allclose(sample.velocity_km_s, [float(v) for v in reference["moon"]["rows"][index][5:8]], atol=1e-12)
        np.testing.assert_allclose(sample.attitude.T @ sample.attitude, np.eye(3), atol=1e-13)
        assert sample.subearth_latitude_deg == pytest.approx(float(reference["orientation"]["rows"][index][4]))


def test_virtual_rate_is_constant_while_real_rate_varies(comparison):
    days = np.linspace(0, 28, 57)
    samples = [comparison.sample(float(day)) for day in days]
    angles = np.unwrap([math.atan2(s.projected_virtual_position_km[1], s.projected_virtual_position_km[0]) for s in samples])
    expected = comparison.mean_angular_rate_rad_s * SECONDS_PER_DAY
    np.testing.assert_allclose(np.diff(angles) / np.diff(days), expected, atol=1e-13)
    assert math.degrees(expected) == pytest.approx(360 / 27.321661)
    real_rates = [s.real_angular_rate_rad_s for s in samples]
    assert max(real_rates) > min(real_rates) * 1.1
    assert samples[0].angular_separation_rad == pytest.approx(0, abs=1e-15)


def test_earth_centred_angle_and_radial_distance_are_independent(comparison):
    for day in (0, 0.27123, 7, 13.7, 28):
        s = comparison.sample(day)
        real, virtual = s.projected_position_km, s.projected_virtual_position_km
        expected = math.atan2(real[0]*virtual[1] - real[1]*virtual[0], np.dot(real, virtual))
        assert s.angular_separation_rad == pytest.approx(expected, abs=1e-14)
        assert np.linalg.norm(real) == pytest.approx(np.linalg.norm(virtual), abs=1e-9)
        assert np.linalg.norm(s.position_km) == pytest.approx(np.linalg.norm(s.virtual_position_km), abs=1e-9)
        np.testing.assert_allclose(comparison.orbit.camera_basis[:, :2].T @ s.virtual_position_km, virtual, atol=1e-9)
        assert np.dot(comparison.orbit.camera_basis[:, 2], s.position_km) == pytest.approx(np.dot(comparison.orbit.camera_basis[:, 2], s.virtual_position_km), abs=1e-9)


def test_phase_changes_only_virtual_bearing_and_no_extrapolation(comparison):
    shifted = build_comparison_track(virtual_phase_rad=0.2)
    for day in (0, 4.2, 28):
        original, offset = comparison.sample(day), shifted.sample(day)
        np.testing.assert_array_equal(original.position_km, offset.position_km)
        assert offset.angular_separation_rad - original.angular_separation_rad == pytest.approx(0.2, abs=1e-14)
    for invalid in (-0.01, 28.01, float("nan")):
        with pytest.raises(ValueError):
            comparison.sample(invalid)


def test_perigee_start_is_radial_minimum_on_positive_x_axis(perigee_track):
    track = perigee_track
    first = track.sample(0)
    assert track.start_epoch.startswith("2026-09-06")
    assert first.projected_position_km[0] > 0
    assert first.projected_position_km[1] == pytest.approx(0, abs=1e-9)
    assert first.geometric_angle_rad == pytest.approx(0, abs=1e-13)
    assert np.dot(first.position_km, first.velocity_km_s)/first.distance_km == pytest.approx(0, abs=1e-10)
    for offset in (-600, 600):
        nearby = track.orbit._state_at_seconds(track.start_seconds+offset)
        assert np.linalg.norm(nearby[:3]) > first.distance_km
    np.testing.assert_allclose(track.orbit.camera_basis.T @ track.orbit.camera_basis, np.eye(3), atol=1e-13)


def test_full_orbit_reaches_360_degrees_and_retains_later_apsidal_shift(perigee_track):
    track = perigee_track
    samples = [track.sample(float(day)) for day in np.linspace(0, track.duration_days, 57)]
    assert np.all(np.diff([s.geometric_angle_rad for s in samples]) > 0)
    assert samples[-1].geometric_angle_rad == pytest.approx(2*math.pi, abs=1e-13)
    assert samples[-1].projected_position_km[1] == pytest.approx(0, abs=1e-9)
    next_perigee = [event for event in track.apses if event['kind']=='perigee' and event['elapsed_days'] > 0][0]
    assert next_perigee['elapsed_days'] < track.duration_days
    assert next_perigee['geometric_angle_deg'] < 350
    # The full-revolution endpoint is allowed to have nonzero radial velocity;
    # it must not be artificially forced to be another perigee.
    final = samples[-1]
    assert abs(np.dot(final.position_km, final.velocity_km_s)/final.distance_km) > 0.001


def test_perigee_selection_preserves_states_mean_rate_and_ephemeris_departures(perigee_track):
    track = perigee_track
    node = len(track.orbit.times_seconds)//2
    day = (track.orbit.times_seconds[node]-track.start_seconds)/SECONDS_PER_DAY
    sample = track.sample(day)
    np.testing.assert_allclose(sample.position_km, track.orbit.states[node,:3], atol=1e-8)
    np.testing.assert_allclose(sample.velocity_km_s, track.orbit.states[node,3:6], atol=1e-12)
    first, last = track.sample(0), track.sample(track.duration_days)
    assert (last.virtual_angle_rad-first.virtual_angle_rad)/track.duration_days == pytest.approx(2*math.pi/27.321661)
    assert first.angular_separation_rad == pytest.approx(0, abs=1e-14)
    apogee = [event for event in track.apses if event['kind']=='apogee'][0]
    assert abs(apogee['angular_separation_deg']) > 0.1
    assert abs(apogee['geometric_angle_deg']-180) > 1


def test_short_cache_rejects_incomplete_revolution():
    path = Path(__file__).resolve().parents[1] / 'data' / 'moon-model-reference.json'
    with pytest.raises(ValueError, match='full revolution'):
        build_perigee_comparison_track(path)
