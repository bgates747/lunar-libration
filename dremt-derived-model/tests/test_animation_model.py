"""Check the 3D geometry that the top-down renderer must preserve."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from animation_model import AnimationTrack, SECONDS_PER_DAY, build_track
import moon_model


@pytest.fixture(scope="module")
def monthly_track():
    return build_track()


def test_inclined_pole_keeps_latitude_and_distinguishes_projected_angle():
    tilt = math.radians(30)
    attitude = np.array([
        [math.cos(tilt), 0.0, math.sin(tilt)],
        [0.0, 1.0, 0.0],
        [-math.sin(tilt), 0.0, math.cos(tilt)],
    ])
    radius = 384400.0
    position = radius / math.sqrt(2) * np.array([1.0, 1.0, 0.0])
    velocity = np.array([-1.0, 1.0, 0.0]) / math.sqrt(2)
    state = np.concatenate((position, velocity, [1.5e8, 0.0, 0.0], [0.0, 30.0, 0.0]))
    track = AnimationTrack(
        times_seconds=np.array([0.0, 3600.0]),
        states=np.array([state, state + 3600 * np.concatenate((velocity, np.zeros(3), [0.0, 30.0, 0.0], np.zeros(3)))]),
        initial_attitude=attitude,
        camera_basis=np.eye(3),
        start_epoch="test",
        end_epoch="test",
    )
    sample = track.sample(0)
    expected_latitude = math.degrees(math.asin(-math.sin(tilt) / math.sqrt(2)))
    expected_longitude = math.degrees(math.atan2(-1.0, -math.cos(tilt)))
    assert sample.subearth_latitude_deg == pytest.approx(expected_latitude)
    assert sample.subearth_longitude_deg == pytest.approx(expected_longitude)
    assert sample.projected_rocking_deg == pytest.approx(-135.0)
    assert abs(sample.projected_rocking_deg - sample.subearth_longitude_deg) > 4
    np.testing.assert_allclose(sample.baseline @ sample.relative, attitude, atol=1e-15)


def test_monthly_geometry_is_proper_3d_at_nodes_and_fractional_times(monthly_track):
    for day in (0.0, 0.123456, 7.501, 16.75, monthly_track.duration_days):
        sample = monthly_track.sample(day)
        for matrix in (sample.attitude, sample.baseline, sample.relative):
            np.testing.assert_allclose(matrix.T @ matrix, np.eye(3), atol=1e-13)
            assert np.linalg.det(matrix) == pytest.approx(1.0, abs=1e-13)
        np.testing.assert_allclose(sample.reconstructed_attitude, sample.attitude, atol=1e-13)
        assert sample.reconstruction_error < 1e-13
        assert sample.projected_rocking_rad == -sample.projected_body_relative_rad
        earth = -sample.position_km / sample.distance_km
        body_earth = sample.attitude.T @ earth
        assert math.sin(math.radians(sample.subearth_latitude_deg)) == pytest.approx(body_earth[2], abs=1e-15)


def test_fractional_state_agrees_with_direct_propagation(monthly_track):
    seconds = 7.123456 * SECONDS_PER_DAY
    index = int(np.searchsorted(monthly_track.times_seconds, seconds) - 1)
    direct = moon_model.propagate(
        monthly_track.states[index],
        np.array([monthly_track.times_seconds[index], seconds]),
        300.0,
        solar=True,
    )[-1]
    interpolated = monthly_track.sample(seconds / SECONDS_PER_DAY)
    # Hourly Hermite interpolation should preserve sub-metre lunar positions
    # and sub-millimetre/s velocities within an interval, not introduce pixels
    # of frame-rate-dependent integration drift.
    assert np.linalg.norm(interpolated.position_km - direct[:3]) < 0.001
    assert np.linalg.norm(interpolated.velocity_km_s - direct[3:6]) < 1e-6


def test_transport_rate_agrees_with_relative_attitude_derivative(monthly_track):
    day = 11.5  # An exact integration node; derivatives from either side agree.
    dt = 1.0
    sample = monthly_track.sample(day)
    before = monthly_track.sample(day - dt / SECONDS_PER_DAY)
    after = monthly_track.sample(day + dt / SECONDS_PER_DAY)
    finite_difference = (after.relative - before.relative) / (2 * dt)
    np.testing.assert_allclose(finite_difference, sample.relative_attitude_rate, atol=1e-11, rtol=0)
    np.testing.assert_allclose(
        sample.relative_rate_vector_rad_s,
        sample.baseline.T @ (
            sample.orientation_rate_rad_s * monthly_track.initial_attitude[:, 2]
            - sample.frame_rate_vector_rad_s
        ),
        atol=1e-18,
    )


def test_interval_limits_and_continuous_orbit_progress(monthly_track):
    assert monthly_track.duration_days == 28
    assert monthly_track.period_days == pytest.approx(27.321661)
    assert monthly_track.sample(0).orbital_progress_rad == pytest.approx(0)
    progress = [monthly_track.sample(day).orbital_progress_rad for day in np.linspace(0, 28, 113)]
    assert np.all(np.diff(progress) > 0)
    assert progress[-1] > 2 * math.pi
    for invalid in (-0.01, 28.01, float("nan")):
        with pytest.raises(ValueError):
            monthly_track.sample(invalid)


def test_only_first_reference_state_and_attitude_drive_model(tmp_path, monthly_track):
    import json

    path = Path(__file__).resolve().parents[1] / "data" / "moon-model-reference.json"
    reference = json.loads(path.read_text())
    # Keep the interval dates but deliberately corrupt every later input state
    # and attitude. The model must still use only first-epoch initialization.
    for key in ("moon", "sun"):
        for row in reference[key]["rows"][1:]:
            row[2:8] = ["0"] * 6
    for row in reference["orientation"]["rows"][1:]:
        row[3:7] = ["0"] * 4
    mutated = tmp_path / "modified-reference.json"
    mutated.write_text(json.dumps(reference))
    rebuilt = build_track(mutated)
    np.testing.assert_array_equal(rebuilt.states, monthly_track.states)
    np.testing.assert_array_equal(rebuilt.initial_attitude, monthly_track.initial_attitude)
