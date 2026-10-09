"""Verify the new analytic track and its relation to the original animation."""
import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from keplerian_comparison import KeplerianComparisonTrack
from pendulum_view import pendulum_view


def test_original_relative_angles_apses_and_symmetry():
    track = KeplerianComparisonTrack()
    T = track.duration_days
    for day in (0, T/7, T/4, T/2, 3*T/4, T):
        sample = track.sample(day)
        original = track.orbit.state(day)
        assert sample.angular_separation_rad == pytest.approx(-original.libration_rad, abs=1e-14)
        assert sample.speed_km_s == pytest.approx(original.speed_km_s, abs=1e-14)
        opposite = track.sample(T-day)
        assert sample.angular_separation_rad == pytest.approx(-opposite.angular_separation_rad, abs=1e-14)
    for event in track.apses:
        assert event["angular_separation_deg"] == pytest.approx(0, abs=1e-12)
    peak = max(abs(track.sample(float(t)).angular_separation_rad)
               for t in np.linspace(0, T, 1001))
    assert math.degrees(peak) == pytest.approx(6.29325, abs=1e-4)


def test_velocity_kepler_invariants_and_pendulum_transform():
    track = KeplerianComparisonTrack()
    a, n = track.orbit.semi_major_axis_km, track.mean_angular_rate_rad_s
    mu = n*n*a**3  # Original model's supplied a/T determine its effective GM.
    angular_momentum = a*a*n*math.sqrt(1-track.orbit.eccentricity**2)
    for day in (0.7, 6.1, 14.2, 23.4):
        s = track.sample(day)
        before, after = track.sample(day-1/86400), track.sample(day+1/86400)
        np.testing.assert_allclose((after.position_km-before.position_km)/2, s.velocity_km_s, atol=1e-9, rtol=0)
        assert np.linalg.norm(np.cross(s.position_km, s.velocity_km_s)) == pytest.approx(angular_momentum, rel=1e-14)
        assert s.speed_km_s**2/2-mu/s.distance_km == pytest.approx(-mu/(2*a), abs=1e-14)
        view = pendulum_view(s)
        assert view.earth_km[0] == pytest.approx(0, abs=1e-9)
        np.testing.assert_allclose(view.rotation.T@(view.real_moon_km-view.earth_km), s.projected_position_km, atol=1e-9)
        assert np.linalg.norm(view.real_moon_km-view.earth_km) == pytest.approx(s.distance_km)
