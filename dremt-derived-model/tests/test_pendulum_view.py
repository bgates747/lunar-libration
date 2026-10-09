"""Check that the middle view changes coordinates without changing motion."""
import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from orbit_comparison import build_perigee_comparison_track, SECONDS_PER_DAY
from pendulum_view import pendulum_view


@pytest.fixture(scope="module")
def track():
    return build_perigee_comparison_track()


def test_fixed_virtual_moon_vertical_tether_and_preserved_geometry(track):
    for day in (0, 3.1, 6, 17.3, track.duration_days):
        sample = track.sample(day)
        view = pendulum_view(sample)
        np.testing.assert_array_equal(view.virtual_moon_km, np.zeros(2))
        assert view.earth_km[0] == pytest.approx(0, abs=1e-9)
        assert view.earth_km[1] > 0
        real_radial = view.real_moon_km-view.earth_km
        virtual_radial = view.virtual_moon_km-view.earth_km
        np.testing.assert_allclose(view.rotation.T @ real_radial, sample.projected_position_km, atol=1e-9)
        np.testing.assert_allclose(view.rotation.T @ virtual_radial, sample.projected_virtual_position_km, atol=1e-9)
        assert np.linalg.norm(view.real_moon_km) == pytest.approx(np.linalg.norm(sample.projected_position_km-sample.projected_virtual_position_km), abs=1e-9)
        angle = math.atan2(real_radial[0]*virtual_radial[1]-real_radial[1]*virtual_radial[0], np.dot(real_radial,virtual_radial))
        assert angle == pytest.approx(sample.angular_separation_rad, abs=1e-14)


def test_rotating_view_velocity_and_inertial_marker_agree_with_finite_differences(track):
    for day in (0.5, 6.123, 17.41):
        sample = track.sample(day)
        view = pendulum_view(sample)
        before = pendulum_view(track.sample(day-1/SECONDS_PER_DAY))
        after = pendulum_view(track.sample(day+1/SECONDS_PER_DAY))
        finite_velocity = (after.real_moon_km-before.real_moon_km)/2
        np.testing.assert_allclose(view.real_velocity_km_s, finite_velocity, atol=1e-8, rtol=0)
        angle_rate = (after.earth_graphic_angle_rad-before.earth_graphic_angle_rad)/2
        assert angle_rate == pytest.approx(-sample.mean_angular_rate_rad_s, abs=1e-14)
        # A fixed x-axis marker in the main view must be recovered by the inverse
        # transform: the Earth's displayed turning comes from camera rotation.
        marker = np.array([math.cos(view.earth_graphic_angle_rad), math.sin(view.earth_graphic_angle_rad)])
        np.testing.assert_allclose(view.rotation.T @ marker, [1,0], atol=1e-14)
    assert pendulum_view(track.sample(0)).real_velocity_km_s[0] > 0
