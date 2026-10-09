import math

import pytest

from moon_libration.orbit import LunarOrbit, TAU


@pytest.fixture
def orbit():
    return LunarOrbit()


def test_perigee_and_apogee_distances(orbit):
    perigee = orbit.state(0, "perigee")
    apogee = orbit.state(0, "apogee")
    assert perigee.distance_km == pytest.approx(orbit.semi_major_axis_km * (1 - orbit.eccentricity))
    assert apogee.distance_km == pytest.approx(orbit.semi_major_axis_km * (1 + orbit.eccentricity))


def test_moon_moves_faster_at_perigee(orbit):
    assert orbit.state(0, "perigee").speed_km_s > orbit.state(0, "apogee").speed_km_s


@pytest.mark.parametrize("start", ["perigee", "apogee"])
def test_reference_marker_is_aligned_at_start_and_after_one_orbit(orbit, start):
    assert orbit.state(0, start).libration_deg == pytest.approx(0, abs=1e-10)
    assert orbit.state(orbit.period_days, start).libration_deg == pytest.approx(0, abs=1e-10)


def test_kepler_equation_residual(orbit):
    for mean in [0.0, 0.1, 1.2, math.pi, 5.9, TAU, 2 * TAU + 0.4]:
        eccentric = orbit.eccentric_anomaly(mean)
        assert eccentric - orbit.eccentricity * math.sin(eccentric) == pytest.approx(mean, abs=1e-12)


def test_expected_libration_amplitude(orbit):
    values = [
        abs(orbit.state(orbit.period_days * i / 2000, "perigee").libration_deg)
        for i in range(2001)
    ]
    assert max(values) == pytest.approx(6.30, abs=0.05)
