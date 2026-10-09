"""Keplerian physics for the simplified top-down Moon model.

The orbit is an ideal two-body ellipse.  The lunar spin is uniform and has
the same period as the orbit.  Inclination, axial tilt, perturbations, and
physical libration are deliberately omitted so the optical longitude effect
is isolated.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

TAU = 2.0 * math.pi


@dataclass(frozen=True)
class OrbitalState:
    elapsed_days: float
    orbit_fraction: float
    mean_anomaly: float
    eccentric_anomaly: float
    true_anomaly: float
    x_km: float
    y_km: float
    distance_km: float
    speed_km_s: float
    spin_angle: float
    earth_bearing: float
    libration_rad: float

    @property
    def libration_deg(self) -> float:
        return math.degrees(self.libration_rad)


@dataclass(frozen=True)
class LunarOrbit:
    semi_major_axis_km: float = 384_400.0
    eccentricity: float = 0.0549
    period_days: float = 27.321661

    @staticmethod
    def _wrap(angle: float) -> float:
        return math.atan2(math.sin(angle), math.cos(angle))

    def eccentric_anomaly(self, mean_anomaly: float) -> float:
        """Solve M = E - e sin(E) with Newton iteration."""
        turns = math.floor(mean_anomaly / TAU)
        reduced_m = mean_anomaly - turns * TAU
        estimate = reduced_m + self.eccentricity * math.sin(reduced_m)
        for _ in range(12):
            residual = estimate - self.eccentricity * math.sin(estimate) - reduced_m
            derivative = 1.0 - self.eccentricity * math.cos(estimate)
            correction = residual / derivative
            estimate -= correction
            if abs(correction) < 1e-14:
                break
        return estimate + turns * TAU

    def state(self, elapsed_days: float, start_apsis: str = "perigee") -> OrbitalState:
        if start_apsis not in {"perigee", "apogee"}:
            raise ValueError("start_apsis must be 'perigee' or 'apogee'")

        fraction = elapsed_days / self.period_days
        progress = TAU * fraction
        initial_mean = 0.0 if start_apsis == "perigee" else math.pi
        mean = initial_mean + progress
        eccentric = self.eccentric_anomaly(mean)

        a = self.semi_major_axis_km
        e = self.eccentricity
        x = a * (math.cos(eccentric) - e)
        y = a * math.sqrt(1.0 - e * e) * math.sin(eccentric)
        distance = math.hypot(x, y)

        # Kepler's third law gives mu for this idealized two-body orbit; the
        # vis-viva equation then gives the instantaneous orbital speed.
        period_seconds = self.period_days * 86_400.0
        mu = 4.0 * math.pi**2 * a**3 / period_seconds**2
        speed = math.sqrt(mu * (2.0 / distance - 1.0 / a))

        true_wrapped = math.atan2(y, x)
        true_unwrapped = true_wrapped + TAU * round((mean - true_wrapped) / TAU)
        earth_bearing = true_unwrapped + math.pi
        initial_body_angle = math.pi if start_apsis == "perigee" else 0.0
        spin_angle = initial_body_angle + progress
        libration = self._wrap(earth_bearing - spin_angle)

        return OrbitalState(
            elapsed_days=elapsed_days,
            orbit_fraction=fraction,
            mean_anomaly=mean,
            eccentric_anomaly=eccentric,
            true_anomaly=true_unwrapped,
            x_km=x,
            y_km=y,
            distance_km=distance,
            speed_km_s=speed,
            spin_angle=spin_angle,
            earth_bearing=earth_bearing,
            libration_rad=libration,
        )
