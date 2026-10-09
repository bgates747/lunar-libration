# SPDX-License-Identifier: MIT
"""Ideal-ellipse comparison using the original author-2d-animation orbit.

Analyst-supplied assumptions: a = 384400 km, e = 0.0549, T = 27.321661 days;
uniform orientation change C = 2π/T, initially facing Earth at perigee.
The virtual Moon shares the actual model Moon's radius and advances at C.
There is no ephemeris input, physical pendulum force, or attitude torque.

DREMT S2's rate identity integrates to ω_s = C − Ω; thus the relative angle
is M − f for these initial conditions. S2 does not determine this orbit or C.
Source: https://docs.google.com/document/d/1h8mrNwdD5oIv_BD0vG98MMT5kxpp-cxi/edit
The local 2026-10-05 snapshot and provenance are documented in dremt-derived-model/README.md.
"""
from dataclasses import dataclass, field
import math
from project_paths import PROJECT_ROOT
import sys

import numpy as np

from orbit_comparison import ComparisonSample

# Reuse the original calculation, rather than fit an ellipse to the ephemeris.
sys.path.insert(0, str(PROJECT_ROOT / "author-2d-animation" / "src"))
from moon_libration.orbit import LunarOrbit


@dataclass(frozen=True)
class KeplerianComparisonTrack:
    orbit: LunarOrbit = field(default_factory=LunarOrbit)
    virtual_phase_rad: float = 0.0

    @property
    def duration_days(self):
        return self.orbit.period_days

    @property
    def mean_angular_rate_rad_s(self):
        return 2 * math.pi / (self.duration_days * 86400)

    @property
    def projected_orbit_km(self):
        return np.array([self.sample(float(t)).projected_position_km
                         for t in np.linspace(0, self.duration_days, 721)])

    @property
    def apses(self):
        result = []
        for kind, day in (("perigee", 0.), ("apogee", self.duration_days/2),
                          ("perigee", self.duration_days)):
            s = self.sample(day)
            result.append({"kind": kind, "elapsed_days": day,
                           "distance_km": s.distance_km,
                           "geometric_angle_deg": math.degrees(s.geometric_angle_rad),
                           "angular_separation_deg": math.degrees(s.angular_separation_rad)})
        return tuple(result)

    def sample(self, elapsed_days):
        if not math.isfinite(elapsed_days) or not 0 <= elapsed_days <= self.duration_days:
            raise ValueError("Elapsed time must lie within the ideal orbit")
        state = self.orbit.state(elapsed_days)
        a, e, E = self.orbit.semi_major_axis_km, self.orbit.eccentricity, state.eccentric_anomaly
        n = self.mean_angular_rate_rad_s
        # Exact time derivative of x=a(cos E−e), y=a√(1−e²) sin E.
        E_rate = n / (1-e*math.cos(E))
        projected = np.array([state.x_km, state.y_km])
        velocity = a*E_rate*np.array([-math.sin(E), math.sqrt(1-e*e)*math.cos(E)])
        f, M = state.true_anomaly, state.mean_anomaly
        virtual_angle = M + self.virtual_phase_rad
        virtual = state.distance_km*np.array([math.cos(virtual_angle), math.sin(virtual_angle)])
        theta = state.spin_angle  # Original model: θ = π + M, a supplied orientation law.
        c, s = math.cos(theta), math.sin(theta)
        attitude = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
        return ComparisonSample(
            elapsed_days=float(elapsed_days), position_km=np.r_[projected, 0.],
            velocity_km_s=np.r_[velocity, 0.], virtual_position_km=np.r_[virtual, 0.],
            attitude=attitude, projected_position_km=projected,
            projected_virtual_position_km=virtual, projected_velocity_km_s=velocity,
            body_marker_angle_rad=theta, real_angle_rad=f, virtual_angle_rad=virtual_angle,
            geometric_angle_rad=f,
            angular_separation_rad=math.atan2(math.sin(virtual_angle-f), math.cos(virtual_angle-f)),
            distance_km=state.distance_km, speed_km_s=float(np.linalg.norm(velocity)),
            mean_angular_rate_rad_s=n,
            real_angular_rate_rad_s=float(np.cross(np.r_[projected, 0.], np.r_[velocity, 0.])[2])/state.distance_km**2,
            subearth_longitude_deg=state.libration_deg, subearth_latitude_deg=0.,
        )
