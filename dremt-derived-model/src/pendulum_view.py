# SPDX-License-Identifier: MIT
"""A rotating, panning view of the same projected ephemeris comparison.

Analyst visualization: rotate by β = −π/2 − α_virtual, then translate the
virtual Moon to a fixed screen anchor. Its Earthward line is exactly vertical.
Translation makes Earth move vertically as the shared radial distance changes;
no length is rescaled. Earth's graphic turns by β to mark fixed inertial axes,
not Earth's physical daily rotation. This adds no pendulum force or constraint.
"""
from dataclasses import dataclass
import math

import numpy as np

from orbit_comparison import ComparisonSample


@dataclass(frozen=True, slots=True)
class PendulumView:
    # Cartesian coordinates in km, relative to the stationary virtual Moon.
    earth_km: np.ndarray
    real_moon_km: np.ndarray
    virtual_moon_km: np.ndarray
    real_velocity_km_s: np.ndarray
    rotation: np.ndarray
    earth_graphic_angle_rad: float
    real_marker_angle_rad: float


def pendulum_view(sample: ComparisonSample) -> PendulumView:
    beta = -math.pi/2 - sample.virtual_angle_rad
    c, s = math.cos(beta), math.sin(beta)
    rotation = np.array([[c, -s], [s, c]])
    real = sample.projected_position_km
    virtual = sample.projected_virtual_position_km
    relative = real - virtual
    radius = float(np.linalg.norm(real))
    radial_rate = float(np.dot(real, sample.projected_velocity_km_s))/radius
    turn = np.array([[0., -1.], [1., 0.]])
    virtual_velocity = radial_rate*virtual/radius + sample.mean_angular_rate_rad_s*(turn @ virtual)
    relative_velocity = (sample.projected_velocity_km_s - virtual_velocity
                         - sample.mean_angular_rate_rad_s*(turn @ relative))
    return PendulumView(
        earth_km=-(rotation @ virtual), real_moon_km=rotation @ relative,
        virtual_moon_km=np.zeros(2), real_velocity_km_s=rotation @ relative_velocity,
        rotation=rotation, earth_graphic_angle_rad=beta,
        real_marker_angle_rad=sample.body_marker_angle_rad + beta,
    )
