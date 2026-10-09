# SPDX-License-Identifier: MIT
"""Three-dimensional model data for a fixed-camera, top-down animation.

DREMT-linked source S1, section 4: B = [radial, transverse, orbit normal] and
relative attitude S = BᵀA. Source S2 gives the planar rate cancellation. The
matrix transport law used below is the analyst's 3D extension of those
coordinates; it supplies no physical torque or transfer of orbital energy.
S1: https://docs.google.com/document/d/1zAOk_j9BZpsHDeSfSsaQUUdbwwat_HxQ/edit
S2: https://docs.google.com/document/d/1h8mrNwdD5oIv_BD0vG98MMT5kxpp-cxi/edit
Both source snapshots were downloaded on 2026-10-05; see dremt-derived-model/README.md.

ANALYST ADDITIONS: the existing moon_model.py Newtonian Sun–Earth–Moon
propagator; empirical first-epoch position, velocity, pole, and meridian phase;
a supplied, NONZERO, uniform rate of orientation change relative to fixed
stars, with a fixed pole and no attitude torque. Later ephemeris rows do not
drive the animation. No physical libration or invented precession is added.

The display camera is fixed along the initial orbital normal. Body and orbit
are still calculated in 3D. Projected rocking is a screen-plane angle, not
literal sub-Earth longitude. Literal longitude and latent latitude are retained
for a future latitude display. All calculations are offline.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np

import moon_model
from project_paths import DATA_DIR, project_path

SECONDS_PER_DAY = 86400.0
DEFAULT_ORIENTATION_RATE = 2 * math.pi / (
    moon_model.ROTATION_PERIOD_DAYS * SECONDS_PER_DAY
)


def _wrapped(angle: float) -> float:
    return (angle + math.pi) % (2 * math.pi) - math.pi


def _heading(vector_2d: np.ndarray) -> float:
    if float(np.dot(vector_2d, vector_2d)) < 1e-24:
        raise ValueError("The projected direction is undefined at this camera angle")
    return math.atan2(float(vector_2d[1]), float(vector_2d[0]))


def _azimuth_rate(vector_2d: np.ndarray, derivative_2d: np.ndarray) -> float:
    return float(
        (vector_2d[0] * derivative_2d[1] - vector_2d[1] * derivative_2d[0])
        / np.dot(vector_2d, vector_2d)
    )


def _rotate_attitude(initial: np.ndarray, rate: float, seconds: float) -> np.ndarray:
    # The empirical initial pole is held fixed. This is an initialized attitude
    # approximation, not an attitude force law derived from either source.
    pole = initial[:, 2]
    angle = rate * seconds
    rotation = (
        math.cos(angle) * np.eye(3)
        + (1 - math.cos(angle)) * np.outer(pole, pole)
        + math.sin(angle) * moon_model.skew(pole)
    )
    return rotation @ initial


@dataclass(frozen=True, slots=True)
class AnimationSample:
    """A complete 3D state plus explicitly labelled screen-plane quantities.

    Matrices A, B, and S map body→inertial, orbital-frame→inertial, and
    body→orbital-frame coordinates, respectively. ``relative_rate_vector_rad_s``
    is expressed in B coordinates; ``frame_rate_vector_rad_s`` is inertial.
    Positive ``projected_rocking_rad`` is Earth-bearing minus marker azimuth.
    The opposite signed ``projected_body_relative_rad`` is the marker's
    orientation relative to that bearing. Neither is literal 3D longitude.
    """

    elapsed_days: float
    position_km: np.ndarray
    velocity_km_s: np.ndarray
    attitude: np.ndarray
    baseline: np.ndarray
    relative: np.ndarray
    reconstructed_attitude: np.ndarray
    relative_attitude_rate: np.ndarray
    projected_position_km: np.ndarray
    projected_velocity_km_s: np.ndarray
    projected_marker: np.ndarray
    body_marker_angle_rad: float
    earth_bearing_angle_rad: float
    projected_body_relative_rad: float
    projected_rocking_rad: float
    orbital_progress_rad: float
    orientation_angle_rad: float
    subearth_longitude_deg: float
    subearth_latitude_deg: float
    distance_km: float
    speed_km_s: float
    orientation_rate_rad_s: float
    projected_body_rate_rad_s: float
    projected_orbital_rate_rad_s: float
    projected_relative_rate_rad_s: float
    frame_rate_vector_rad_s: np.ndarray
    relative_rate_vector_rad_s: np.ndarray
    reconstruction_error: float

    @property
    def projected_rocking_deg(self) -> float:
        return math.degrees(self.projected_rocking_rad)

    @property
    def projected_body_relative_deg(self) -> float:
        return math.degrees(self.projected_body_relative_rad)


@dataclass(frozen=True, slots=True)
class AnimationTrack:
    """A propagated orbit with exact attitude evaluation at arbitrary times.

    Fractional orbital states use cubic Hermite interpolation of position and
    velocity. Attitude is evaluated analytically, and B/S are recomputed from
    the interpolated state. Rotation matrices are never linearly interpolated.
    Times outside the finite reference interval are rejected, not extrapolated.
    ``camera_basis`` columns are camera x, camera y, and viewing normal in
    inertial coordinates; project a vector with ``camera_basis[:, :2].T``.
    """

    times_seconds: np.ndarray
    states: np.ndarray
    initial_attitude: np.ndarray
    camera_basis: np.ndarray
    start_epoch: str
    end_epoch: str
    orientation_rate_rad_s: float = DEFAULT_ORIENTATION_RATE

    def __post_init__(self) -> None:
        if (
            self.times_seconds.ndim != 1
            or len(self.times_seconds) < 2
            or self.states.shape != (len(self.times_seconds), 12)
        ):
            raise ValueError("A track requires at least two times and 12-component states")
        if self.times_seconds[0] != 0 or np.any(np.diff(self.times_seconds) <= 0):
            raise ValueError("Track times must start at zero and increase strictly")
        if not np.all(np.isfinite(self.states)) or not np.all(np.isfinite(self.times_seconds)):
            raise ValueError("Track times and states must be finite")
        if not math.isfinite(self.orientation_rate_rad_s):
            raise ValueError("The orientation rate must be finite")
        for name, matrix in (("initial attitude", self.initial_attitude), ("camera basis", self.camera_basis)):
            if (
                matrix.shape != (3, 3)
                or not np.all(np.isfinite(matrix))
                or not np.allclose(matrix.T @ matrix, np.eye(3), atol=1e-12, rtol=0)
                or not math.isclose(float(np.linalg.det(matrix)), 1.0, abs_tol=1e-12)
            ):
                raise ValueError(f"The {name} must be a proper orthonormal 3D basis")

    @property
    def duration_days(self) -> float:
        return float(self.times_seconds[-1] / SECONDS_PER_DAY)

    @property
    def period_days(self) -> float:
        return 2 * math.pi / self.orientation_rate_rad_s / SECONDS_PER_DAY

    @property
    def times_days(self) -> np.ndarray:
        return self.times_seconds / SECONDS_PER_DAY

    @property
    def projected_orbit_km(self) -> np.ndarray:
        return self.states[:, :3] @ self.camera_basis[:, :2]

    @property
    def projected_orbital_progress_rad(self) -> np.ndarray:
        projected = self.projected_orbit_km
        angles = np.unwrap(np.arctan2(projected[:, 1], projected[:, 0]))
        return angles - angles[0]

    def _state_at_seconds(self, seconds: float) -> np.ndarray:
        if not math.isfinite(seconds) or seconds < 0 or seconds > self.times_seconds[-1]:
            raise ValueError(f"Elapsed time must be between 0 and {self.duration_days:g} days")
        if seconds == self.times_seconds[-1]:
            return self.states[-1].copy()
        index = int(np.searchsorted(self.times_seconds, seconds, side="right") - 1)
        t0, t1 = self.times_seconds[index:index + 2]
        left, right = self.states[index:index + 2]
        span = float(t1 - t0)
        u = (seconds - t0) / span
        state = np.empty(12)
        for offset in (0, 6):
            p0, v0 = left[offset:offset + 3], left[offset + 3:offset + 6]
            p1, v1 = right[offset:offset + 3], right[offset + 3:offset + 6]
            state[offset:offset + 3] = (
                (2 * u**3 - 3 * u**2 + 1) * p0
                + (u**3 - 2 * u**2 + u) * span * v0
                + (-2 * u**3 + 3 * u**2) * p1
                + (u**3 - u**2) * span * v1
            )
            state[offset + 3:offset + 6] = (
                (6 * u**2 - 6 * u) * p0 / span
                + (3 * u**2 - 4 * u + 1) * v0
                + (-6 * u**2 + 6 * u) * p1 / span
                + (3 * u**2 - 2 * u) * v1
            )
        return state

    def sample(self, elapsed_days: float) -> AnimationSample:
        seconds = float(elapsed_days) * SECONDS_PER_DAY
        state = self._state_at_seconds(seconds)
        attitude = _rotate_attitude(self.initial_attitude, self.orientation_rate_rad_s, seconds)
        # S1's quaternion composition, expressed as rotation matrices. In 3D
        # the relative correction is a full rotation, not one additive angle.
        geometric_rate = moon_model.acceleration(state, solar=True)
        baseline, frame_rate = moon_model.moving_frame(state, geometric_rate)
        relative = baseline.T @ attitude
        recovered = baseline @ relative
        omega = self.orientation_rate_rad_s * self.initial_attitude[:, 2]
        # Analyst's transport extension: Ṡ = [Bᵀ(ω−Ω_frame)]× S. In the
        # fixed-axis planar limit, ω_relative = C−Ω. This is coordinate motion.
        relative_rate = baseline.T @ (omega - frame_rate)
        matrix_rate = moon_model.skew(relative_rate) @ relative

        position, velocity = state[:3].copy(), state[3:6].copy()
        projection = self.camera_basis[:, :2].T
        projected_position = projection @ position
        projected_velocity = projection @ velocity
        projected_marker = projection @ attitude[:, 0]
        projected_marker_rate = projection @ np.cross(omega, attitude[:, 0])
        marker_angle = _heading(projected_marker)
        earth_bearing = _heading(-projected_position)
        body_relative = _wrapped(marker_angle - earth_bearing)
        orbital_rate = _azimuth_rate(projected_position, projected_velocity)
        body_rate = _azimuth_rate(projected_marker, projected_marker_rate)
        body_earth = attitude.T @ (-position / np.linalg.norm(position))

        # Unwrap the current orbital azimuth using a neighbouring stored node.
        # This reports progress without resetting the counter at ±π.
        stored_progress = self.projected_orbital_progress_rad
        node = min(int(np.searchsorted(self.times_seconds, seconds)), len(stored_progress) - 1)
        origin_angle = _heading(self.projected_orbit_km[0])
        raw_progress = _heading(projected_position) - origin_angle
        progress = raw_progress + 2 * math.pi * round((stored_progress[node] - raw_progress) / (2 * math.pi))
        return AnimationSample(
            elapsed_days=float(elapsed_days),
            position_km=position,
            velocity_km_s=velocity,
            attitude=attitude,
            baseline=baseline,
            relative=relative,
            reconstructed_attitude=recovered,
            relative_attitude_rate=matrix_rate,
            projected_position_km=projected_position,
            projected_velocity_km_s=projected_velocity,
            projected_marker=projected_marker,
            body_marker_angle_rad=marker_angle,
            earth_bearing_angle_rad=earth_bearing,
            projected_body_relative_rad=body_relative,
            projected_rocking_rad=-body_relative,
            orbital_progress_rad=progress,
            orientation_angle_rad=self.orientation_rate_rad_s * seconds,
            subearth_longitude_deg=math.degrees(math.atan2(float(body_earth[1]), float(body_earth[0]))),
            subearth_latitude_deg=math.degrees(math.asin(float(np.clip(body_earth[2], -1, 1)))),
            distance_km=float(np.linalg.norm(position)),
            speed_km_s=float(np.linalg.norm(velocity)),
            orientation_rate_rad_s=self.orientation_rate_rad_s,
            projected_body_rate_rad_s=body_rate,
            projected_orbital_rate_rad_s=orbital_rate,
            projected_relative_rate_rad_s=body_rate - orbital_rate,
            frame_rate_vector_rad_s=frame_rate,
            relative_rate_vector_rad_s=relative_rate,
            reconstruction_error=float(np.max(np.abs(recovered - attitude))),
        )


def build_track(
    reference_path: Path | str | None = None,
    sample_seconds: float = 3600.0,
    max_step_seconds: float = 300.0,
) -> AnimationTrack:
    """Propagate the cached 28-day interval, using only the first state/pole.

    Remaining reference rows determine the interval length, not model states.
    Node spacing controls display interpolation; RK4 step controls integration.
    """
    if not math.isfinite(sample_seconds) or sample_seconds <= 0:
        raise ValueError("Sample spacing must be a positive finite number")
    if not math.isfinite(max_step_seconds) or max_step_seconds <= 0:
        raise ValueError("Integration step must be a positive finite number")
    path = project_path(reference_path) if reference_path is not None else (DATA_DIR / "moon-model-reference.json")
    reference = json.loads(path.read_text())
    moon_rows = reference["moon"]["rows"]
    duration = (float(moon_rows[-1][0]) - float(moon_rows[0][0])) * SECONDS_PER_DAY
    if duration <= 0:
        raise ValueError("The reference interval must have positive duration")
    times = np.arange(0.0, duration, sample_seconds)
    times = np.append(times, duration)
    initial = np.array(
        [float(value) for value in moon_rows[0][2:8]]
        + [float(value) for value in reference["sun"]["rows"][0][2:8]]
    )
    states = moon_model.propagate(initial, times, max_step_seconds, solar=True)
    initial_orientation = np.array([float(value) for value in reference["orientation"]["rows"][0][3:7]])
    initial_attitude = moon_model.apparent_orientation(states[:1], np.array([0.0]), initial_orientation)[4][0]
    camera_basis = moon_model.moving_frame(initial)
    return AnimationTrack(
        times_seconds=times,
        states=states,
        initial_attitude=initial_attitude,
        camera_basis=camera_basis,
        start_epoch=moon_rows[0][1].strip(),
        end_epoch=moon_rows[-1][1].strip(),
    )
