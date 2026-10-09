# SPDX-License-Identifier: MIT
"""Ephemeris Moon and a uniformly advancing virtual Moon, viewed from Earth.

DREMT's S1/S2 define a baseline along the body-to-host-COM line, but do not
specify this virtual Moon or a pendulum suspension mechanism. This comparison
is an analyst visualization requested by the repository author.
S1: https://docs.google.com/document/d/1zAOk_j9BZpsHDeSfSsaQUUdbwwat_HxQ/edit
S2: https://docs.google.com/document/d/1h8mrNwdD5oIv_BD0vG98MMT5kxpp-cxi/edit

REAL MOON: all cached NASA/JPL Horizons position/velocity and orientation rows
drive this display. It is an ephemeris playback, not a model prediction.
VIRTUAL MOON ASSUMPTIONS: constant projected angular rate 2π/27.321661 days;
initial bearing equal to the real Moon (plus an optional phase offset); same
instantaneous projected radius and out-of-plane height as the real Moon.
Thus its 3D Earth distance also matches. Its tether length varies; no force,
rigid arm, physical pivot, or virtual body orientation is inferred.

The perigee-start camera remains fixed along the initial orbital normal, with
its +x direction pointing to the first perigee (geometric angle 0°). Real 3D states,
body attitude and sub-Earth latitude are retained for a later latitude view.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime, timedelta
import json
import math
from pathlib import Path

import numpy as np

from animation_model import AnimationTrack, SECONDS_PER_DAY
import moon_model
from project_paths import DATA_DIR, project_path

# Analyst-supplied mean sidereal orbital period, not a fitted monthly rate or
# the body's instantaneous orientation rate about its tilted pole.
MEAN_SIDEREAL_PERIOD_DAYS = 27.321661
DEFAULT_MEAN_MOTION = 2 * math.pi / (MEAN_SIDEREAL_PERIOD_DAYS * SECONDS_PER_DAY)


@dataclass(frozen=True, slots=True)
class ComparisonSample:
    elapsed_days: float
    position_km: np.ndarray
    velocity_km_s: np.ndarray
    virtual_position_km: np.ndarray
    attitude: np.ndarray
    projected_position_km: np.ndarray
    projected_virtual_position_km: np.ndarray
    projected_velocity_km_s: np.ndarray
    body_marker_angle_rad: float
    real_angle_rad: float
    virtual_angle_rad: float
    geometric_angle_rad: float
    angular_separation_rad: float
    distance_km: float
    speed_km_s: float
    mean_angular_rate_rad_s: float
    real_angular_rate_rad_s: float
    subearth_longitude_deg: float
    subearth_latitude_deg: float


@dataclass(frozen=True, slots=True)
class EphemerisComparisonTrack:
    orbit: AnimationTrack
    orientation_rows: np.ndarray
    mean_angular_rate_rad_s: float = DEFAULT_MEAN_MOTION
    virtual_phase_rad: float = 0.0
    start_seconds: float = 0.0
    stop_seconds: float | None = None
    reference_start_jd: float | None = None
    apses: tuple[dict, ...] = ()

    def __post_init__(self):
        if self.orientation_rows.shape != (len(self.orbit.times_seconds), 4):
            raise ValueError("One longitude/latitude/pole row is required per ephemeris epoch")
        if not np.all(np.isfinite(self.orientation_rows)):
            raise ValueError("Orientation rows must be finite")
        if not math.isfinite(self.mean_angular_rate_rad_s) or self.mean_angular_rate_rad_s <= 0:
            raise ValueError("The mean angular rate must be positive and finite")
        if not math.isfinite(self.virtual_phase_rad):
            raise ValueError("The virtual phase must be finite")
        stop = self.orbit.times_seconds[-1] if self.stop_seconds is None else self.stop_seconds
        if not 0 <= self.start_seconds < stop <= self.orbit.times_seconds[-1]:
            raise ValueError("Playback interval must lie within the cached ephemeris")

    @property
    def duration_days(self):
        stop = self.orbit.times_seconds[-1] if self.stop_seconds is None else self.stop_seconds
        return float((stop - self.start_seconds) / SECONDS_PER_DAY)

    @property
    def projected_orbit_km(self):
        stop = self.orbit.times_seconds[-1] if self.stop_seconds is None else self.stop_seconds
        mask = (self.orbit.times_seconds > self.start_seconds) & (self.orbit.times_seconds < stop)
        positions = np.vstack((self.orbit._state_at_seconds(self.start_seconds)[:3],
                               self.orbit.states[mask, :3], self.orbit._state_at_seconds(stop)[:3]))
        return positions @ self.orbit.camera_basis[:, :2]

    @property
    def start_epoch(self):
        if self.reference_start_jd is None:
            return self.orbit.start_epoch
        return _calendar_epoch(self.reference_start_jd + self.start_seconds / SECONDS_PER_DAY)

    @property
    def end_epoch(self):
        if self.reference_start_jd is None:
            return self.orbit.end_epoch
        return _calendar_epoch(self.reference_start_jd + self.start_seconds / SECONDS_PER_DAY + self.duration_days)

    def sample(self, elapsed_days):
        if not math.isfinite(elapsed_days) or not 0 <= elapsed_days <= self.duration_days:
            raise ValueError(f"Elapsed time must be between 0 and {self.duration_days:g} days")
        elapsed_seconds = float(elapsed_days) * SECONDS_PER_DAY
        seconds = self.start_seconds + elapsed_seconds
        # Reuse the validated cubic Hermite position/velocity interpolation;
        # the stored nodes here are Horizons outputs, not integrated states.
        state = self.orbit._state_at_seconds(seconds)
        position, velocity = state[:3].copy(), state[3:6].copy()
        camera = self.orbit.camera_basis
        projected = camera[:, :2].T @ position
        projected_velocity = camera[:, :2].T @ velocity
        radius = float(np.linalg.norm(projected))
        initial = self.projected_orbit_km[0]
        origin = math.atan2(initial[1], initial[0])
        real_angle = math.atan2(projected[1], projected[0])
        virtual_angle = origin + self.virtual_phase_rad + self.mean_angular_rate_rad_s * elapsed_seconds
        virtual_projected = radius * np.array([math.cos(virtual_angle), math.sin(virtual_angle)])
        # Retain the hidden height as well as the screen-plane radius. The
        # displayed comparison changes azimuth alone, not radial distance.
        virtual_position = camera[:, :2] @ virtual_projected + camera[:, 2] * np.dot(camera[:, 2], position)
        separation = (virtual_angle - real_angle + math.pi) % (2 * math.pi) - math.pi
        real_rate = float(projected[0]*projected_velocity[1] - projected[1]*projected_velocity[0]) / radius**2
        # Report continuous geometric progress in fixed axes, not a longitude
        # reset at each later apsis. Actual apsidal precession is retained.
        full_projection = self.orbit.projected_orbit_km
        unwrapped = np.unwrap(np.arctan2(full_projection[:, 1], full_projection[:, 0]))
        node = min(int(np.searchsorted(self.orbit.times_seconds, seconds)), len(unwrapped)-1)
        real_unwrapped = real_angle + 2*math.pi*round((unwrapped[node]-real_angle)/(2*math.pi))
        first_node = min(int(np.searchsorted(self.orbit.times_seconds, self.start_seconds)), len(unwrapped)-1)
        origin_unwrapped = origin + 2*math.pi*round((unwrapped[first_node]-origin)/(2*math.pi))
        geometric_angle = real_unwrapped - origin_unwrapped

        # Longitude and pole RA are circular data: unwrap before interpolation.
        values = []
        for column in range(4):
            data = self.orientation_rows[:, column]
            if column in (0, 2):
                data = np.degrees(np.unwrap(np.radians(data)))
            values.append(float(np.interp(seconds, self.orbit.times_seconds, data)))
        longitude, latitude, ra, dec = values
        ra, dec, longitude_radians = map(math.radians, (ra, dec, longitude))
        north = np.array([math.cos(dec)*math.cos(ra), math.cos(dec)*math.sin(ra), math.sin(dec)])
        earth = -position / np.linalg.norm(position)
        equatorial_earth = earth - np.dot(earth, north)*north
        equatorial_earth /= np.linalg.norm(equatorial_earth)
        x = math.cos(longitude_radians)*equatorial_earth - math.sin(longitude_radians)*np.cross(north, equatorial_earth)
        attitude = np.column_stack((x, np.cross(north, x), north))
        marker = camera[:, :2].T @ x
        return ComparisonSample(
            elapsed_days=float(elapsed_days), position_km=position, velocity_km_s=velocity,
            virtual_position_km=virtual_position, attitude=attitude,
            projected_position_km=projected, projected_virtual_position_km=virtual_projected,
            projected_velocity_km_s=projected_velocity,
            body_marker_angle_rad=math.atan2(marker[1], marker[0]),
            real_angle_rad=real_angle, virtual_angle_rad=virtual_angle,
            geometric_angle_rad=geometric_angle,
            angular_separation_rad=separation,
            distance_km=float(np.linalg.norm(position)), speed_km_s=float(np.linalg.norm(velocity)),
            mean_angular_rate_rad_s=self.mean_angular_rate_rad_s, real_angular_rate_rad_s=real_rate,
            subearth_longitude_deg=(longitude + 180) % 360 - 180, subearth_latitude_deg=latitude,
        )


def build_comparison_track(reference_path=None, virtual_phase_rad=0.0):
    path = project_path(reference_path) if reference_path else (DATA_DIR / "moon-model-reference.json")
    reference = json.loads(path.read_text())
    moon = reference["moon"]["rows"]
    sun = reference["sun"]["rows"]
    orientation = reference["orientation"]["rows"]
    if len(moon) != len(sun) or len(moon) != len(orientation):
        raise ValueError("Reference state and orientation grids must have equal lengths")
    if any(float(m[0]) != float(s[0]) for m, s in zip(moon, sun)):
        raise ValueError("Moon and Sun vector epochs must match")
    times = np.array([(float(row[0])-float(moon[0][0]))*SECONDS_PER_DAY for row in moon])
    states = np.array([[float(v) for v in m[2:8]+s[2:8]] for m, s in zip(moon, sun)])
    orientation_rows = np.array([[float(v) for v in row[3:7]] for row in orientation])
    initial_attitude = moon_model.apparent_orientation(states[:1], np.array([0.0]), orientation_rows[0])[4][0]
    orbit = AnimationTrack(times, states, initial_attitude, moon_model.moving_frame(states[0]),
                           moon[0][1].strip(), moon[-1][1].strip())
    return EphemerisComparisonTrack(orbit, orientation_rows, virtual_phase_rad=virtual_phase_rad)


def _calendar_epoch(julian_day):
    # JD is TDB here. Format its calendar representation without falsely
    # labelling this ephemeris time as UTC or local civil time.
    epoch = datetime(2000, 1, 1, 12) + timedelta(days=julian_day-2451545.0)
    return epoch.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3] + " TDB"


def _find_apses(orbit):
    """Find radial-velocity sign changes and refine them within Hermite spans."""
    def radial_velocity(seconds):
        state = orbit._state_at_seconds(seconds)
        return float(np.dot(state[:3], state[3:6]) / np.linalg.norm(state[:3]))

    values = [radial_velocity(float(time)) for time in orbit.times_seconds]
    roots = []
    for index in range(len(values)-1):
        left, right = values[index:index+2]
        if left == 0:
            if index > 0 and values[index-1]*right < 0:
                roots.append(("perigee" if values[index-1] < 0 else "apogee", float(orbit.times_seconds[index])))
            continue
        if left*right > 0:
            continue
        kind = "perigee" if left < 0 else "apogee"
        lo, hi = map(float, orbit.times_seconds[index:index+2])
        for _ in range(45):
            mid = (lo+hi)/2
            value = radial_velocity(mid)
            if left*value <= 0:
                hi = mid
            else:
                lo, left = mid, value
        root = (lo+hi)/2
        if not roots or abs(root-roots[-1][1]) > 1e-6:
            roots.append((kind, root))
    return roots


def _first_revolution_end(orbit, start):
    """First complete projected revolution in fixed axes after a given time."""
    projected = orbit.projected_orbit_km
    angles = np.unwrap(np.arctan2(projected[:, 1], projected[:, 0]))
    start_position = orbit.camera_basis[:, :2].T @ orbit._state_at_seconds(start)[:3]
    start_angle = math.atan2(start_position[1], start_position[0])
    node = min(int(np.searchsorted(orbit.times_seconds, start)), len(angles)-1)
    start_unwrapped = start_angle + 2*math.pi*round((angles[node]-start_angle)/(2*math.pi))
    target = start_unwrapped + 2*math.pi
    crossings = np.flatnonzero((orbit.times_seconds > start) & (angles >= target))
    if not len(crossings):
        return None
    index = int(crossings[0])
    lo, hi = map(float, orbit.times_seconds[index-1:index+1])
    for _ in range(45):
        mid = (lo+hi)/2
        position = orbit.camera_basis[:, :2].T @ orbit._state_at_seconds(mid)[:3]
        angle = math.atan2(position[1], position[0])
        unwrapped = angle + 2*math.pi*round((target-angle)/(2*math.pi))
        if unwrapped < target:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def build_perigee_comparison_track(reference_path=None, virtual_phase_rad=0.0):
    """Start at the first perigee with a full geometric revolution available.

    The camera's fixed +x axis is the first perigee direction. This choice of
    coordinate origin and playback interval does not fit the uniform virtual
    rate, force symmetry, or reset the later perigee to 0°. The end is the
    first return to the initial bearing after 360°, not a second radial minimum.
    """
    path = project_path(reference_path) if reference_path else (DATA_DIR / "animation-ephemeris.json")
    raw = build_comparison_track(path, virtual_phase_rad)
    roots = _find_apses(raw.orbit)
    selected = None
    for kind, start in roots:
        if kind != "perigee":
            continue
        camera = moon_model.moving_frame(raw.orbit._state_at_seconds(start))
        orbit = replace(raw.orbit, camera_basis=camera)
        stop = _first_revolution_end(orbit, start)
        if stop is not None:
            selected = (orbit, start, stop)
            break
    if selected is None:
        raise ValueError("Cache needs a full revolution after a perigee; run python dremt-derived-model/scripts/fetch_animation_ephemeris.py")
    orbit, start, stop = selected
    reference = json.loads(path.read_text())
    jd = float(reference["moon"]["rows"][0][0])
    track = replace(raw, orbit=orbit, start_seconds=start, stop_seconds=stop, reference_start_jd=jd)
    apses = []
    for kind, time in roots:
        if start <= time <= stop:
            elapsed = (time-start)/SECONDS_PER_DAY
            sample = track.sample(elapsed)
            apses.append({"kind": kind, "elapsed_days": elapsed, "epoch": _calendar_epoch(jd+time/SECONDS_PER_DAY),
                          "distance_km": sample.distance_km,
                          "geometric_angle_deg": math.degrees(sample.geometric_angle_rad),
                          "angular_separation_deg": math.degrees(sample.angular_separation_rad)})
    return replace(track, apses=tuple(apses))
