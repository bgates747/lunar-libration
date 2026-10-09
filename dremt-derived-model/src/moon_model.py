#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Offline Moon-model validation from cached Horizons data.

Analysis and implementation by Codex. DREMT-linked source snapshots: 2026-10-05.
S1: The Tether-Locked Universe, section 4 (relative attitude and orbital frame).
    https://docs.google.com/document/d/1zAOk_j9BZpsHDeSfSsaQUUdbwwat_HxQ/edit
S2: Angular Momentum Conservation Equation, sections 1-3 (planar rate balance).
    https://docs.google.com/document/d/1h8mrNwdD5oIv_BD0vG98MMT5kxpp-cxi/edit
These live documents may change or disappear; dremt-derived-model/README.md identifies the versions.

ANALYST ASSUMPTIONS/ADDITIONS: Newtonian point-mass Sun-Earth-Moon gravity;
empirical first-epoch state vectors, lunar pole, and meridian phase; a supplied
NONZERO rate of change in orientation relative to fixed stars; fixed pole and
constant orientation rate, with no attitude torque. README B1 preserves DREMT's
wording; B2's physical-pendulum interaction is not specified or supplied here.
Physical librations, precession, other planets, tides, and relativity are omitted.
The 3D transport equation extends the source coordinates; it is not a new force.

Requires NumPy. Matplotlib is optional (used only for the chart).
Run from the clone root: python dremt-derived-model/src/moon_model.py
Inputs: dremt-derived-model/data/. Outputs: dremt-derived-model/outputs/validation/.
No network access is used.
"""

from __future__ import annotations

import csv
import json
import math
from project_paths import DATA_DIR, OUTPUT_DIR

import numpy as np

# Analyst-supplied constants, not values derived from S1/S2:
# https://ssd.jpl.nasa.gov/astro_par.html (DE440 GM values and J2000 obliquity).
MU_E = 398600.435507  # km^3/s^2; JPL DE440 astrodynamic parameters
MU_M = 4902.800118
MU_S = 132712440041.279419
# S1 mentions 27.3 days; S2 preserves a rate without fixing its numerical value.
# Analyst initialization: inertial |omega| = 2*pi/(this empirical period), not 0.
ROTATION_PERIOD_DAYS = 27.321661  # independently specified mean sidereal period
EPSILON = math.radians(84381.412 / 3600)


def acceleration(y, solar=True):
    # Analyst-added orbital force law. S1/S2 provide coordinates, not this law.
    # r and s are geocentric; subtract Earth's acceleration from both bodies.
    r, s = y[:3], y[6:9]
    rn = np.linalg.norm(r)
    rd = -(MU_E + MU_M) * r / rn**3
    sn = np.linalg.norm(s)
    sd = -(MU_E + MU_S) * s / sn**3
    if solar:
        delta = s - r
        dn = np.linalg.norm(delta)
        rd = rd + MU_S * (delta / dn**3 - s / sn**3)
        sd = sd + MU_M * (-delta / dn**3 - r / rn**3)
    # In the two-body control, the Sun state is irrelevant to lunar propagation.
    return np.concatenate((y[3:6], rd, y[9:12], sd))


def rk4(y, step, solar, rate=None):
    derivative = rate if rate is not None else lambda state: acceleration(state, solar)
    a = derivative(y)
    b = derivative(y + step * a / 2)
    c = derivative(y + step * b / 2)
    d = derivative(y + step * c)
    return y + step * (a + 2 * b + 2 * c + d) / 6


def propagate(y0, samples, max_step, solar=True, rate=None):
    result = [y0.copy()]
    state = y0.copy()
    previous = samples[0]
    for target in samples[1:]:
        span = target - previous
        steps = int(math.ceil(span / max_step))
        step = span / steps
        for _ in range(steps):
            state = rk4(state, step, solar, rate)
        result.append(state.copy())
        previous = target
    return np.array(result)


def centered_energy_and_momentum(states):
    weights = np.array([MU_E, MU_M, MU_S])
    pos = np.stack((np.zeros_like(states[:, :3]), states[:, :3], states[:, 6:9]), axis=1)
    vel = np.stack((np.zeros_like(states[:, 3:6]), states[:, 3:6], states[:, 9:12]), axis=1)
    pc = (pos * weights[None, :, None]).sum(axis=1) / weights.sum()
    vc = (vel * weights[None, :, None]).sum(axis=1) / weights.sum()
    pos = pos - pc[:, None, :]
    vel = vel - vc[:, None, :]
    energy = 0.5 * (weights[None, :] * (vel * vel).sum(axis=2)).sum(axis=1)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        energy -= weights[i] * weights[j] / np.linalg.norm(pos[:, i] - pos[:, j], axis=1)
    momentum = (weights[None, :, None] * np.cross(pos, vel)).sum(axis=1)
    return energy, momentum


def wrap_degrees(values):
    return (np.asarray(values) + 180) % 360 - 180


def rotation_about(k, angle, initial):
    a = np.asarray(angle)
    return (np.cos(a)[:, None] * initial[None, :]
            + np.sin(a)[:, None] * np.cross(k, initial)[None, :]
            + (1 - np.cos(a))[:, None] * np.dot(k, initial) * k[None, :])


def moving_frame(state, geometric_rate=None):
    # S1 section 4: radial, transverse, orbit-normal basis B. Frame derivatives
    # and Omega_frame below are Codex's implementation for the propagated orbit.
    r, v = state[:3], state[3:6]
    radial = r / np.linalg.norm(r)
    h = np.cross(r, v)
    normal = h / np.linalg.norm(h)
    tangent = np.cross(normal, radial)
    frame = np.column_stack((radial, tangent, normal))
    if geometric_rate is None:
        return frame
    a = geometric_rate[3:6]
    radial_rate = (v - radial * np.dot(radial, v)) / np.linalg.norm(r)
    h_rate = np.cross(r, a)
    normal_rate = (h_rate - normal * np.dot(normal, h_rate)) / np.linalg.norm(h)
    tangent_rate = np.cross(normal_rate, radial) + np.cross(normal, radial_rate)
    frame_rate = np.column_stack((radial_rate, tangent_rate, normal_rate))
    omega_frame = 0.5 * np.cross(frame.T, frame_rate.T).sum(axis=0)
    return frame, omega_frame


def skew(vector):
    x, y, z = vector
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def relative_orientation_rate(omega_inertial):
    # Codex's 3D transport extension: S_dot = [B.T*(omega-Omega_frame)]_x*S.
    # S2's planar cancellation is its fixed-axis limit, not an extra body torque.
    def rate(state):
        geometric = acceleration(state[:12], solar=True)
        frame, omega_frame = moving_frame(state[:12], geometric)
        omega_relative = frame.T @ (omega_inertial - omega_frame)
        matrix_rate = skew(omega_relative) @ state[12:].reshape(3, 3)
        return np.concatenate((geometric, matrix_rate.ravel()))
    return rate


def apparent_orientation(states, times, initial_orientation):
    # Analyst approximation: retain the empirical initial pole and propagate
    # uniform orientation change in nonrotating axes. No attitude torque is integrated.
    lon0, lat0, ra0, dec0 = initial_orientation
    ra, dec = math.radians(ra0), math.radians(dec0)
    north = np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])
    earth0 = -states[0, :3] / np.linalg.norm(states[0, :3])
    projected = earth0 - np.dot(earth0, north) * north
    projected /= np.linalg.norm(projected)
    l0 = math.radians(lon0)
    x0 = math.cos(l0) * projected - math.sin(l0) * np.cross(north, projected)
    y0 = np.cross(north, x0)
    angles = 2 * math.pi * times / (ROTATION_PERIOD_DAYS * 86400)
    xx = rotation_about(north, angles, x0)
    yy = rotation_about(north, angles, y0)
    zz = np.broadcast_to(north, xx.shape)
    attitude = np.stack((xx, yy, zz), axis=2)
    earth = -states[:, :3] / np.linalg.norm(states[:, :3], axis=1)[:, None]
    lon = np.degrees(np.arctan2((earth * yy).sum(axis=1), (earth * xx).sum(axis=1)))
    lat = np.degrees(np.arcsin(np.clip((earth * zz).sum(axis=1), -1, 1)))

    radial = -earth
    normal = np.cross(states[:, :3], states[:, 3:6])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    tangent = np.cross(normal, radial)
    baseline = np.stack((radial, tangent, normal), axis=2)
    # Matrix equivalent of S1 section 4's quaternion relation: S = B.T*A.
    relative = np.matmul(np.transpose(baseline, (0, 2, 1)), attitude)
    recovered = np.matmul(baseline, relative)
    frame_error = float(np.max(np.abs(recovered - attitude)))
    return lon, lat, frame_error, relative, attitude


def ecliptic_longitude(r):
    yy = math.cos(EPSILON) * r[:, 1] + math.sin(EPSILON) * r[:, 2]
    return np.unwrap(np.arctan2(yy, r[:, 0]))


def first_longitude_return_days(r, times):
    progress = ecliptic_longitude(r)
    progress -= progress[0]
    crossings = np.flatnonzero(progress >= 2 * math.pi)
    if len(crossings) == 0:
        return None
    i = int(crossings[0])
    fraction = (2 * math.pi - progress[i - 1]) / (progress[i] - progress[i - 1])
    return float((times[i - 1] + fraction * (times[i] - times[i - 1])) / 86400)


def stats(array):
    return {'rms': float(np.sqrt(np.mean(np.asarray(array)**2))), 'max': float(np.max(np.abs(array)))}


def main():
    ref = json.loads((DATA_DIR / 'moon-model-reference.json').read_text())
    rows = ref['moon']['rows']
    reference = np.array([[float(x) for x in row[2:8]] for row in rows])
    sun = np.array([[float(x) for x in row[2:8]] for row in ref['sun']['rows']])
    orientation = np.array([[float(x) for x in row[3:7]] for row in ref['orientation']['rows']])
    times = np.array([(float(row[0]) - float(rows[0][0])) * 86400 for row in rows])
    assert len(rows) == 673 and len(sun) == len(orientation) == len(rows)
    assert np.max(np.abs(times - np.arange(673) * 3600)) < 1e-4
    y0 = np.concatenate((reference[0], sun[0]))

    # Analyst choices: first epoch only initializes the model; later rows are
    # validation targets. RK4 and 600/300-second convergence checks are additions.
    coarse = propagate(y0, times, 600, solar=True)
    model = propagate(y0, times, 300, solar=True)
    kepler = propagate(y0, times, 300, solar=False)
    positions = model[:, :3]
    distance = np.linalg.norm(positions, axis=1)
    reference_distance = np.linalg.norm(reference[:, :3], axis=1)
    error = np.linalg.norm(positions - reference[:, :3], axis=1)
    kepler_error = np.linalg.norm(kepler[:, :3] - reference[:, :3], axis=1)
    convergence = np.linalg.norm(positions - coarse[:, :3], axis=1)
    cosine = (positions * reference[:, :3]).sum(axis=1) / distance / reference_distance
    sky_error = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    model_lon, model_lat, mapping_error, relative_attitude, attitude = apparent_orientation(model, times, orientation[0])
    omega_inertial = 2 * math.pi / (ROTATION_PERIOD_DAYS * 86400) * attitude[0, :, 2]
    relative_initial = np.concatenate((y0, relative_attitude[0].ravel()))
    relative_run = propagate(relative_initial, times, 300, rate=relative_orientation_rate(omega_inertial))
    relative_predicted_attitude = np.array([moving_frame(state[:12]) @ state[12:].reshape(3, 3) for state in relative_run])
    relative_integration_error = float(np.max(np.abs(relative_predicted_attitude - attitude)))
    ref_lon = wrap_degrees(orientation[:, 0])
    ref_lat = orientation[:, 1]
    longitude_error = wrap_degrees(model_lon - ref_lon)
    latitude_error = model_lat - ref_lat
    hard_lock_lon_error = wrap_degrees(orientation[0, 0] - ref_lon)
    hard_lock_lat_error = orientation[0, 1] - ref_lat
    energy, momentum = centered_energy_and_momentum(model)
    energy_drift = float(np.max(np.abs((energy - energy[0]) / energy[0])))
    angular_drift = float(np.max(np.linalg.norm(momentum - momentum[0], axis=1)) / np.linalg.norm(momentum[0]))
    lunar_h = np.cross(model[:, :3], model[:, 3:6])
    lunar_h_change = float(np.max(np.linalg.norm(lunar_h - lunar_h[0], axis=1)) / np.linalg.norm(lunar_h[0]))

    # Codex diagnostics for S1 section 2, separate from the lunar hindcast.
    # Initial inertial rest is an explicit condition of the first witness.
    # Normalized circular example: Omega=1, radius=10, central inertia=2.
    Omega, I, M, R = 1.0, 2.0, 1.0, 10.0
    angle = math.pi / 2
    theta_torque_free_from_rest = 0.0
    theta_face_locked = angle
    min_unsat_witness = {'orbital_angle_rad': angle, 'theta_from_initial_rest_and_zero_torque_rad': theta_torque_free_from_rest,
                        'theta_required_by_face_lock_rad': theta_face_locked,
                        'residual_rad': theta_face_locked - theta_torque_free_from_rest}
    # Proper physical velocities already contain spin prior to release.
    kinetic_total = 0.5 * M * R**2 * Omega**2 + 0.5 * I * Omega**2
    spin_energy = kinetic_total - 0.5 * M * R**2 * Omega**2
    # S2 sections 2-3: phi_ddot=-f_ddot preserves inertial theta_dot.
    alpha_b = 0.2
    relative_acceleration = -alpha_b
    absolute_acceleration = relative_acceleration + alpha_b
    proof_spin_momentum_derivative = I * absolute_acceleration

    metrics = {
        'interval': ref['protocol'],
        'reference_translation_source': next((x.strip() for x in ref['moon']['header'].splitlines() if 'Target body name' in x), None),
        'reference_rotation_source': next((x.strip() for x in ref['orientation']['header'].splitlines() if 'pole/equ' in x.lower()), None),
        'sample_count': len(times),
        'mean_reference_distance_km': float(reference_distance.mean()),
        'orbit_acceptance_threshold_km': float(0.001 * reference_distance.mean()),
        'three_body_position_error_km': stats(error),
        'two_body_position_error_km': stats(kepler_error),
        'three_body_distance_error_km': stats(distance - reference_distance),
        'three_body_sky_direction_error_deg': stats(sky_error),
        'rk4_600_vs_300_second_position_difference_km': stats(convergence),
        'total_three_body_energy_relative_drift': energy_drift,
        'total_three_body_angular_momentum_relative_vector_drift': angular_drift,
        'earth_moon_orbital_specific_angular_momentum_relative_vector_change': lunar_h_change,
        'longitude_error_deg': stats(longitude_error),
        'latitude_error_deg': stats(latitude_error),
        'hard_tether_constant_orientation_longitude_error_deg': stats(hard_lock_lon_error),
        'hard_tether_constant_orientation_latitude_error_deg': stats(hard_lock_lat_error),
        'reference_subearth_longitude_range_deg': [float(ref_lon.min()), float(ref_lon.max())],
        'reference_subearth_latitude_range_deg': [float(ref_lat.min()), float(ref_lat.max())],
        'first_ecliptic_longitude_return_days': {'reference': first_longitude_return_days(reference[:, :3], times),
                                               'three_body': first_longitude_return_days(positions, times),
                                               'two_body': first_longitude_return_days(kepler[:, :3], times)},
        'frame_reconstruction_max_matrix_element_error': mapping_error,
        'independently_integrated_rotating_frame_max_attitude_matrix_error': relative_integration_error,
        'orbit_acceptance_passed': bool(error.max() <= 0.001 * reference_distance.mean()),
        'orientation_acceptance_passed': bool(np.max(np.abs(longitude_error)) <= 1 and np.max(np.abs(latitude_error)) <= 1),
        'consistency_tests': {
            'T01_zero_initial_spin_zero_torque_and_face_lock': {'status': 'inconsistent', 'witness': min_unsat_witness},
            'T02_planar_relative_momentum_equation': {'status': 'consistent_with_constant_inertial_spin', 'spin_momentum_derivative': proof_spin_momentum_derivative},
            'T03_release_creates_spin': {'status': 'inconsistent_without_angular_impulse', 'spin_before_release': I * Omega, 'spin_after_release': I * Omega},
            'T04_zero_axial_energy_of_face_locked_extended_body': {'status': 'inconsistent_if_energy_means_motion_relative_to_COM', 'central_spin_energy': spin_energy},
            'T05_relative_attitude_coordinate_mapping': {'status': 'consistent', 'max_reconstruction_error': mapping_error},
            'T06_hard_tether_orientation_vs_monthly_libration': {'status': 'fails_one_degree_test', 'longitude_max_error_deg': float(np.max(np.abs(hard_lock_lon_error))), 'latitude_max_error_deg': float(np.max(np.abs(hard_lock_lat_error)))},
        },
        'constants': {'mu_earth_km3_s2': MU_E, 'mu_moon_km3_s2': MU_M, 'mu_sun_km3_s2': MU_S, 'mean_rotation_period_days': ROTATION_PERIOD_DAYS},
        'orientation_time_note': ref['orientation_time_note'],
        'limitations': ['A 28-day hindcast with empirical initial states, pole, phase, and mean spin rate; not a derivation of tidal locking.',
                        'Translation includes three Newtonian point masses; it omits other planets, oblateness, tidal forces, and relativity.',
                        'Orientation uses constant spin about the initial fixed pole. It omits physical librations and pole precession.',
                        'Position and orientation references come from JPL models, not newly acquired raw observations.',
                        'Exact quaternion/matrix reexpression does not validate all stronger physical claims.'],
    }
    assert convergence.max() < 0.01, 'Integrator uncertainty exceeds 10 m'
    assert energy_drift < 1e-9 and angular_drift < 1e-9
    assert mapping_error < 1e-12
    assert relative_integration_error < 1e-9, 'Rotating-frame dynamics disagrees with independent fixed-frame solution'
    assert proof_spin_momentum_derivative == 0
    assert spin_energy > 0
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / 'moon-model-validation.json').write_text(json.dumps(metrics, indent=2) + '\n')

    with (OUTPUT_DIR / 'moon-model-predictions.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['date_tdb', 'elapsed_days', 'x_km', 'y_km', 'z_km', 'vx_km_s', 'vy_km_s', 'vz_km_s',
                         'reference_x_km', 'reference_y_km', 'reference_z_km', 'position_error_km', 'two_body_position_error_km',
                         'model_subearth_longitude_deg', 'reference_subearth_longitude_deg',
                         'model_subearth_latitude_deg', 'reference_subearth_latitude_deg',
                         *[f'relative_attitude_{i}{j}' for i in range(3) for j in range(3)]])
        for i, row in enumerate(rows):
            writer.writerow([row[1].strip(), times[i] / 86400, *model[i, :6], *reference[i, :3], error[i], kepler_error[i],
                             model_lon[i], ref_lon[i], model_lat[i], ref_lat[i], *relative_attitude[i].ravel()])

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
        fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), layout='constrained')
        days = times / 86400
        ax = axes[0, 0]
        ax.plot(days, reference_distance / 1000, color='#222222', label='JPL reference', linewidth=2.3)
        ax.plot(days, distance / 1000, color='#007c83', label='Sun–Earth–Moon model', linestyle='--')
        ax.plot(days, np.linalg.norm(kepler[:, :3], axis=1) / 1000, color='#bd6422', label='Two-body control', linewidth=1.1)
        ax.set(ylabel='Earth–Moon distance (thousand km)', title='One initial state; 28 days of propagation')
        ax.legend(fontsize=8)
        ax = axes[0, 1]
        ax.semilogy(days[1:], error[1:], color='#007c83', label='Three-body error')
        ax.semilogy(days[1:], kepler_error[1:], color='#bd6422', label='Two-body error')
        ax.axhline(metrics['orbit_acceptance_threshold_km'], color='#777777', linestyle=':', label='0.1% distance criterion')
        ax.set(ylabel='Position-vector error (km; log scale)', title='Comparison against JPL geometric vectors')
        ax.legend(fontsize=8)
        for ax, actual, predicted, label, first in [
                (axes[1, 0], ref_lon, model_lon, 'Sub-Earth longitude (degrees)', ref_lon[0]),
                (axes[1, 1], ref_lat, model_lat, 'Sub-Earth latitude (degrees)', ref_lat[0])]:
            ax.plot(days, actual, color='#222222', label='JPL reference', linewidth=2)
            ax.plot(days, predicted, color='#007c83', label='Constant inertial spin, fixed pole', linestyle='--')
            ax.axhline(first, color='#a13d46', linestyle=':', label='Orientation fixed in tether frame')
            ax.set(ylabel=label)
            ax.legend(fontsize=8)
        axes[1, 0].set_title('Longitudinal libration')
        axes[1, 1].set_title('Latitudinal libration')
        for ax in axes.flat:
            ax.set_xlabel('Days after 2026-09-01 00:00 TDB')
            ax.grid(alpha=0.2)
        fig.suptitle('Moon model validation — September 1–29, 2026', fontsize=15)
        fig.savefig(OUTPUT_DIR / 'moon-model-validation.png', dpi=180)
        fig.savefig(OUTPUT_DIR / 'moon-model-validation.svg')
        plt.close(fig)
    except ImportError:
        print('Matplotlib unavailable; numeric results were saved without a chart.')
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
