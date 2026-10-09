# Earth-centred angular comparison

[`animate.py`](../src/animate.py) uses the visual style of `author-2d-animation/`: a dark
background, Earth–Moon trajectory, instantaneous velocity arrow and readout
panel. The angle is measured at **Earth’s centre of mass** between two radial
lines. This version replaces the earlier body-centred attitude comparison.
The earlier [`animation_model.py`](../src/animation_model.py) numerical model remains available separately.
The pendulum view is on the left, the top-down orbital view is in the centre,
and the numerical dashboard is on the right.
All rendered versions use the title **Pendular Libration Model**.

## Ideal Keplerian version

Select `--orbit keplerian` to use the original `author-2d-animation` calculation
directly through [`keplerian_comparison.py`](../src/keplerian_comparison.py). Both model folders must remain in
their sibling locations. This version uses no ephemerides. Its supplied inputs
are a = 384,400 km, e = 0.0549 and T = 27.321661 days. As in the original,
these a/T values determine an effective μ = 4π²a³/T²; this is an idealized
orbit, not a fit using the exact measured Earth–Moon gravitational parameter.

**Supplied orientation and phase:** C = n = 2π/T is nonzero, and θ(0) = π
faces Earth at starting perigee. Thus θ(t) = π + nt. These are analyst inputs,
not initial conditions determined by the coordinate rate identity.

Kepler's equation M = E − e sin E gives E(t), where M = nt. The model Moon
has x = a(cos E − e), y = a√(1 − e²) sin E. The virtual Moon shares its
instantaneous radius and has bearing M. Their Earth-centred separation is
Δ = M − f, the negative of the original animation's displayed angle f − M.
With zero initial virtual offset, the Moons align at both perigees and apogee;
the peak separations are equal and opposite, approximately ±6.293253°.
The curve is nearly sinusoidal but is evaluated from Kepler's equation.

This implements the relative-motion consequence of S2's ω̇ₛ = −Ω̇ with
the specified C and phase. The rotation and translation used in the left view
are unchanged. No physical pendulum suspension, restoring torque or transfer
of orbital energy into rocking is supplied. The ephemeris-derived orientation
and latent latitude data described below apply only to the ephemeris version;
the ideal version is planar with zero latitude and uniform orientation.

```bash
.venv/bin/python dremt-derived-model/src/animate.py --orbit keplerian
```

Its default output is `dremt-derived-model/outputs/animation/pendular_libration_keplerian_1080p.mp4`.

## What it displays

The solid grey **real Moon** follows all the cached NASA/JPL Horizons
Earth-centred position and velocity samples in [`animation-ephemeris.json`](../data/animation-ephemeris.json).
Hourly position/velocity nodes are interpolated with cubic Hermite polynomials
for smooth playback. This is an ephemeris playback, not a prediction produced
by DREMT’s equations or by the separately supplied Newtonian numerical model.
Its cyan tether connects its centre directly to Earth’s COM. Its only arrow
originating at the Moon is the green instantaneous velocity vector, projected
into the camera plane at a fixed scale.

The **virtual Moon** is a magenta dotted circle with a magenta dotted line
connecting its centre to Earth’s COM. It has no orientation markers, surface
features or velocity arrow. Its projected bearing advances uniformly:

αᵥ(t) = αᵣ(0) + δ₀ + nt,

n = 2π / (27.321661 days) ≈ 13.17635849°/day.

This is the supplied mean sidereal orbital rate, rather than an instantaneous
body-orientation rate about the tilted lunar pole. By default δ₀ = 0: the two
Moons begin at the same bearing. This phase choice is a visualization
convention, not a determination of the Moon’s mean anomaly or mean longitude.
Use `--virtual-phase-deg` to supply another initial offset. Changing that offset
shifts the separation without changing either Moon’s angular-rate history.
Consequently the separation need not be symmetric about zero.

**Radial assumption:** the virtual Moon shares the real Moon’s instantaneous
projected Earth distance. Its hidden height above the camera plane also
matches the real Moon, so their full 3D Earth distances agree. Their tether
lengths vary with the ephemeris distance; neither tether is a rigid arm of
constant length. This shared radial history isolates the angular comparison.

The gold arc is centred on Earth’s COM. Its signed angle is

Δ(t) = wrap[αᵥ(t) − αᵣ(t)],

with counterclockwise positive and the displayed value wrapped to −180°…180°.
It measures **positional angular separation**, not rocking about the Moon’s
COM and not the Moon’s body-frame libration longitude. The two displayed rates
are azimuth rates in the same fixed camera plane.

## Left pendulum view

The left panel displays the **same sampled motion at exactly the same body
and distance scales** as the main view. It is a rotation within the displayed
plane, together with a translation, rather than a second simulation or a new
orbital plane. The magenta virtual Moon is stationary and its line to Earth
is perfectly vertical. The real Moon moves from side to side relative to it.

[`pendulum_view.py`](../src/pendulum_view.py) rotates by β(t) = −π/2 − αᵥ(t) and anchors the virtual Moon:

r′ᵣ = R(β)(rᵣ − rᵥ),

r′ᵥ = 0,

r′Earth = −R(β)rᵥ.

All vectors in these expressions are the main view's two-dimensional projected
vectors. The translation keeps the virtual Moon’s **centre**, as well as its
bearing, fixed. Because the radial distance varies, Earth pans slightly up and
down in this view. Distances and the angular separation at Earth’s COM remain
unchanged. Both radial lines meet at Earth's displayed COM; there are no
body-orientation arrows or velocity arrows in the left panel.

The identical Earth graphic is fixed in the main view and turns by β in the
left view. Its surface glyphs mark the **fixed inertial axes** of the main
view; this apparent turning does not model Earth's physical daily rotation.
The real Moon’s decorative crater cues receive the same view rotation. The
virtual Moon remains only a dotted magenta outline and dotted radial line.
The apogee label is placed inside the main orbital path to make room.

This is a coordinate visualization of the relative motion. It does not supply
a physical pendulum suspension, restoring torque or energy-transfer mechanism.

## Perigee start and geometric coordinates

The renderer selects the **first perigee in the cache with a complete projected
revolution available after it**. Perigee is found as a negative-to-positive
zero crossing of radial velocity, (r · v)/|r|, refined within the hourly cubic
Hermite interpolation. It is not rounded to the closest hourly table row.

The fixed camera basis is constructed from r and v at that perigee. Its +x
axis points from Earth to the starting perigee, so the Moon begins on the right
at **geometric angle 0°**, with the virtual Moon aligned. Its normal is the
initial orbital normal. The geometric-angle readout then increases continuously
until **360°**, the first full return to that initial bearing.

The end is not forced to be the next perigee. The real Moon’s varying apsidal
direction is preserved: the bundled cache gives a later perigee at about
328.37°, before the full revolution ends. Likewise, the actual apogee need
not appear at exactly 180°, and the virtual Moon need not align at later apses.
The mean virtual rate is unchanged; it is not fitted to the selected interval.
No angular extrema are symmetrized and no final alignment is imposed.

With the bundled cache, playback begins September 6, 2026 at
20:43:53.921 TDB and reaches 360° on October 4, 2026 at 02:25:35.309 TDB,
covering approximately 27.237285 days. The marked apogee is September 19;
the next perigee is October 1. These are interpolated geometric events in the
cached ephemeris, not UTC timestamps or new observations.

## Three-dimensional data and scope

The camera looks along the starting-perigee orbital normal and stays fixed
throughout the selected revolution. [`orbit_comparison.py`](../src/orbit_comparison.py) retains the real 3D position and
velocity, virtual 3D position, ephemeris-derived body attitude, and sub-Earth
longitude and latitude. A later latitude display can use these quantities;
latitude is not displayed here. The real Moon’s crater glyphs are decorative
surface cues, not a rendered terrain map. Its pole and meridian phase are
interpolated from the cached orientation rows, unwrapping circular angles.
The reference file records the negligible TT/TDB timing difference between
its vector and orientation grids.

Earth is kept at the displayed geocentric origin. Body radii are enlarged
separately from orbital distances, retaining the Earth/Moon radius ratio.

## Provenance and analyst choices

DREMT’s baseline along the body-to-host-COM line is documented in
[S1, sections 1 and 4](https://docs.google.com/document/d/1zAOk_j9BZpsHDeSfSsaQUUdbwwat_HxQ/edit)
and [S2, section 1](https://docs.google.com/document/d/1h8mrNwdD5oIv_BD0vG98MMT5kxpp-cxi/edit).
The [main README](../README.md) records the October 5, 2026 document snapshots;
the live links may change or disappear. The longer animation cache contains
hourly September 1–October 10, 2026 data, retrieved October 6, 2026. It records
its retrieval time, requests, coordinate conventions and source headers.
The original 28-day [`moon-model-reference.json`](../data/moon-model-reference.json) and numerical model-validation
inputs are unchanged. Requests use the documented
[NASA/JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html).

**Analyst choices:** the virtual Moon, its independently supplied mean rate,
initial phase, shared radial history, perigee-based fixed camera and Earth-centred angle
comparison were added at the repository author’s direction. They are not a
physical pendulum force law derived from the source documents. Drawing both
lines through Earth’s COM does not establish a material hinge, suspension
constraint, torque or orbital-to-rocking energy transfer. No private source
archives or full blog comments are required by the renderer.

## Render

Run all commands from the **`lunar-libration/` clone root**. Install
`dremt-derived-model/requirements/requirements-animation.txt` in the same
environment as the model:

```bash
python3 -m venv .venv
.venv/bin/pip install -r dremt-derived-model/requirements/requirements-animation.txt
.venv/bin/python dremt-derived-model/src/animate.py
```

The default is a 1920×1080, 60 fps, 33-second silent MP4: 30 seconds for the
selected full revolution, with a 1.5-second hold at each endpoint. It is written to
`dremt-derived-model/outputs/animation/dremt_pendulum_left_1080p.mp4`, alongside a JSON record of
source, settings, selected epochs, actual apses, phase, mean rate and angular range. Generated artifacts are
excluded from version-control selection. A shorter preview:

```bash
.venv/bin/python dremt-derived-model/src/animate.py --resolution 854,480 --fps 24 \
  --duration 5 --hold 0.25 --output dremt-derived-model/outputs/animation/preview.mp4
```

Use `--frame-day 0 --output dremt-derived-model/outputs/animation/perigee.png` for a starting still.
Relative `--output` and `--reference` paths are resolved from the clone root,
independently of the current working directory. Explicit absolute paths are
also accepted when supplied by the caller. Use `--reference PATH` to select
another compatible cache containing a full
revolution after a perigee. To explicitly refresh the longer bundled cache:

```bash
.venv/bin/python dremt-derived-model/scripts/fetch_animation_ephemeris.py
```

The renderer itself performs no downloads. Use
`--ffmpeg /path/to/ffmpeg` to select an existing encoder; otherwise the script
uses an installed FFmpeg or the `imageio-ffmpeg` bundled encoder. Rendering is
offline once dependencies are installed. The calculation modules are usable
without Pillow or the encoder package.
