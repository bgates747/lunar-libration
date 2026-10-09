# Moon longitudinal libration

A top-down, physically timed demonstration of the Moon's optical longitudinal
libration. The orbit follows Kepler's equation while the Moon rotates uniformly
once per sidereal orbital period. Axial tilt, latitude libration, perturbations,
and physical libration are intentionally omitted.

## Model

- Semi-major axis: 384,400 km
- Eccentricity: 0.0549
- Sidereal period: 27.321661 days
- Earth fixed at one focus
- Earth and Moon are exaggerated for legibility, but retain their real 3.67:1
  radius ratio; their displayed distance uses a separate scale

At each frame, mean anomaly advances uniformly. Newton's method solves
`M = E - e sin(E)`, and the focus-centered position is calculated from the
eccentric anomaly. The instantaneous speed comes from the vis-viva equation.

## Set up and render

Run these commands from the **`lunar-libration/` clone root**. Python 3.11+
and FFmpeg are required. All Python packages are isolated in the
project venv; no system Python packages are changed.

```bash
python3 -m venv .venv
.venv/bin/pip install -e './author-2d-animation[test]'
.venv/bin/python -m pytest author-2d-animation/tests
.venv/bin/moon-libration --start perigee
```

The default render is 1920×1080 at 60 fps and is written to
`author-2d-animation/outputs/moon_libration_1080p.mp4`. A quick proof render can be made with:

```bash
.venv/bin/moon-libration \
  --duration 5 --hold 0.25 --fps 24 --resolution 854,480 \
  --output author-2d-animation/outputs/moon_libration_preview.mp4
```

Frames are streamed directly to FFmpeg, so no temporary image sequence is kept.

Relative output paths are resolved from the clone root, independently of the
current working directory. Generated files stay in `author-2d-animation/outputs/`.
