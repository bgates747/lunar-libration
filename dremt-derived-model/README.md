# Moon model in a moving Earth–Moon frame

This derivation and implementation were prepared by Codex. The repository author
identified the sources and requested the analysis. The model retains the sources’
relative-orientation coordinates and adds a specified orbital force law and
empirical initial conditions.

The [provisional findings pending further information](REPORT-2026-10-09.md) present the
published equations with the relevant exact comments. The preceding textual
audit appears as Appendix A. A
[dedicated information-needed section](REPORT-2026-10-09.md#10-information-needed-to-resolve-the-provisional-findings)
is linked from the relevant mathematical and textual passages.
The findings distinguish the valid coordinate description and testable
constant-rate approximation from the unestablished pendulum-energy explanation.

Dated findings are retained as separate files in this directory. The
[October 9, 2026 snapshot](REPORT-2026-10-09.md) is preserved unchanged.
Future revisions will use new dated filenames (with a distinct suffix for any
additional same-day version), retain the earlier files, and be listed here.

## Source provenance

The primary sources are these two DREMT-linked documents:

| ID | Document | Downloaded |
|---|---|---|
| S1 | [The Tether-Locked Universe: A Unified Kinematic & Dynamic Synthesis](https://docs.google.com/document/d/1zAOk_j9BZpsHDeSfSsaQUUdbwwat_HxQ/edit) | October 5, 2026 |
| S2 | [Angular Momentum Conservation Equation in the Tether-Locked Framework](https://docs.google.com/document/d/1h8mrNwdD5oIv_BD0vG98MMT5kxpp-cxi/edit) | October 5, 2026 |

The download dates use the UTC calendar dates of the local archive files’
modification timestamps. They identify the captured versions, rather than an
exact download instant. These are live Google Docs links: their contents may
change, access may be withdrawn, and the documents may disappear. This analysis
refers to the October 5 snapshots, not necessarily their future contents.

Two comments on Roy Spencer’s blog are cited below to clarify DREMT’s own
description of the combined motion and his claim about its physical cause.
Their timestamps are reproduced exactly as displayed by the site; the site
does not specify their timezone. These comment links may also be ephemeral.

The full documents and their equation images are not distributed with this
project. Short attributed excerpts identify the particular claims being
examined; the derivation and numerical implementation are original commentary.
The [MIT licence](../LICENSE) covers our original code and documentation, not the
quoted source material or the linked works. Attribution and excerpt length do
not guarantee a fair-use outcome; the [U.S. Copyright Office](https://www.copyright.gov/fair-use/)
describes the contextual assessment.

## Analyst assumptions and additions

These are Codex’s modeling choices, not results derived from S1 or S2:

- The orbit is propagated with Newtonian gravity between point-mass Sun, Earth,
  and Moon. Other planets, nonspherical gravity, tides, and relativity are omitted.
- The Moon’s orientation changes relative to the fixed stars at a constant,
  **nonzero rate**, using the supplied empirical period **27.321661 days**,
  about its initial fixed pole. Physical librations and pole precession are omitted.
- Only the first-epoch Moon and Sun state vectors, lunar pole, and meridian phase
  initialize the trajectory and attitude. Later reference samples are comparison
  targets; they do not drive the propagation or fit the parameters.
- The planar derivation assumes a fixed rotation axis and constant central
  moment of inertia I > 0. The 3D implementation transports attitude matrices
  explicitly instead of applying a planar scalar equation to arbitrary axes.
- RK4 integration, 300-second steps, a 600-second convergence comparison, and
  the acceptance criteria in `dremt-derived-model/config/moon-model-test-protocol.json` are analyst choices.

The gravitational parameters and J2000 ecliptic obliquity are taken from
[JPL’s astrodynamic parameters](https://ssd.jpl.nasa.gov/astro_par.html).
`dremt-derived-model/data/moon-model-reference.json` contains NASA/JPL Horizons request parameters,
source metadata, and numerical samples. Its recorded retrieval time is
October 5, 2026 at 21:55:19 UTC. These numerical ephemeris data have a separate
provenance from the DREMT-linked documents.

## Derivation from the source definitions

### DREMT’s description of the combined motion

DREMT wrote on [September 24, 2026 at 10:42 AM](https://www.drroyspencer.com/2026/09/uah-v6-1-global-temperature-update-for-august-2026-0-68-deg-c/#comment-1755253)
(B1):

> Basically, “orbit without spin”, in an ellipse, is motion where the same side of the body always faces along the “tether line” of gravity from the orbiting body to the orbited body. In addition to this, there is a secondary motion, a slight side to side “wobbling” of the orbiting body on its own axis, relative to the “tether line”. Put these motions together, and it perfectly matches the movement of the moon, where the change in orientation proceeds at a constant rate relative to the fixed stars, and the position changes at a variable rate.

The quantity initialized in this implementation is the **rate of change in
orientation relative to the fixed stars**. DREMT describes the secondary
wobbling relative to the tether line and the combined orientation change
relative to the stars. These are distinct angular rates. The derivation below
keeps that distinction explicit instead of assigning his word “spin” to both.

### Planar rate balance

Use Ω = ω_b = ḟ for the orbital-line rate, α_b = Ω̇, θ for the body’s orientation
in nonrotating axes, and φ = θ − f for its relative angle. The source’s relative
rate is ω_s = φ̇. Here I is the body’s central moment of inertia.

S2, section 2, states:

> In an inertial stellar frame, absolute angular momentum is modeled as I(ω_b + ω_s).

Its displayed momentum equation, in sections 2–3, is

```text
Iω̇_s = −Iα_b, or Iω̇_s + IΩ̇ = 0.
```

For constant I, Codex’s integration of that equation gives

```text
Ω + ω_s = C
θ̇ = C
θ(t) = θ₀ + C(t − t₀)
φ(t) = θ₀ + C(t − t₀) − f(t)
```

The integration constant C is the combined rate of change in orientation
relative to the fixed stars, θ̇. It need not equal
ω_s at the initial epoch. S1, section 1, mentions the conventional 27.3-day
period; S2’s differential equation does not itself determine a numerical C.
This implementation supplies C = 2π/(27.321661 × 86,400 s). For an ideal closed
orbit with bounded relative rocking and zero net relative turns, C = 2π/T
also follows from that additional periodicity condition.

This interpretation agrees with S1, section 3:

> The relative back-and-forth shaking is an artifact of the tether's changing velocity.

S1, section 5, also attributes a driving role to the effective inertial torque.
That causal interpretation must be distinguished from the coordinate balance.
In the same planar approximation, the physical torque about the lunar centre
of mass is

```text
τ_C = Iθ̈ = I(φ̈ + f̈)
Iφ̈ = τ_C − If̈
```

Thus S2’s cancellation φ̈ = −f̈ is the τ_C = 0 case. The term −If̈ is required
in the accelerating frame; it does not establish an additional interaction or
an exchange of orbital energy with central rotational energy.

DREMT also wrote on [September 28, 2026 at 5:40 AM](https://www.drroyspencer.com/2026/09/uah-v6-1-global-temperature-update-for-august-2026-0-68-deg-c/#comment-1755462)
(B2):

> These torques would cause the moon to physically rock back and forth on its attachment point like a giant pendulum in space.

In that comment he identifies inertial effects across the Moon’s asymmetric
mass distribution and proposes a diversion of orbital energy into the rocking.
The sources do supply an effective torque formula, −Iα_b, and DREMT explicitly
invokes distributed inverse-square gravity. What they do not supply is an
evaluated material torque, work, storage and orbital-reaction account deriving
the proposed transmission from that interaction. The effective frame term alone
does not supply those results. Zero resultant torque can coexist with opposing
physical torques and power transmission; an intended lunar interaction pair
would still have to be identified and calculated. [REPORT-2026-10-09.md](REPORT-2026-10-09.md)
examines the exact statements and their qualifications.

Relative rocking φ̈ ≠ 0 does not by itself imply an additional physical central
torque: under the source rate balance, it coexists with θ̈ = 0. To establish
B2’s stronger causal interpretation, one must specify the material forces,
evaluate their moments about the lunar centre of mass, and demonstrate their
work. This implementation adds Newtonian orbital forces but does not invent
an attitude-driving force to complete the pendulum claim.

### Three-dimensional coordinate reconstruction

S1, section 4, defines

```text
q_spin(t) = q_baseline⁻¹(t) ⊗ q_actual(t)
```

Its basis is built from the orbital position and velocity. Using the code’s
notation, the same basis is

```text
eᵣ = r/‖r‖
eₕ = (r × v)/‖r × v‖
eₜ = eₕ × eᵣ
B = [eᵣ eₜ eₕ]
```

The matrix equivalent is S = BᵀA and A = BS, where A maps body axes into
inertial axes and S maps them into the moving orbital frame. These are source
coordinate definitions expressed with matrices, not an added force law.

Codex’s 3D transport implementation computes the moving-frame angular velocity
Ω_frame from the basis derivatives, then integrates

```text
ω_rel = Bᵀ(ω_inertial − Ω_frame)
Ṡ = [ω_rel]×S
```

Here [u]× maps a vector w to u × w. The independently supplied inertial rate is
ω_inertial = Ck₀, with fixed initial pole k₀. Reconstructing A = BS can therefore
be compared with a separate attitude calculation using constant change in
orientation relative to nonrotating axes.

### Added orbital dynamics

Let r be the Moon’s geocentric position, s the Sun’s geocentric position, and
μ_E, μ_M, μ_S their gravitational parameters. Codex adds the equations obtained
by subtracting Earth’s Newtonian acceleration from each body’s acceleration:

```text
r̈ = −(μ_E + μ_M)r/‖r‖³ + μ_S[(s − r)/‖s − r‖³ − s/‖s‖³]
s̈ = −(μ_E + μ_S)s/‖s‖³ + μ_M[(r − s)/‖r − s‖³ − r/‖r‖³]
```

Both relative states are propagated together. The changing orbital-line rate
then comes from r(t) and v(t); no separate libration waveform is imposed.
These dynamics are additions to the source coordinates. They are not inferred
from the quaternion identity or the local momentum equation.

## Files and execution

All commands below run from the **`lunar-libration/` clone root**. Source modules
locate bundled files relative to that root, independently of the working
directory. Relative `--reference` and `--output` arguments use the same root.

```text
dremt-derived-model/
├── README.md
├── REPORT-2026-10-09.md
├── src/           # propagation, attitude, comparison geometry and renderer
├── data/          # cached NASA/JPL reference and animation ephemerides
├── config/        # validation protocol and acceptance criteria
├── requirements/  # core, plotting, animation and test dependencies
├── docs/          # animation explanation and rendering instructions
├── scripts/       # explicit network cache refresh
└── tests/         # numerical and geometry checks
```

- [Model and validation](src/moon_model.py) use the
  [28-day reference data](data/moon-model-reference.json) and
  [test protocol](config/moon-model-test-protocol.json).
- [Animation instructions](docs/ANIMATION.md) distinguish ephemeris playback
  from the ideal Keplerian comparison and document the rendering assumptions.
- [Dependency files](requirements/requirements.txt) keep optional plotting,
  animation and test packages separate from the core NumPy calculation.

Python 3.11+ is recommended:

```bash
python3 -m venv .venv
.venv/bin/pip install -r dremt-derived-model/requirements/requirements.txt
.venv/bin/python dremt-derived-model/src/moon_model.py
```

For charts, also install
`dremt-derived-model/requirements/requirements-plot.txt`. Once dependencies are
installed, validation is offline. Generated JSON, CSV and optional PNG/SVG
charts go to `dremt-derived-model/outputs/validation/`. Rendered videos, stills
and rendering metadata go to `dremt-derived-model/outputs/animation/` by default.
The entire generated `outputs/` directory is ignored by version control.
It is created on demand and is not required in a clone.

Run the existing checks with:

```bash
.venv/bin/pip install -r dremt-derived-model/requirements/requirements-test.txt
.venv/bin/python -m pytest dremt-derived-model/tests
```

## What works

The relative-coordinate transformation is consistent and reconstructs the same
body orientation. Independently integrating the moving-frame attitude agrees
with the calculation using a constant rate of change in orientation relative
to the fixed stars to approximately 1.75 × 10⁻¹⁴ per matrix element over this interval.

With the analyst-supplied gravity law and initial conditions, the model meets
the specified position and orientation criteria for September 1–29, 2026:

| Check | Maximum error against NASA/JPL reference |
|---|---|
| Three-dimensional position | 1.3989 km |
| Geocentric sub-Earth longitude | 0.02247° |
| Geocentric sub-Earth latitude | 0.06965° |

These angular figures are prediction errors, not libration amplitudes. The
comparison uses 673 hourly samples. The two-body control has a maximum position
error of approximately 16,535 km, demonstrating the importance of including
solar perturbations for this particular monthly comparison.

The scalar rate balance works as a description of constant change in orientation
relative to the fixed stars, viewed from an orbital frame whose rate varies.
The 3D implementation also handles latitude geometry using the supplied lunar
pole and orbital state.

## What does not work

**The central unresolved issue is a mismatch between the causal explanation
and what the equations establish.** The relative rocking is well defined;
the proposed material transmission and pendulum-energy conversion are not
derived by the published balance. This conclusion preserves DREMT’s explicit
non-inertial-frame qualification and his invocation of distributed gravity.

The coordinate identity and local rate balance do not by themselves predict
an orbit, determine the observed numerical rate of change in orientation
relative to the fixed stars, or derive tidal synchronization. They preserve
θ̇ = C; its numerical value requires initialization or an additional condition.
They require additional physical dynamics, initial conditions, and attitude
information. S2’s central-momentum expression alone is not a proof of conservation
of the entire multibody system’s orbital angular momentum and energy.

S1, section 2, claims that maintaining a fixed stellar orientation requires
continuous counter-torque. Under its stated zero-central-torque premise, an
initially nonrotating body instead retains that orientation without such a
counter-torque. Exact face-locking in an eccentric orbit requires
θ̈ = f̈ and therefore τ_C = If̈ wherever f̈ ≠ 0. A body already initialized in
synchronous rotation can remain face-locked without torque in the ideal circular
case; that does not establish S1’s proposed uniquely torque-free baseline.

Zero net relative turns do not imply zero inertial rotation or zero central
rotational kinetic energy. With the nonzero C used here,

```text
K_C = ½I(Ω + ω_s)² = ½IC² > 0
```

The physical-pendulum claim in B2 also remains unestablished. The sources give
the effective rotating-frame torque −Iα_b and name distributed gravity, but
do not derive the claimed material torque/work/potential/reaction account from
them. This gap concerns the causal claim, not the existence of the relative
rocking predicted by the coordinate equation. Treating ½Iω_s² alone as the
inertial central energy omits the baseline and cross terms. Constant central
kinetic energy does not rule out power transmission between opposing material
interactions; those individual interactions and powers remain unspecified for
the claimed lunar mechanism. The coordinate reconstruction does not establish
the stronger physical and energy interpretations in S1, S2, and B2.

Finally, this is a 28-day hindcast against JPL models, not a validation against
new raw observations or a long-term lunar theory. A constant star-relative
orientation rate and a fixed pole omit physical librations and precession.
The computation does not test the
18.6-year Cassini-state claim in S1 or justify extending the reported monthly
accuracy beyond the tested interval.
