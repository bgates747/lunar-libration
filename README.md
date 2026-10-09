# lunar-libration

Two separate models of lunar motion and longitudinal libration.

- [`author-2d-animation/`](author-2d-animation/README.md): a two-dimensional
  Keplerian orbit and uniform lunar rotation, with a renderer for the resulting
  optical longitudinal libration.
- [`dremt-derived-model/`](dremt-derived-model/README.md): a three-dimensional
  Sun–Earth–Moon approximation expressed in a moving Earth–Moon frame and
  compared with NASA/JPL reference data over 28 days.

The [provisional findings pending further information](dremt-derived-model/REPORT-2026-10-09.md)
include the equations, exact attributed comments, a textual-analysis appendix,
and a dedicated section spelling out the information needed to test its
stronger physical claims, with references beside the relevant analysis.

Each model has its own setup and execution instructions. Neither model is a
spacecraft navigation ephemeris. The second model combines the relative-frame
coordinate description with standard Newtonian dynamics and empirical initial
conditions; successful reproduction of the reference trajectory does not
establish all physical claims associated with that coordinate description.

## Working directory and layout

This directory is the repository and clone root. Run documented commands here;
all project input/output paths are relative to this root or located from it.
Neither model needs the author's private workspace. The derived model separates
[source](dremt-derived-model/src/moon_model.py),
[data](dremt-derived-model/data/moon-model-reference.json),
[configuration](dremt-derived-model/config/moon-model-test-protocol.json),
[dependencies](dremt-derived-model/requirements/requirements.txt),
[animation documentation](dremt-derived-model/docs/ANIMATION.md) and tests.
Generated artifacts belong in each model's ignored `outputs/` directory.

## Licence

The original code and documentation are licensed under the [MIT License](LICENSE).
Attributed third-party excerpts, linked source documents, and NASA/JPL source
data retain their respective rights and source terms; the MIT licence does not
relicense them.
