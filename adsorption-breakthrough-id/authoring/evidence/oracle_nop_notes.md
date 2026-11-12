# Oracle / nop / fairness notes (Anwes tips)

## Visible condition coverage

Every deciding condition used in hidden evaluation appears at least once in visible samples:

- `V1_base` — clean baseline
- `V2_dense_pack` — dense packing class (also `H1`)
- `V3_dispersion_regime` — high-Pe / large `D_ax` (also `H2`)
- `V4_sensor_faults` — mid-bed NaN/spikes (also `H3`)
- `V5_high_cin` — elevated feed (also `H4`)
- `V6_near_linear` — small-K isotherm (also `H5`)

Hidden cases reuse those classes with different numbers (not reseeds; no visible case ID reused).

## Labels

Velocity field is `interstitial_velocity_m_s` = PDE advection speed $u$ (not a mislabeled superficial velocity).

## Standard-method grading

`L_min` is graded by true-model outlet feasibility (+ oversize ceiling), not by matching one proprietary binary-search length. Parameter tolerances were checked with the reference solver and an independent alternate-grid feasibility check (`independent_calibration.py`).

## Oracle / nop

- Oracle: pytest suite passes with `solution/column_solver.py` installed as the artifact.
- Nop: missing `/workspace/output/column_solver.py` → reward 0.
- Timeout: hard 120 s per hidden case fails the task.
