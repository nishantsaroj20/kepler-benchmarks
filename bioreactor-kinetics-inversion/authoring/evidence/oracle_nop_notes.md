# Oracle / nop / baseline evidence

Local Python validation (Docker daemon unavailable in this environment; images/Dockerfiles are authored for Harbor).

## Oracle

- Outputs from the reference solver path under `output/`.
- `pytest tests/test_evaluation.py`: **5 passed** → score **1**.

Representative recovered parameters vs ground truth:

| key | estimate | truth | error |
| --- | --- | --- | --- |
| mu_max | 0.49995 | 0.50 | 0.01% rel |
| K_S | 0.750 | 0.80 | 6.2% rel |
| K_I | 18.32 | 22.0 | 16.7% rel |
| P_crit | 32.90 | 30.0 | 9.7% rel |
| t_switch | 11.79 | 12.0 | 0.21 h abs |
| alpha | 0.4966 | 0.50 | 0.68% rel |
| beta | 0.0984 | 0.10 | 1.6% rel |

## Nop

- Empty `output/` (missing JSON artifacts): failures → score **0**.

## Naive baselines

- Zero feed + far-off parameters: fails parameter gates and product-mass gate.
- NaN parameter values: rejected by finite-value checks before tolerances.
- Constant feed 0.12 L/h: fails mass gate (mass ≈ 154 g vs min 249.33 g).
- Independent calibration (`independent_calibration.py`): reference feed mass 277.03 g (≥ 249.33); constant-feed mass 154.24 g fails; CSV-mean fake kinetics fail primary relative-error gates.

## Independent reimplementation

`authoring/evidence/independent_calibration.py` re-integrates the true Radau model without importing `solution/solver.py`, confirming reference feasibility and baseline failure modes used in `verification_explanation`.
