# Oracle / nop / naive checks

Host-side validation (Docker daemon unavailable in this environment):

- Oracle (`solution/solver.py` → `output/adsorption_identification.json`): pytest `tests/test_evaluation.py` → **4/4 pass** (reward 1).
- Nop (missing artifact): schema test fails → reward 0.
- Naive baselines that fail at least one gate:
  - constant isotherm / crude transport guess
  - copying `lab_bed_length_m` as `L_min_m`
  - near-zero `D_ax`
  - NaN payload

Independent calibration (`independent_calibration.py`) re-integrates the true model on a separate MOL grid, confirms the reference `L_min` sits on the purity cliff, and shows a 8% undersized bed exceeds the outlet limit.
