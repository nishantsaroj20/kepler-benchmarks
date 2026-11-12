# adsorption-breakthrough-id

Write a timed program that identifies Langmuir–LDF–dispersion parameters from fixed-bed breakthroughs and computes a minimum design packed length. The verifier re-runs the program on hidden, deliberately shifted cases.

## Difficulty

Visible samples are synthetic Langmuir–LDF–dispersion breakthroughs with sparse noisy outlet and mid-bed ports. The submitted program is re-run under a 120 s / 1 CPU limit on hidden cases that deliberately shift material properties, dispersion regime, sensor integrity, feed range, and isotherm linearity. Difficulty comes from ill-conditioned `D_ax` versus `k_LDF` identification, robustness to mid-bed faults, and a design length that fails when parameters are only roughly right, all under a tight per-run compute budget.

## Reference solution

Seed Langmuir parameters from PDE-consistent completed-breakthrough uptake on both experiments, then multi-start log-space least-squares of all four unknowns against front-weighted outlet and cleaned mid-bed residuals using a method-of-lines BDF column model. Binary-search the minimum packed length for the case design duty with a small safety margin, keeping the MOL grid coarse enough to finish each case within the 120 s per-run limit on one CPU.

## Verification

The verifier executes `/workspace/output/column_solver.py` via subprocess on five sealed hidden cases under OMP/OpenBLAS thread pins with a hard 120 s timeout per case; timeout or nonzero exit yields reward 0. Outputs are checked for finite positive keys and relative tolerances (12% `qm`, 25% `K`, 50% `k_LDF`, 45% `D_ax`, 20% `L_min` with a 0.95 undersize floor). Hidden shifts are built in `authoring/provenance/generate_all_cases.py` and documented in `authoring/evidence/hidden_shift_notes.json`; no visible case appears in the hidden set.
