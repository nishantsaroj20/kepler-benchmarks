# adsorption-breakthrough-id

Write a timed program that identifies Langmuir–LDF–dispersion parameters from fixed-bed breakthroughs and computes a feasible design packed length. The verifier re-runs the program on held-out cases from the same condition classes shown in the visible samples.

## Difficulty

All breakthrough curves and design duties are **synthetic**: they are produced by a seeded method-of-lines generator (`authoring/provenance/generate_all_cases.py`) of a 1-D Langmuir–LDF–axial-dispersion column, then sparsified and noise-corrupted. The shapes are realistic for laboratory / pilot fixed-bed adsorber breakthrough tests (outlet and optional mid-bed ports). In practice this inverse sizing problem is handled by adsorption process engineers and separation scientists during adsorbent screening and packed-bed design. Difficulty is scientific: ill-conditioned `D_ax` versus `k_LDF` coupling, robust identification when mid-bed sensors fault, and a purity-constrained design length under a 120 s / 1 CPU program-under-test budget.

## Reference solution

Seed Langmuir parameters from PDE-consistent completed-breakthrough uptake on both experiments, then multi-start log-space least-squares of all four unknowns against front-weighted outlet residuals (and cleaned mid-bed points when available) using a method-of-lines BDF column model. Binary-search a feasible packed length for the case design duty with a small safety margin, keeping the MOL grid coarse enough to finish each case within the 120 s per-run limit on one CPU.

## Verification

The verifier copies each sealed case into a temporary workdir and runs `/workspace/output/column_solver.py` there under thread pins, `python -s`, a new process group, and (when root) an unprivileged user so the child cannot read `/tests/hidden` or write `/logs/verifier`. Hard timeout is 120 s per case. Parameter relative tolerances are 12% `qm`, 25% `K`, 50% `k_LDF`, 45% `D_ax`; `L_min` must be feasible on the true model and not trivially oversized. An independent alternate-structure estimator (`authoring/evidence/independent_estimator.py`: LSODA N=48, equal-weight outlet+mid, log10 grid search — not the oracle BDF/front-weighted least_squares) recovered all five hidden cases within those gates. Measured relative errors from `authoring/evidence/independent_estimator_report.json`: H1 (`qm` 0.008, `K` 0.020, `k_LDF` 0.273, `D_ax` 0.000), H2 (0.000, 0.028, 0.078, 0.121), H3 (0.001, 0.013, 0.026, 0.079), H4 (0.003, 0.017, 0.000, 0.000), H5 (0.004, 0.004, 0.048, 0.073); every case also met true-model `L_min` feasibility.
