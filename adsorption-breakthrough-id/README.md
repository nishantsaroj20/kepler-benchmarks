# adsorption-breakthrough-id

Write a timed program that identifies Langmuir–LDF–dispersion parameters from fixed-bed breakthroughs and computes a feasible design packed length. The verifier re-runs the program on held-out cases from the same condition classes shown in the visible samples.

## Difficulty

Visible samples already expose every deciding condition class (dense packing, high-Pe dispersion, mid-bed faults, elevated feed, near-linear isotherm). Held-out cases reuse those classes with different numbers. Difficulty is scientific: ill-conditioned `D_ax` versus `k_LDF` coupling, robust identification from sparse noisy breakthroughs, and a purity-constrained design length under a 120 s / 1 CPU program-under-test budget—not from hidden condition types or unwritten formulas.

## Reference solution

Seed Langmuir parameters from PDE-consistent completed-breakthrough uptake on both experiments, then multi-start log-space least-squares of all four unknowns against front-weighted outlet residuals (and cleaned mid-bed points when available) using a method-of-lines BDF column model. Binary-search a feasible packed length for the case design duty with a small safety margin, keeping the MOL grid coarse enough to finish each case within the 120 s per-run limit on one CPU.

## Verification

The verifier executes `/workspace/output/column_solver.py` via subprocess on five sealed held-out cases under OMP/OpenBLAS thread pins with a hard 120 s timeout per case; timeout or nonzero exit yields reward 0. Parameter relative tolerances (12% `qm`, 25% `K`, 50% `k_LDF`, 45% `D_ax`) are checked against provenance truth. `L_min` is accepted by true-model outlet feasibility with a modest slack and an oversize ceiling, so a correct standard method need not match one proprietary binary-search length. Independent calibration notes are in `authoring/evidence/`.
