# adsorption-breakthrough-id

Write a timed program that identifies Langmuir–LDF–dispersion parameters from fixed-bed breakthroughs and computes a feasible design packed length. The verifier re-runs the program on held-out cases from the same condition classes shown in the visible samples.

## Difficulty

All breakthrough curves and design duties are **synthetic**: they are produced by a seeded method-of-lines generator (`authoring/provenance/generate_all_cases.py`) of a 1-D Langmuir–LDF–axial-dispersion column, then sparsified and noise-corrupted. The shapes are realistic for laboratory / pilot fixed-bed adsorber breakthrough tests (outlet and optional mid-bed ports). In practice this inverse sizing problem is handled by adsorption process engineers and separation scientists during adsorbent screening and packed-bed design. Difficulty is scientific: ill-conditioned `D_ax` versus `k_LDF` coupling, robust identification when mid-bed sensors fault, and a purity-constrained design length under a 120 s / 1 CPU program-under-test budget.

## Reference solution

Seed Langmuir parameters from PDE-consistent completed-breakthrough uptake on both experiments, then multi-start log-space least-squares of all four unknowns against front-weighted outlet residuals (and cleaned mid-bed points when available) using a method-of-lines BDF column model. Binary-search a feasible packed length for the case design duty with a small safety margin, keeping the MOL grid coarse enough to finish each case within the 120 s per-run limit on one CPU.

## Verification

The verifier copies each sealed case into a temporary workdir and runs `/workspace/output/column_solver.py` there under thread pins, `python -s`, a new process group (always `killpg` after completion), and — when the verifier is root — an unprivileged `nobody` user with **no silent root fallback**. Hard timeout is 120 s per case (no network). Parameter relative tolerances are 12% `qm`, 25% `K`, 50% `k_LDF`, 45% `D_ax`. Reward is derived from CTRF/JUnit passed/failed/skipped counts (expect 10 passed, 0 failed, 0 skipped), not pytest exit code alone.

`L_min` gates use two calibrated constants (see `authoring/evidence/independent_calibration_report.json`):

- **1.30 `L_min` oversize ceiling** — reported `L_min` ≤ 1.30 × reference design length (reject padded beds). Measured max oracle `L_est/L_ref` ≈ **1.049**.
- **1.03 outlet-feasibility slack** — on the true parameters, max outlet at reported `L_min` ≤ 1.03 × purity limit under the **same BDF / N=64** check the verifier uses (plus an n=80 cross-check). Measured max outlet/limit ≈ **0.762**.

Visible-parameter cheats (copy any visible sample’s params at L×1.0/1.1/1.2) fail all five hidden cases under those gates (`authoring/evidence/visible_cheat_discrimination.json`).

An independent alternate-structure estimator (`authoring/evidence/independent_estimator.py`: BDF N=40 equal-weight outlet+mid, log10 least_squares — not the oracle’s N=32 front-weighted ln-space fit) was executed on the **actual verifier path** (isolated subprocess + BDF/N=64 gates + 120 s). All five hidden cases passed. Measured relative errors and wall times from `authoring/evidence/independent_estimator_report.json`: H1 (`qm` 0.009, `K` 0.021, `k_LDF` 0.045, `D_ax` 0.059, **36.2 s**), H2 (0.011, 0.008, 0.065, 0.092, **19.6 s**), H3 (0.006, 0.025, 0.239, 0.094, **32.7 s**), H4 (0.001, 0.021, 0.006, 0.164, **42.3 s**), H5 (0.074, 0.089, 0.192, 0.220, **24.1 s**).
