# Executable pilot plan

This is the plan that was encoded in `params.py` and executed by `scripts/run_fresh.py`. The frozen JSON was written at the start of the run, before `run_campaign`. Thresholds were not changed afterward.

## Command

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python scripts/run_fresh.py --output runs/2026-11-08-stage0-4
```

The script exits if the output directory already exists. It does not read another run directory. Official output is `runs/2026-11-08-stage0-4` on commit `595ff19f1b24955dfbf15423ab5f889342f49fa6`.

## Before any controller search

1. Catalog checks for SLICOT, the Hamiltonian bracket, and a bounded-real certificate: `1/(s+1)`, feedthrough 2.5, second-order peak `1/(2ζ√(1-ζ²))`, cancellation to `1/(s+2)`, narrow resonance of damping `1e-4`, and unstable `1/(s-1)` kept as L-infinity only.
2. Build four equilibria and both Jacobians. Abort the scientific claim if the paths disagree beyond the recorded tolerances. They did not.
3. Write `FROZEN_REQUIREMENTS.json`.

## Search budget that was actually run

Threads pinned to 1. Wall clock of the official run: 86.4 s, of which the campaign was 43.0 s.

| Family | What was run |
|---|---|
| Zero law | `Dk = 0` with unused poles at -100, and the same law with `Ak = 0` |
| Static diagonal | 12 positive magnitudes from 0.05 to 8, plus -0.4, -1, and 0 |
| Tied diagonal PI grid | `kp ∈ {0.1,0.2,0.4,0.8,1.5,3}`, `ki ∈ {0.5,1,2,4,8,16,40}`, both blocks equal |
| Wrong-sign PI | `kp = -0.4`, `ki = -4` |
| Tied PI differential evolution | 2 seeds, 1101 and 1102, about 500 evaluations each |
| Tied PI with full 2-by-2 `C` and `D` | 1 seed, 1201, 1947 evaluations, iteration cap hit |
| General 32-parameter blocks | 4 Powell starts, seeds 1301–1304, 25 evaluations each, uniform initial guess in `[-5, 5]` |
| Nominal-only | same PI grid scored on `nominal` only, then replayed on all cases |
| Per-case PI | same grid scored separately on each case |
| Decoupled design | cross-converter `A` and `Bu` blocks zeroed, PI grid on that model, best replayed on the coupled plants |
| Centralized full order | SLICOT `sb10ad` per case, then the returned controller recomputed with this repository's closed-loop formula and `AB13DD` |
| State-feedback LMI | feasibility at `γ = 1` per case, Clarabel. Optimistic bound only |

Search cost for the grids and optimizers was SLICOT plus the physical stability gate. Every best simple controller cited in the verdict was replayed with the refined frequency bound, the Hamiltonian report, the formula-versus-interconnection check, JSON reload, and the bounded-real certificate.

The 25-evaluation Powell runs are a weak global search in 32 parameters. They are recorded as failures to find a stabilizer. They are not a proof that the 2-state family is empty. The PI grid already lies inside that family.

## Decision rule, fixed before the search

- STOP if a zero law, a static diagonal law, or a tied diagonal PI meets the pass rule on every case under SLICOT and the refined peak, with the bounded-real certificate at `γ = 1` for the cited witness.
- Do not tighten `γ`, the stability floor, the weights, or the case list after that happens.
- PROCEED only if that simple family misses and a shared admissible controller still passes both paths with margin. That branch was not taken.
- RECONSIDER if the simple witness cannot be cross-checked, or if even full-order synthesis cannot meet `γ ≤ 1`.
- Do not call the specification infeasible only because a short random search failed.
- MATLAB `hinfstruct` stays UNVERIFIED if it cannot be run.
- Harbor, the future verifier, and the submission-screen difficulty band stay UNVERIFIED.

## Numerical allowance

Agreement tolerance `max(0.03, 0.01 · scale)` is the spread allowed between SLICOT and the refined peak. It is separate from the physical pass level of 1. The witness levels 0.95 and -1 rad/s are comfort margins, not replacement specifications.
