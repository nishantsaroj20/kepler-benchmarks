# Requirements inventory

This is a checklist for a later human-written `instruction.md`. It is not that instruction. No agent-facing prose is authored here.

Run used for every number: `pilots/converter-robust-control/runs/2026-11-08-stage0-4`. Source `595ff19f1b24955dfbf15423ab5f889342f49fa6`.

## Identity

| Item | Value |
|---|---|
| Future task name | `afterquery/converter-robust-control` |
| Domain / field / subcategory | `engineering-sciences` / `electrical-engineering` / Electrical & Computer Engineering |
| Plant set | `crc-two-gfl-lcl-v1` |
| Schema | `1.0.0` |
| Future artifact | `/app/output/controller.json` |
| Agent timeout, if a bundle is ever built | 28800 s, from the guide, not from this pilot |
| Verifier network | no-network, separate image, if a bundle is ever built |
| Harbor, nop, oracle, probes | UNVERIFIED. Not run. |
| Current submission-screen difficulty band | UNVERIFIED. The supplied PDF says both "at least 0 and at most 7" and that zero solves are treated as unsolvable. Swift's old band is not assumed. |

## Public cases

`nominal`, `weak_grid`, `high_power`, `asymmetric`. Same controller. Case data and equilibria are in the plant JSON files. Biases `u_bias_pu` are applied outside the controller.

## Controller contract

Keys, dimensions, signs, size, symlink rule, and `σ_max(D) ≤ 15`: see the specification, section 8. Realizations are not unique. `Ak` need not be Hurwitz.

## Pass rules that were actually enforced

1. Physical closed loop: `abscissa + uncertainty ≤ -0.5`.
2. Weighted closed loop clearly Hurwitz: `abscissa + uncertainty ≤ -0.001`.
3. SLICOT `AB13DD` value and the refined frequency lower bound agree within `max(0.03, 0.01 · scale)` and their maximum is `≤ 1`.
4. `σ_max(D) ≤ 15`.
5. For the witness flag only: both norm estimates `≤ 0.95`, physical `abscissa + uncertainty ≤ -1`, and the bounded-real lemma certifies `γ = 1` with a negative recomputed LMI eigenvalue.

A finite SLICOT number on an unstable loop is not a pass. The catalog system `1/(s-1)` returns L-infinity 1 and is recorded as infinite H-infinity.

## Traceability

| Requirement | Trusted data | Check in this pilot | Passing example | Failing example |
|---|---|---|---|---|
| Schema and block structure | `schema.py` | `load_controller` on `hand_pi_controller.json` | reloaded PI | not re-tested as a full adversarial suite; that suite belongs to a later verifier |
| `D_yu = 0` | plant JSON | max entry is 0 in all four files | all cases | none found |
| Two model paths | equilibria | `path_comparison.json` | relative Jacobian error ~1e-11 | none found |
| Physical stability margin | scaled physical `A, Bu, Cy` | `physical_closed_loop` | PI `kp=0.4`, `ki=4`, abscissa about -2.9 to -3.0 | zero law realized with `Ak = 0` (poles at 0); wrong-sign PI `kp=ki=-0.4/-4` |
| Stacked weighted H-infinity ≤ 1 | generalized plants | SLICOT plus refined peak | same PI, values 0.305 to 0.319 | zero law, SLICOT about 10.7 to 11.6; static gains, best about 6.2 |
| Gain bound | `D` only | singular value | PI `D` singular value 0.4 | not the binding failure in this run |
| Same controller on every case | all four plants | `pass_all` | tied PI grid, many entries | a controller is not accepted on one case alone |
| Bounded-real certificate at 1 | Clarabel residual | witness flag | the same PI, LMI eigenvalue about -0.07 to -0.08 | zero law, not certified; narrow catalog resonance not certified at 1.02 times its peak |
| Rounding neighborhood | gains of the hand PI | 12 perturbations at 1e-6, 1e-4, 1e-3 relative | all 12 still pass, worst upper stays 0.319 | none in that neighborhood |

## Explicit non-requirements

No `instruction.md`. No `task.toml`. No Docker image. No oracle script. No sealed tests. No claim that the four cases cover a continuum or a nonlinear transient.
