# Pilot report

**Verdict: STOP.**

A tied diagonal PI meets every frozen requirement on all four operating cases. Several neighboring gains do too. Zero control and static proportional control do not. Full-order H-infinity does better still, so the specification is feasible, but feasibility is not in question and the decentralized order limit is not what makes the problem hard. Thresholds were not tightened. No Harbor bundle should be built on this specification.

Official run: `runs/2026-11-08-stage0-4`. Git `595ff19f1b24955dfbf15423ab5f889342f49fa6`. Wall clock 86.4 s. Python 3.12.3, numpy 2.2.6, scipy 1.16.2, slycot 0.7.1, cvxpy 1.7.3, clarabel 0.11.1. One thread requested. Machine `Linux x86_64`, 4 CPUs reported by `os.cpu_count`.

## Gates

| Gate | Result |
|---|---|
| 0 Application | Pass as an engineering description. Two grid-following LCL inverters, one local outer controller, four public operating points. |
| 1 Model | Pass. Equilibria, power balance, and the two Jacobians agree. `D_yu = 0` by the output map. |
| 2 Interface | Pass. Schema and closed-loop formula checked against a second interconnection. |
| 3 Evaluation | Pass for the controllers cited below. SLICOT and the refined peak agree. Bounded-real lemma certifies `γ < 1`. Hamiltonian bisection is reliable on the smooth catalog and on the cited PI, and is inactive when its probe misses a narrow peak. |
| 4 Pilot | STOP. Ordinary tied PI clears the frozen specification with margin. |

## Cited controller

Both blocks:

```
A = 0,  B = I,  C = diag(4, 4),  D = diag(0.4, 0.4)
```

`σ_max(D) = 0.4`. JSON reload of `hand_pi_controller.json` reproduces the pass.

| Case | Abscissa | Uncertainty | SLICOT | Refined peak | Hamiltonian upper | BRL at 1 |
|---|---|---|---|---|---|---|
| nominal | -2.893 | 5.1e-4 | 0.31755 | 0.31755 | 0.31750 | certified, LMI -0.0756 |
| weak_grid | -3.001 | 5.0e-4 | 0.30533 | 0.30533 | 0.30524 | certified, LMI -0.0801 |
| high_power | -2.974 | 5.1e-4 | 0.31943 | 0.31943 | 0.31930 | certified, LMI -0.0763 |
| asymmetric | -2.895 | 4.8e-4 | 0.31793 | 0.31793 | 0.31787 | certified, LMI -0.0726 |

Twelve relative perturbations of these gains, four each at `1e-6`, `1e-4`, and `1e-3`, all still pass. The worst replayed upper bound stays 0.3195.

A second tied PI, `kp = 0.1`, `ki = 16`, is the grid minimizer of the worst-case SLICOT value (0.1004). SLICOT and the refined peak match to about `1e-7`. The bounded-real lemma again certifies `γ = 1`. The Hamiltonian probe is inactive on that loop; the pass does not use an inactive probe as an upper bound. Physical abscissae are about -7.85.

The PI grid contains many other passes (`kp` from 0.1 to 0.8 and `ki` from 2 to 40, with the high-`kp` edge closer to the limit). This is a region, not a single coefficient vector.

## Families

| Family | Result |
|---|---|
| Zero law, stable realization | Stable (abscissa about -7.86) and fails performance. Worst SLICOT 11.61. Bounded-real test at `γ = 1` does not certify. |
| Zero law with `Ak = 0` | Fails internal stability. Poles sit at 0. |
| Static diagonal | Best stable gain near 0.80 still has SLICOT 6.16. Gains from about 2 upward lose stability on at least one case. |
| Wrong-sign PI | Unstable. Abscissa about +8.5. |
| Tied PI grid and two DE seeds | Pass. DE worst-case SLICOT about 0.090. |
| Cross-term PI, one long DE | Pass, worst-case SLICOT 0.086. Iteration limit reached; still a pass, not a better problem. |
| Nominal-only selection | The nominal winner `kp=0.1`, `ki=16` also passes the other three cases. |
| Per-case PI | Each case has a tied PI near 0.08 to 0.10. Sharing is not the obstacle. |
| Decoupled-model design | Best decoupled grid point replays on the coupled plants at 0.105 and passes. Interaction is visible in the model and is not binding at this specification. |
| 32-parameter Powell, 4 seeds, 25 evaluations | All four remain without a finite SLICOT value (stability not achieved). Too short to search that dimension. Not evidence of infeasibility. |
| `sb10ad` full order, inadmissible | Independent recomputed norms 0.071, 0.054, 0.113, 0.076 by case. Reported `sb10ad` gammas are 0.065, 0.055, 0.055, 0.053. The high-power pair differs (0.055 reported versus 0.113 recomputed). The recomputed number is the one used here. All are far below 1, and the physical loops pass. |
| State-feedback LMI at `γ = 1` | Feasible on every case. Optimistic only. |

## Evaluation limits, stated as limits

- The narrow catalog resonance, damping `1e-4` at 200 rad/s, is recovered by SLICOT as 5000.000025. The frequency probe that feeds the Hamiltonian test misses it (lower estimate 6.68) and the test is marked inactive. A coarse grid is not an H-infinity norm. The bounded-real solver did not certify 2% above that peak. On the converter loops the peaks are broad and the refined search matches SLICOT.
- `python-control` was not imported. PyPI returns a quarantined project page with no files. The `linfnorm` documentation was read.
- MATLAB `hinfstruct` was not run.
- Harbor oracle, nop, quality check, and any agent probe are UNVERIFIED. Not applicable until a bundle exists, and this verdict says not to build one.
- The current new-task difficulty screen was not opened. The PDF's "0 to 7" sentence and its "solved never is unsolvable" sentence disagree. No band is assumed.

## What was not concluded

The plant is not "uncontrollable" and the specification is not "infeasible." Full-order design and ordinary PI both succeed. Adding cases, states, or tighter weights after this result would be changing the problem to chase difficulty. That was not done.

If a harder commissioning problem is wanted later, it needs a separate engineering reason: for example a resonance that the factory current loop does not already damp, or a requirement the outer PI cannot see. That would be a new specification and a new pilot, not a silent edit of this one.
