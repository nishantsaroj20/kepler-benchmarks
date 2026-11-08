# Stage 1–3 report: thermomechanical schedule pilot

> Superseded. The PASS table in this file overstated what the runs showed. The corrected record is `correction.md` and `results/corrected/`. This file is kept so the earlier claims are not erased.



The v1 notch-forecast task was solved by the work it asked for. A redesign that still only asks for an inverse fit will not add difficulty. This pilot changed the decision to a four-variable loading schedule under coupled force, temperature, gauge and time constraints, then ran that decision. Two schedules pass every public scenario with margin. The fastest ramp, the fastest proportional ramp, and a short hold at full actuator speed do not. Conduction is what makes the moderate schedule legal on a warm start: the same durations are 334.5 K with conduction and 373.8 K if the bar is adiabatic.

Frontier difficulty is **UNVERIFIED**. No Harbor run and no model probe were started. Stages 4–8 were not started. No `instruction.md` and no submission ZIP were written.

## Diagnosis

The supplied low-effort screen (claude-opus-5-5, low, 1 of 3) does not explain the xhigh result. The same model identifier at xhigh solved 8 of 8 and missed the displayed 0–5 band. Trials 1–7, as pasted, are independent Numba or JAX axisymmetric F-bar solvers of multiplicative J2 plasticity with Hencky elasticity and an eight-parameter Swift–Voce rate/temperature law, plus mesh and time checks. Reported physical errors are about 0.13–0.4%. Runtimes were 62, 286, 88, 130, 100, 95, 99 and 144 minutes, all inside eight hours. Trial 8’s narrative is truncated; its method is not inferred. The original ZIP hash is `0c154a40d1f6d594dab71b35147d6abf3fb66a69589d6593c2d491d4b95fc43d` (46 files). Nothing in that ZIP was modified.

Tighter notches, more of the same printed fields, or a ban on acceleration would grade the solver the agents already wrote. The additional decision has to be a process choice that a good card does not finish by itself.

## What was run

Primary solver: `src/coupled_current.py`. Current-configuration F-bar Q4, incremental elastic Finger tensor, closed-form in-plane eigenpair, staggered lumped backward-Euler conduction. Independent cross-check: `src/reference_pk1.py`. Total Lagrangian, `Cp^{-1}` state, `numpy.linalg.eigh`, Piola residual, bisection return, uniform mesh. It does not import the primary solver.

Command, from `pilots/swift-thermomechanical-schedule/src`, with `OMP_NUM_THREADS=1`:

```text
python3 run_pilot.py
```

Logged in `logs/pilot.log`. Python 3.12.3, NumPy 2.1.3, SciPy 1.14.1, 4 CPUs. Calibration seed `20261108`. Public calibration CSV SHA-256 `d43146d85f1789f1db43065ef4e9bcc5d30ee685298700f0f59b4e26b6baf12b`. Results CSV SHA-256 `af9861d976ee2a9f14b680cd5511bb1b9c2526595559448b2c2a41f02d2d2b4c`.

Checkpoint wall times: calibration generation 11 s, fit 224 s, schedule sweep 46 s, mesh/time refinement 87 s, perturbations 62 s, S0 thermal modes 36 s, independent cross-check 310 s, fitted-card design 71 s. The warm adiabatic supplement was an extra 28 s. A deeper 2.6 mm extension was tried first and dropped because the neck radius was still mesh-sensitive; those runs are not in the frozen table.

## Schedule and limits

Half bar, `R0 = 3 mm`, `L = 10 mm`, gauge material point at `Z = 4 mm`. Fixed end increments 0.70, 0.75 and 0.55 mm. Hold, if any, at the stage-2 boundary. Smoothstep ramps, so `v_max = 2δ/t` and `a_max = 2πδ/t²`. Actuator caps 8 mm/s and 40 mm/s² bind the floors at 0.34, 0.35 and 0.30 s. Inertia at those caps is a few pascals.

Scenarios, all public: S0 at 293 K with 0.015 mm imperfection, S1 at 333 K, S2 with 0.030 mm imperfection. Grip temperature 293 K. Lateral surface insulated.

Limits the pilot then froze, after seeing the coarse gaps and before calling refinement a pass:

| Constraint | Limit | What it did |
|---|---|---|
| Neck-band average temperature, whole cycle | ≤ 345 K | Rejects FAST, CR, a 4 s full-speed hold, and M1 on S1 |
| Cycle time | ≤ 8 s | Rejects SLOW (`τ = 14.5 s`) |
| Final gauge opening | 2.40–2.95 mm | Rejects FAST, CR and the 4 s hold on S0 (openings 2.29–2.33 mm) |
| Minimum radius | ≥ 2.20 mm | Did not bind |
| Force | ≤ 14000 N | Did not bind (observed peaks 11.5–12.1 kN) |

Uniform stretch of the 2 mm end motion would open the gauge by 1.60 mm. Computed openings are 2.3–2.7 mm, so the gauge constraint is a structural result.

## Schedules

Decision mesh 4×8, 8 substeps per ramp, truth card.

| Schedule | (t1, t2, t3, h) s | τ s | S0 T / r / gauge | S1 T / r / gauge | S2 T / r / gauge | All pass |
|---|---|---|---|---|---|---|
| FAST | 0.34, 0.35, 0.30, 0 | 0.99 | 326.8 / 2.512 / 2.332 | 368.3 / 2.397 / 2.636 | 329.5 / 2.444 / 2.464 | no |
| CR | 0.382, 0.409, 0.300, 0 | 1.09 | 326.8 / 2.509 / 2.338 | 367.7 / 2.390 / 2.650 | 329.4 / 2.441 / 2.472 | no |
| H4 | 0.34, 0.35, 0.30, 4 | 4.99 | 316.5 / 2.530 / 2.287 | 352.5 / 2.408 / 2.620 | Newton failed | no |
| M1 | 1.0, 1.1, 0.9, 0 | 3.00 | 321.6 / 2.453 / 2.467 | 347.9 / 2.343 / 2.738 | 322.8 / 2.391 / 2.581 | no |
| M2 | 1.5, 1.6, 1.3, 1 | 5.40 | 316.1 / 2.428 / 2.525 | 334.5 / 2.337 / 2.747 | 316.7 / 2.368 / 2.631 | yes |
| M3 | 2.5, 2.6, 2.0, 0 | 7.10 | 314.0 / 2.409 / 2.557 | 332.4 / 2.339 / 2.726 | 314.5 / 2.345 / 2.670 | yes |
| SLOW | 5.0, 5.5, 4.0, 0 | 14.50 | 306.5 / 2.426 / 2.509 | 330.5 / 2.394 / 2.590 | 306.7 / 2.354 / 2.640 | no |

M2 is comfortably inside (5 K, 1 s, 0.05 mm). M3 passes but sits 0.9 s under the time cap, so it is feasible and not marked comfortable. Temperature and cycle time are the two constraints that cut the set from opposite sides. Gauge cuts the fast nominal runs as well. Force and radius do not.

M2’s ten neighbors (each duration ×0.85 and ×1.15, clipped to the box, and all four durations ×0.9 and ×1.1) pass S0, S1 and S2. The accepted set is a region.

## Conduction

S0, same mesh, adiabatic return versus coupled conduction:

| Schedule | Coupled neck max | Adiabatic neck max | Difference |
|---|---|---|---|
| FAST | 326.83 K | 327.18 K | 0.35 K |
| H4 | 316.48 K | 327.18 K | 10.7 K |
| M2 | 316.07 K | 337.20 K | 21.1 K |
| SLOW | 306.47 K | 339.67 K | 33.2 K |

A one-second ramp does not care about conduction. A hold or a slow ramp does. Nodal `k = 0` on the same coupled code recovers the adiabatic neck temperature to 2.0 K (FAST) and 0.9 K (SLOW). The `k = 0` runs of H4 and M2 failed Newton at `u_end = 1.45 mm` and were left failed.

Warm start, executed afterwards and stored in `results/thermal_warm_supplement.json`:

| Schedule | Coupled | Adiabatic |
|---|---|---|
| M2 | 334.45 K, gauge 2.747 mm, passes | 373.82 K, gauge 2.409 mm, fails the 345 K cap |
| SLOW | 330.54 K | 375.81 K |
| H4 | 352.51 K | 363.95 K |
| FAST | 368.33 K | 363.95 K |

M2 is legal only with conduction. That is a feasibility change of 39 K, against a refinement scatter below 0.2 K.

## Calibration

Three noisy tests, seed 20261108, noise 20 N, 0.002 mm, 0.4 K, 0.003 mm. Generated on 5×12. Fit on 3×6 from two starts, public rows only.

| Card | K | ε0 | n | α | Q | β | C | m |
|---|---|---|---|---|---|---|---|---|
| Generating | 680 | 0.015 | 0.18 | 0.50 | 220 | 14 | 0.012 | 0.80 |
| Fit start 0 | 669.5 | 0.0218 | 0.187 | 0.606 | 200 | 18.6 | 0.0103 | 0.921 |
| Fit start 1 | 727.1 | 0.0189 | 0.209 | 0.464 | 186 | 18.9 | 0.0142 | 0.861 |

The constants are not the same card. On M2 / S1 the three cards give neck maxima 334.52, 335.02 and 335.17 K, and all three pass M2, fail FAST on S0 and S1, and fail SLOW. No pass/fail disagreement on a completed case. Stage-2 costs (66 and 94) are a few noise standard deviations on the coarse fit mesh, not a claim of unique identification.

## Checks

Material point, separate scalar return: plastic consistency error 0 MPa, `det Fp = 1`, objectivity error `3×10⁻¹¹` MPa. Small-strain uniaxial Kirchhoff stress is 0.10% under `E` times engineering strain, as a Hencky measure should be.

Conduction benchmark, quenched warm bar, negligible plasticity: neck average versus the 1D cosine series differs by at most 1.26 K (at 4 s). Heat balance on that run is `1.6×10⁻⁹` N·mm. Coupled plastic runs in the sweep report lumped heat imbalance at roundoff.

M2 refinement, base 4×8 / 8 substeps versus 5×12 and versus 16 substeps: temperature change ≤ 0.16 K, force ≤ 0.10%, radius ≤ 0.023 mm, gauge ≤ 0.027 mm, pass label unchanged. FAST and SLOW likewise keep their fail labels. FAST S0 gauge stays under 2.40 mm on all three grids (2.332, 2.352, 2.298).

Independent 3×6 cross-check: isothermal and adiabatic force agree to at most 0.058%. Adiabatic neck temperature differs by 6.14 K (FAST) and 4.39 K (M2); radius by up to 0.073 mm; gauge by up to 0.174 mm. Those gaps are larger than the primary solver’s refinement scatter. They do not flip M2’s coupled decision, whose warm margin is 10.5 K and whose gauge sits 0.12–0.20 mm inside the window, but they are the residual disagreement and they are on a coarse mesh. The second solver does not repeat the conduction matrix.

## Failures left as failures

Newton died at `u_end = 1.45 mm` on: H4 / S2 in the sweep; H4 / S1 on the 5×12 mesh; H4 and M2 with `k = 0`; and most fitted-card H4 reruns at 6 substeps. H4 is already illegal on the runs that finished (S0 gauge 2.287 mm, S1 temperature 352.5 K). The step size was not quietly reduced until those points passed.

## Gates

| Gate | Verdict |
|---|---|
| Robust feasible region (M2 comfortable, M3 also legal) | PASS |
| FAST, proportional CR, and 4 s full-speed hold are not legal | PASS |
| Temperature and cycle time both cut the set; gauge cuts the fast nominal runs; force and radius do not | PASS |
| Conduction changes M2 on S1 from pass to fail | PASS |
| Mesh and time refinement, separately, keep the labels | PASS |
| Independent mechanics; adiabatic temperature within 6.1 K on 3×6 | PASS, with the gap above |
| Two fitted cards agree with the generating card on pass/fail | PASS |
| Series conduction ≤ 1.26 K; material-point identities | PASS |
| Pilot wall time is minutes, not an unbounded search | PASS |
| Frontier solve rate | UNVERIFIED |

## Not done

No verifier, no oracle/nop, no quality probe, no final ZIP, no `instruction.md`. The requirements a human author would still have to write are in `REQUIREMENTS_INVENTORY.md`. Equations, units and the validated numbers are in `specification.md`.
