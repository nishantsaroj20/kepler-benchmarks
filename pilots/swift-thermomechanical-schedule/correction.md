# Corrected pilot

This file supersedes the PASS table in `individualreport.md` and the augmented `results/gate_verdicts.json`. Those files are unchanged on purpose. The pre-correction proportional run is in `results/precorrection_immutable/proportional_family.json`. The corrected solver’s fresh outputs are in `results/corrected/`.

## Verdict

The one-variable no-hold family passes. A moderate cycle time, with stage durations in proportion to the three displacement increments and with no hold, meets every published limit on every scenario. The four-variable schedule was not optimized, and this pilot does not show that optimizing the four durations is necessary. M2 and M3 in the first report were hand-picked.

Do not lower the 8 s cap or tighten the temperature, gauge, noise, or mesh tolerances to push this family out. Conduction is physically real: the same warm 4 s schedule is 340.2 K with `k = 45` and 376.1 K with `k = 0` in the same integrator, so only the conducting case passes 345 K. That is a material effect. It is not a four-dimensional design problem. The practical choice is the duration.

Frontier difficulty is **UNVERIFIED**. No model probe was run. Recommendation: reconsider this redesign as a way to make the task hard. A later instruction should not be written as if the agent must search `(t1, t2, t3, h)`.

## 1. Proportional family, before and after the solver fix

Family, not a 4-D optimum: `(t1, t2, t3, h) = (0.35, 0.375, 0.275, 0) * tau` for `tau = 4, 5, 6, 7` s, so

`(1.4, 1.5, 1.1, 0)`, `(1.75, 1.875, 1.375, 0)`, `(2.1, 2.25, 1.65, 0)`, `(2.45, 2.625, 1.925, 0)`.

Pre-correction solver, 4×8, 8 substeps, `k = 45`, stagger cap 2, which could still commit an unconverged stagger. Command: `python3 reproduce_precorrection_proportional.py` from `src/`. Twelve of twelve passed, zero rejected steps, 20.5 s. Warm `tau = 4` neck maximum 340.24 K. That matches the independent user run (about 48 s there, same temperature).

Corrected solver, same mesh and substeps, but a step is kept only when Newton meets its residual tolerance and the stagger change is at most 0.25 K. Rejected attempts restore displacement, plastic strain, the elastic Finger tensor, the deformation gradient, and the nodal temperature. The base labels did not change: twelve of twelve still pass, zero rejections. Warm `tau = 4` is again 340.24 K, gauge 2.7493 mm, minimum radius 2.3356 mm, peak force 11541 N.

Margins on that warm 4 s case, which is the tightest temperature: 4.76 K under 345 K, 2.60 s under 8 s, 0.136 mm over the 2.20 mm radius floor, about 2460 N under 14000 N, and 0.349 mm over the 2.40 mm gauge floor. Full accepted-step histories, including the t = 0 row, are `results/corrected/histories/base_P4_*.npz` and the siblings for P5–P7 and the hand-picked M2.

For warm `tau = 6` and `tau = 7` the cycle maximum is 333.00 K at t = 0, because the initial temperature is now included. The maximum after t = 0 is lower. Both still pass. The first report’s 332.6 K and 332.4 K were the process maxima with t = 0 omitted.

## 2. What the solver fix changed, and what it did not

`src/test_acceptance.py` checks three things. Restoring a snapshot is bitwise on displacement, plastic strain, Finger tensor, deformation gradient, nodal temperature, and Gauss-point temperature. A Newton patched to fail after corrupting those arrays does not leave the corruption in place once the caller restores. A clean 4 s run records `coupling_converged` on every accepted step and starts from neck temperature 333 K.

The fix does not make the hold at `u_end = 1.45 mm` robust. These corrected runs still stop there, and the step is not committed:

- M2, nominal, 4 substeps, `k = 45`
- M2, all three scenarios, 8 substeps, `k = 0`

Finer time steps (16 substeps) pass M2 at `k = 45`. The `k = 0` hold failure was not retuned.

## 3. Refinement of the simple family

P4, P7, and hand-picked M2, all three scenarios, unless noted. Conductivity 45 unless noted.

Warm P4 neck maximum, the number nearest the cap:

| Discretization | Neck max (K) | Gauge (mm) | Min radius (mm) | Pass |
|---|---:|---:|---:|---|
| 3×6, 8 substeps | 339.39 | 2.718 | 2.361 | yes |
| 4×8, 4 substeps | 340.99 | 2.811 | 2.293 | yes |
| 4×8, 8 substeps | 340.24 | 2.749 | 2.336 | yes |
| 4×8, 16 substeps | 339.90 | 2.727 | 2.352 | yes |
| 5×12, 8 substeps | 340.06 | 2.777 | 2.312 | yes |
| 5×12, 16 substeps | 339.75 | 2.754 | 2.331 | yes |
| stagger tolerance 0.50 K and 0.05 K | 340.24–340.25 | 2.749–2.751 | 2.335–2.336 | yes |

The temperature spread on this one schedule is 1.24 K. The gauge spread is 0.084 mm. The radius spread is 0.068 mm. None of those moves crosses a limit. P7’s warm cycle maximum stays at the initial 333 K on every completed primary grid. M2’s warm temperature stays 334.2–335.4 K on the completed grids. The one missing M2 row is the failed 4-substep nominal hold, not a silent pass.

## 4. Conduction, same integrator

Only `k` changes. The heat operator, end Dirichlet condition (`T = 293 K` on `Z = L` for every `k`, including zero), mesh, substeps, card, and stagger tolerance stay the same. At `k = 0` there is no conductive flux; plastic heat is still projected onto element nodes, and the share assigned to grip nodes leaves the free-node balance. That is the discrete zero-conductivity limit of this operator, not the separate Gauss-point adiabatic return used in the first report.

| Case | k = 45 neck (K) | k = 0 neck (K) | k = 45 pass | k = 0 pass |
|---|---:|---:|---|---|
| P4 warm, primary 4×8 / 8 | 340.24 | 376.09 | yes | no |
| P7 warm, primary 4×8 / 8 | 333.00 | 378.59 | yes | no |
| P4 warm, independent 4×8 / 4 | 339.36 | 371.31 | yes | no |
| P4 warm, independent on the primary’s nodes, 4 substeps | 340.04 | 374.65 | yes | no |
| M2 warm, independent 3×6 / 4 | 334.24 | 368.92 | yes | no |

Primary M2 at `k = 0` did not finish. The independent Piola solver did, and it shows the same warm failure. So the conduction effect on the warm moderate schedule is not an artifact of one residual. The first report’s adiabatic-versus-coupled comparison changed the integrator as well as `k`. Those files are still in `results/precorrection_immutable/`. They are not the controlled comparison.

## 5. The old 6.1 K and 0.174 mm gap

That gap was adiabatic, nominal, 3×6, four substeps, two different meshes: primary gauge 2.455 mm, reference gauge 2.379 mm, temperatures 336.2 K and 331.8 K. Those two gauges do sit on opposite sides of 2.40 mm. Force agreed to 0.06%, which does not certify the gauge.

On the conducting warm 4 s problem, with the reference solver given the primary node coordinates, the gauges are 2.712 mm (reference, 4 substeps) and 2.749 mm (primary, 8 substeps). Both are above 2.40 mm. Neck temperatures are 340.04 K and 340.24 K. Peak forces are 11557 N and 11541 N. The large gap tracked the coarse mismatched mesh and the separate adiabatic return, not a remaining constitutive bias on the decision mesh. The coarse independent mesh can still put a `k = 0` nominal gauge at 2.377 mm, under the floor. That case is not the conducting design point, and it is reported rather than adjusted.

## 6. Calibration

Observations were regenerated with the corrected solver. Seed 20261108. Generating grid 5×12 / 8. Fit grid 4×8 / 6. All five channels enter the residual. The denominator for each channel is the hypotenuse of the noise standard deviation and the RMS of the fit-grid versus generating-grid difference.

That RMS is not uniform. Slow-test force differs by 2467 N at the single sample `t = 1 s`, where a coarse accepted step straddles the elastic rise, and by at most 92 N at the other nine samples. Temperature discrepancies are 0.09–0.23 K, inside the 0.4 K noise. Gauge discrepancies are several times the 0.002 mm noise. The standardized RMS of a fit is therefore not a statement that force was fit to 20 N.

Four starts were attempted. Two returned real residuals:

| Start | Standardized RMS | K | alpha | C | m |
|---|---:|---:|---:|---:|---:|
| 0, from the generating card | 0.743 | 683.5 | 0.408 | 0.0135 | 0.784 |
| 2, from K=900, alpha=0.25 | 0.746 | 743.8 | 0.330 | 0.0137 | 0.777 |
| Generating card |  | 680 | 0.50 | 0.012 | 0.80 |

`alpha` and `K` move. `C` and `m` stay close. On P4, all three scenarios, both fitted cards and the generating card pass. Warm neck maxima are 340.23 K, 340.24 K, and 340.24 K.

Start 1 and start 3 never produced a real residual: every evaluation was the solver-failure penalty (standardized RMS 20). Replaying the start-3 parameter vector, which is not a fit, violates the gauge ceiling and the radius floor. That is a bad card failing the design, not evidence that accepted cards disagree. A full profile or bootstrap was not computed. The claim that is supported is narrower: the two fits that actually match the data, and the generating card, do not change the simple-family pass.

## 7. Corrected gates

| Gate | Verdict | Replaces |
|---|---|---|
| Pre-correction proportional reproduction, 12/12, warm 4 s at 340.24 K | PASS | new |
| Corrected acceptance rule; rollback tests; base labels unchanged | PASS | the old commit-on-failed-stagger |
| Simple no-hold family passes on the base mesh and on the P4 refinement set | PASS | the claim that simple baselines fail |
| Four-variable schedule hardening is required | FAIL | the old “robust feasible region / meaningful optimization” PASS |
| Hand-picked M2 still passes where the solve finishes | PASS | keep, but label it manual |
| Same-integrator `k` comparison flips warm P4 and warm P7 | PASS | the old mixed adiabatic/coupled claim |
| Primary M2 at `k = 0` | FAIL | Newton stops at the hold; not retuned |
| Independent conducting solver agrees on the warm flip | PASS | the old independence PASS that rested on force |
| Shared-node gauge/temperature gap is a few hundredths of a millimetre and 0.2 K | PASS | the unresolved 0.174 mm adiabatic gap, which remains a coarse-mesh result |
| Noise-consistent local fits do not flip P4 | PASS | the old two-card claim, now with the RMS and the failed starts stated |
| Full uncertainty ensemble | UNVERIFIED | two starts never evaluated a real residual |
| Continuous inter-step peak certificate beyond accepted steps and the refinement table | UNVERIFIED | extrema are the accepted steps, including t = 0 |
| Frontier difficulty | UNVERIFIED | unchanged |
| Overall: use this table, not the first report’s PASS list | SUPERSEDES | `individualreport.md`, `results/gate_verdicts.json` |

## Commands

From `pilots/swift-thermomechanical-schedule/src`, with one thread:

```text
python3 reproduce_precorrection_proportional.py
python3 test_acceptance.py
python3 run_corrected_pilot.py
```

The third command died in the joint fit when a failed residual had the wrong length (`EXIT:1` in `/tmp/corrected_stdout.txt`). Histories from that run were already on disk. The finish command, after the residual-length fix, was:

```text
FIT_STARTS=0,2,3 python3 finish_corrected.py
```

That exited 0. Python 3.12.3, NumPy 2.1.3, SciPy 1.14.1. Calibration seed 20261108. It re-read the saved histories and re-ran calibration; it did not trust the old `results/*.json` cache.
