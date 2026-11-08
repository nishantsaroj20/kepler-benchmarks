# Reviewed corrected Swift pilot

## Verdict

**The simple no-hold proportional family passes the tested constraints. Four-variable schedule hardening is not demonstrated; reconsider this redesign as a difficulty mechanism.** This is stronger evidence than the first pilot, but not complete validation of all numerical, uncertainty and reproducibility requirements. Frontier difficulty remains **UNVERIFIED**. No limits should be tightened merely to exclude the successful simple schedules.

This review supersedes overbroad gate wording and arithmetic in `correction.md` and the supplied corrected gate JSON; it does not modify them. It is an evidence review of the completed pilot, not an additional simulation run. No simulation, scientific implementation, repository checkout/edit, verifier or ZIP creation occurred during this audit.

## Sources and preservation

Let **P** be `/workspace/swift-pilot-correction-review/extracted/pilots/swift-thermomechanical-schedule/`.

- Full model: `P/specification.md` plus `P/SPECIFICATION_AMENDMENT.md`. Equations, units, geometry and public limits remain those of the proposed model. `src/problem_spec.py` is byte-identical to the first pilot; the specification only adds a warning referring to the amendment. No numerical limit was tightened.
- Corrected execution: `P/logs/corrected_pilot.log`; `P/results/corrected/assembled_rows.json`, `summary.csv`, `histories/*.npz`; `src/coupled_current.py`, `reference_pk1.py`, `run_corrected_pilot.py`, `finish_corrected.py`, `test_acceptance.py`.
- Before-fix evidence: `P/results/precorrection_immutable/proportional_family.json` and its historical reproduction script.
- Original correction report, amendment, corrected gate JSON and summary CSV match the separately supplied assets byte-for-byte.
- Research tarball SHA-256: `360474f93a0fd0e53aa2c33447a7007bd0e68aabbbaba417156c6ac483089ce1`.

All extracted originals are preserved; their hashes are in `original_artifact_sha256.json`. Derived audit tables are outside P. Numerical checks here only reduced stored arrays and subtracted published limits; no model code was imported or executed.

## 1. Principal result: an easy one-variable family exists

The schedules are `(t1,t2,t3,h) = (.35 tau,.375 tau,.275 tau,0)` for tau in {4,5,6,7} s. These retain the specified three smoothstep ramps; they are not continuous constant-rate motion and are not four-dimensional optimized schedules.

**12/12 before-fix rows and 12/12 corrected base rows pass all three scenarios**, on 4×8 elements / 8 substeps per ramp, k=45, generating card. The corrected base family has zero logged rejected attempts. The old defect therefore did not cause these particular pass labels.

The unchanged limits are neck-band average T≤345 K, tau≤8 s, minimum outer radius≥2.20 mm, axial force≤14000 N and final gauge∈[2.40,2.95] mm. The averages are not pointwise peak temperatures. S0 is nominal at 293 K; S1 starts at 333 K; S2 has the larger imperfection at 293 K. There is no combined warm/large-imperfection scenario.

| Schedule | tau (s) | S0 Tmax (K) | S1 Tmax (K) | S2 Tmax (K) | Worst radius margin (mm) | Worst force margin (N) | Nearest gauge margin across scenarios (mm) |
|---|---:|---:|---:|---:|---:|---:|---:|
| P4 | 4 | 319.823 | 340.239 | 320.750 | .13563 | 2111.57 | .11781 |
| P5 | 5 | 317.790 | 333.770 | 318.518 | .13334 | 2116.98 | .14576 |
| P6 | 6 | 315.936 | 333.000 | 316.528 | .13410 | 2117.21 | .15953 |
| P7 | 7 | 314.262 | 333.000 | 314.779 | .13824 | 2112.09 | .16107 |

Cycle-time headrooms are respectively 4, 3, 2 and 1 s. The CSV `per_scenario_proportional_margins.csv` provides every output and all six margins for each of the 12 base cases, including both gauge bounds and the nearer bound, with the exact source history filename.

### Corrected arithmetic for warm P4

Warm P4 uses `(1.4,1.5,1.1,0)` s. Its stored results and margins are:

| Quantity | Output | Headroom |
|---|---:|---:|
| Neck Tmax | 340.238681 K | 4.761319 K |
| Cycle time | 4.000000 s | **4.000000 s**, not 2.60 s |
| Minimum radius | 2.335628 mm | .135628 mm |
| Peak force | 11540.668497 N | 2459.331503 N |
| Final gauge | 2.749260 mm | .349260 mm above floor; **.200740 mm below ceiling**, the nearer limit |

The original correction's 2.6 s margin belonged to M2, not P4. Its gauge-floor margin was arithmetically correct but was not the nearest gauge margin. P4 does not meet the optional 5 K “comfortable” temperature margin, though it passes the actual 345 K requirement. Several simple schedules passing does not establish a continuum of feasible tau values or global optimality, but it directly defeats the claim that a four-variable search is necessary. No schedule optimizer ran.

## 2. Acceptance fix: supported, with clearly limited regression coverage

Static inspection supports the intended fix (`coupled_current.py:517–536,583–628,697–714`): snapshots include u, ep, elastic Finger tensor, both stored deformation-gradient components, nodal T and Gauss-point T. A failed Newton or unconverged stagger returns no payload; the run loop restores the snapshot and predictor history before halving the attempted step. Only a successful payload is committed. Hitting the stagger cap without meeting its temperature tolerance rejects the attempt. Force is now evaluated before committing the incremental constitutive state, rather than reapplying the update after commit.

The stored stdout records the four acceptance tests passing. Their scope is snapshot bitwise restoration, injected Newton/stagger failures followed by **explicit caller restoration**, and a clean run recording acceptance. They do not independently exercise every corrupted-state failure through the complete run-loop retry path. The actual P7/S2 16-substep run records three rejected attempts and completes, providing additional real retry evidence; that is not an exhaustive fault-injection test.

This audit checked all **91 NPZ observable histories**: finite arrays, strictly increasing saved times, t=0 included, and all assembled extrema agree with their source arrays. For histories storing primary residual/coupling diagnostics, no saved stagger delta exceeds its filename tolerance. Largest saved deltas were .249771 K for .25 K tolerance, .040766 K for .05 K, and .273173 K for .50 K. The largest saved primary Newton residual divided by max(1, absolute reaction) is 9.95×10⁻⁷, consistent with the 10⁻⁶ mechanical criterion. These diagnostics are self-reported solver quantities, not an independent recomputation of equilibrium. In particular, `post_residual_N` aliases the accepted Newton residual rather than separately verifying mechanics at the newly updated heat field; coupling accuracy remains tolerance-controlled.

The histories now include the initial 333 K neck state: warm P6/P7 maxima occur at t=0. They store complete accepted observable histories, not every internal state tensor or every rejected attempt. Peak times are retained or recoverable, but continuous extrema between steps remain uncertified.

**Remaining solver failures are real:** M2/S0 at 4 ramp substeps and k=45 fails at the hold; all three primary M2 k=0 cases fail there. They are unresolved numerical failures, not evidence that the physical schedules are infeasible. Their absence from assembled successful rows must not silently count as a pass.

## 3. Refinement: substantially broader coverage, not an asymptotic certificate

Base=4×8/8; mesh-only variants=3×6/8 and 5×12/8; time-only=4×8/4 and 4×8/16.

- P4, P7 and M2 have mesh/time variants in **all three scenarios**, except the explicitly failed M2/S0 time case.
- Combined 5×12/16 covers P4 in all three scenarios; P7 and M2 combined refinement covers **warm only**.
- Coupling-tolerance .05/.50 K comparisons cover warm P4 and warm M2 only. P5/P6 have no corrected refinement set, and no corrected four-dimensional neighborhood was refined.
- There are 20 completed primary k=45 rows for P4, 16 for P7, and 17 for M2 in these base/refinement/coupling sets: all pass sampled constraints. Complete margins for these 53 rows are in `reviewed_refinement_margins.csv`. The missing M2 case remains documented separately here and in the gate JSON.

Warm P4 illustrates the independent discretization changes:

| Mesh / ramp substeps | Tmax (K) | Peak force (N) | Gauge (mm) | Minimum radius (mm) |
|---|---:|---:|---:|---:|
| 3×6 / 8 | 339.3866 | 11547.42 | 2.71815 | 2.36051 |
| 4×8 / 4 | 340.9941 | 11499.71 | 2.81137 | 2.29270 |
| 4×8 / 8 | 340.2387 | 11540.67 | 2.74926 | 2.33563 |
| 4×8 / 16 | 339.8979 | 11552.17 | 2.72725 | 2.35196 |
| 5×12 / 8 | 340.0639 | 11534.33 | 2.77673 | 2.31195 |
| 5×12 / 16 | 339.7513 | 11549.08 | 2.75389 | 2.33058 |

The CSV/history files are authoritative for full precision. Across the supplied warm P4 primary variants, temperature range is 339.3866–340.9941 K: **1.6075 K spread**, not 1.24 K. Gauge range is 2.71815–2.81137 mm: **.09322 mm**, not .084 mm; radius spread is .06781 mm. The coupled-tolerance repeats leave Tmax near 340.24–340.25 K and do not change labels.

Across all 20 P4 rows, the smallest headrooms are 4.0059 K, .09270 mm radius, 2109.40 N force and .08312 mm nearest gauge bound. Thus no supplied primary P4 refinement crosses a limit, including S2 and the combined cases. Across all 16 P7 rows the corresponding minima are 12 K, .10180 mm, 2107.75 N and .12010 mm.

This is credible **discretization-sensitivity/label-stability evidence**, including extrema at accepted internal times rather than sparse report times. It is not proof of asymptotic convergence: grids remain coarse, no convergence order/error bound is established, mesh/time changes sometimes oppose each other, and a continuous interstep extremum certificate is absent. The first pilot's robust-neighborhood uncertainty/refinement cross-product is still not supplied by these corrections.

## 4. Conduction: the controlled comparison is now present

The primary paired runs use the same nodal heat operator, 4×8/8 settings, card, schedule and stagger criterion; only k changes. Cold-grip Dirichlet T=293 K is retained even at k=0. This is a **discrete zero-conductivity nodal limit**, not the previous separate Gauss-point adiabatic update; heat projected to constrained grip nodes is not retained in the free-node balance.

| Warm case | k=45 Tmax (K) | k=0 Tmax (K) | Outcome |
|---|---:|---:|---|
| Primary P4 | 340.2387 | 376.0888 | Temperature pass → fail |
| Primary P7 | 333.0000 | 378.59 | Temperature pass → fail |
| Independent P4, own 4×8 / 4 | 339.3611 | 371.3080 | Temperature pass → fail |
| Independent P4, primary nodes / 4 | 340.0393 | 374.6527 | Temperature pass → fail; values logged |
| Independent M2, 3×6 / 4 | 334.24 | 368.92 | Temperature pass → fail |

For primary P4, the k-only difference is 35.8501 K, much larger than the tested temperature sensitivity. A material conduction-related feasibility effect is supported for these tested warm cases. No search for the best achievable cycle time with/without conduction was performed.

## 5. Independent conducting verification: scope and remaining gaps

The reference module now contains its own conducting heat assembly and stagger procedure (`reference_pk1.py:277–351,394–434`), using independent PK1 mechanics and no import of the primary solver. **It is not only a warm comparison:** 3×6/4 independent k=45 and k=0 histories cover P4 and M2 in all three scenarios. All six conducting cases pass the sampled force, temperature, radius, gauge and time constraints. The refined own-mesh and shared-node comparisons, however, cover only warm P4. P5/P6/P7 and fitted-card designs have no corresponding independent conducting coverage.

The shared-node warm k=45 check gives reference gauge 2.71250 mm, Tmax 340.0393 K and peak force 11556.85 N, versus primary base gauge 2.74926 mm, 340.2387 K and 11540.67 N. Both pass. **This does not establish that every old bias is resolved.** The reference used 4 requested substeps/ramp with adaptive subdivisions, .5 K coupling tolerance and a 10⁻⁵ residual factor; primary base used 8, .25 K and 10⁻⁶. The shared-node independent log reports 73 saved states, versus 25 primary base states. The reference own-mesh warm 4×8 case has 81 states. Mesh matching did not match integration histories.

Even comparing the primary's 4-substep variant rather than its 8-substep base gives about .099 mm gauge and .955 K temperature differences from the shared-node reference. The raw correction's conclusion that the old gap was not a remaining constitutive bias is therefore stronger than the controlled evidence permits. The present warm conducting label agreement is useful; complete bias attribution and independently converged all-scenario/all-constraint accuracy remain **UNVERIFIED**. Shared-node outputs are in the execution log but no shared-node NPZ histories are supplied.

Both solvers still derive boundary heat from their own assembled heat residual. No independent boundary-flux or external-work/elastic-energy/plastic-dissipation audit is supplied. Independent coupled energy verification remains **UNVERIFIED**.

## 6. Calibration: joint fits improve coverage, but .74 is not a noise-only certificate

The corrected calibration regenerates noisy records on 5×12/8 and jointly fits all eight parameters using all five observation channels on 4×8/6. This addresses the old partial-residual problem. Successful local fits 0 and 2 have standardized RMS .743005 and .745647. Their P4 replays, and the generating card, pass all three scenarios: **9/9** sampled cases. Start 0 begins at the generating card. Start 2 supplies a genuinely displaced start. These are two local fits, not a full uncertainty ensemble.

The residual denominator is `sqrt(noise_SD² + discretization_RMS²)` per test/channel. The discretization estimate is measured at the **generating card**, not established uniformly over the parameter box. CAL_SLOW force has discretization RMS **823.028 N** and maximum absolute discrepancy **2466.897 N**; its effective denominator is about **823.27 N**, roughly **41 times** the stated 20 N measurement noise. The raw report attributes the maximum to an elastic-rise/interpolation sample at t=1 s. This is an unresolved sampling/discretization problem, not additional measurement noise, and the one channel-wide denominator downweights the other slow-force samples too. Gauge RMS discrepancy .00576–.01176 mm also exceeds the .002 mm noise. Prediction coupling settings (.5 K / cap 3) additionally differ from discrepancy-generation settings (.25 K / cap 4).

Accordingly, “noise-consistent” must mean only a low residual under this **ad hoc discrepancy-augmented scaling**; it does not certify measurement-noise consistency, statistical coverage or satisfactory early-time interpolation. The RMS also includes 15 exact initial entries among 135 residual components. No residual-coverage likelihood test, numerical-error uncertainty validation, profile or bootstrap is present.

Starts 1 and 3 finish with the constant solver-failure penalty RMS 20, not accepted fits. In the resumed run only 0,2,3 were selected; start 1 belongs to the earlier aborted run. Start 3's changed returned vector and success flag must not be interpreted as a fitted card or evidence that every intermediate evaluation failed—the saved artifacts do not contain that evaluation history. Its design failures are not uncertainty counterexamples. Fitted cards were replayed only on base P4, not across the refinement/independent/neighborhood cross-product. Thus broad uncertainty-preserved numerical feasibility remains **UNVERIFIED**.

## 7. Reproduction and provenance boundary

Historical commands, from P/src with one thread:
```
python3 reproduce_precorrection_proportional.py
python3 test_acceptance.py
python3 run_corrected_pilot.py
FIT_STARTS=0,2,3 python3 finish_corrected.py
```
Recorded environment: Python 3.12.3, NumPy 2.1.3, SciPy 1.14.1; seed 20261108. These commands were not run by this audit. The first reproduction script imports whichever `coupled_current.py` is present; it must historically precede the fix and does not reconstruct old source from this frozen corrected tree.

The full corrected driver **aborted during calibration** with an 80-versus-135 residual-length error (stdout EXIT:1). Its sweep/refinement/independent observable histories were already saved. The delivered source fixes the penalty length to 135. The finish command then exited 0, assembled 79 successful non-card history rows and four logged primary failures, and reran calibration/card replays. The archive now has 91 observable-history files, including 12 card replays (three of these are the rejected fit3 card).

That recovery is legitimate partial-run evidence, **not a clean all-stages success after the fix**. `finish_corrected.py` reads `/tmp/corrected_stdout.txt`, which is not at that path in the archive; the separate stdout attachment contains its historical contents. It reuses saved simulation histories and does not rerun those stages. Same-node results also remain log-only. No clean end-to-end run of the final source is recorded. Clean-run reproducibility is **UNVERIFIED**, not PASS.

## Reviewed gates

PASS is restricted to its named stored-evidence scope. FAIL denotes a tested numerical failure or a contradicted difficulty requirement, not automatic physical infeasibility.

| Gate | Verdict |
|---|---|
| Recorded pre-/post-fix proportional family 12/12 | PASS |
| Acceptance fix / scoped regression evidence / saved diagnostic tolerances | PASS |
| P4/P7 refinement label stability within stated coverage | PASS |
| Four-variable search required or demonstrated harder | FAIL |
| Same-integrator k-only warm feasibility flips | PASS |
| Independent conducting coarse P4/M2 all-scenario label checks | PASS |
| Shared-node warm label agreement, mismatched time settings disclosed | PASS |
| Complete independent converged all-scenario/all-constraint accuracy | UNVERIFIED |
| M2 nominal coarse hold solver robustness | FAIL |
| Primary M2 k=0 completion | FAIL |
| Two local discrepancy-scaled fit P4 replays | PASS |
| Noise-only calibration consistency / full uncertainty coverage | UNVERIFIED |
| Uncertainty × refinement × neighborhood numerical stability | UNVERIFIED |
| True asymptotic spatial/temporal convergence | UNVERIFIED |
| Initial-inclusive accepted-step histories / reported extrema consistency | PASS |
| Continuous interstep extrema certificate | UNVERIFIED |
| Independent coupled energy audit | UNVERIFIED |
| Clean all-stages rerun after final fixes | UNVERIFIED |
| Frontier difficulty | UNVERIFIED |

Machine-readable scopes and failures are in `reviewed_corrected_gate_verdicts.json`; history integrity checks are in `history_evidence_checks.json`.

## Final recommendation

**RECONSIDER**, on the evidence already obtained: the simple one-variable no-hold family succeeds with unchanged requirements, and both the corrected primary and selected independent conducting checks support that result. Conduction is material; necessity of four-variable schedule optimization and increased frontier difficulty are not demonstrated. This conclusion does not require inventing new physics, tightening limits, claiming a clean run that did not happen, or proceeding to a final verifier/ZIP. No such additional work was undertaken.
