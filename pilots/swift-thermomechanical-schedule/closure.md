# Closure: Swift thermomechanical schedule pilot

Stop this candidate. The final evidence audit is the source of gate wording. This note does not replace that audit and does not rerun any simulation.

## What stands

Conduction changes the warm result by tens of kelvin, enough to flip the 345 K neck-band limit. The same no-hold proportional schedules, `(t1, t2, t3, h) = (0.35, 0.375, 0.275, 0) τ` at `τ = 4, 5, 6, 7` s, pass all three public scenarios under the unchanged requirements, both before the acceptance-rule fix and after it. Twelve of twelve stored base rows pass on each solver. Warm `τ = 4` is 340.239 K at the neck, with 4.000 s of cycle-time headroom, not 2.60 s. The nearer gauge margin is 0.201 mm below the 2.95 mm ceiling. Four-duration schedule hardening is not demonstrated. M2 and M3 were hand-picked. No schedule optimizer ran.

The 8 s cap stays. Simple families stay allowed. No extra unknowns are added to push a passing schedule out.

## What the audit leaves unresolved

These are not passes. The raw correction overclaimed several of them.

- Frontier difficulty is UNVERIFIED. No fresh probe was run.
- Primary M2 still fails at the hold: nominal, 4 substeps, `k = 45`; and all three scenarios at `k = 0`. Those are unfinished solves, not proof the schedules are physically illegal.
- Full calibration uncertainty is UNVERIFIED. Two local fits (standardized RMS 0.743 and 0.746) replay P4 successfully. That RMS uses a discrepancy-augmented denominator. On the slow record the force denominator is about 823 N against 20 N of measurement noise, driven by the `t = 1` s sample. This is not a noise-only certificate. Starts that returned only the solver-failure penalty are not accepted cards.
- Clean end-to-end reproduction of the final source is UNVERIFIED. `run_corrected_pilot.py` exited 1. `finish_corrected.py` exited 0 by reusing saved histories and a historical `/tmp` log.
- Independently converged accuracy on every scenario and constraint is UNVERIFIED. Coarse conducting checks and one warm shared-node comparison agree on labels, but the time-step settings do not match, and shared-node histories were not saved as npz. The old 0.174 mm gauge split is not thereby closed.
- Independent coupled energy audit is UNVERIFIED. Both heat balances are internal residuals of the solver that produced them.
- Asymptotic mesh and time convergence, and continuous extrema between accepted steps, are UNVERIFIED. Warm P4’s primary temperature spread is 1.607 K and its gauge spread is 0.093 mm. Labels held. That is sensitivity evidence, not an error bound.

## Supersession, without rewriting evidence

The imported audit supersedes overbroad sentences in `correction.md` and in `results/corrected/gate_verdicts.json`. Those files are left byte-for-byte as committed in `e60bab59690a5bc02624f3fa1fb3fdb02a1e181b`. `individualreport.md` remains the first report, already marked superseded by `correction.md`, and now also superseded by this closure and by `provenance/imported_final_review_2026-11-08/`. Nothing in `results/`, `logs/`, or `src/` from that commit was edited for this closure.

The review copies are static imports. They are not new run data.

## Freeze point

Evidence files are frozen by the SHA-256 list in `freeze_manifest.sha256`, checked from this directory. The manifest does not hash itself. Exclusions are listed in that file’s header comments. The git commit that adds this closure, the manifest, and the imported review is the retrieval commit for the whole pilot tree, including failed runs and the earlier reports. The repository is not an access-control system. The original task ZIP was not opened for editing. Its SHA-256, rechecked as a static file and not repackaged, is `0c154a40d1f6d594dab71b35147d6abf3fb66a69589d6593c2d491d4b95fc43d`.

Parent evidence commit before this documentation-only freeze: `e60bab59690a5bc02624f3fa1fb3fdb02a1e181b`.
