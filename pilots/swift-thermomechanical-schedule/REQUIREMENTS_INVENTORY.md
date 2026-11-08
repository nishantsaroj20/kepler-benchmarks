# Requirements inventory for a later human-written instruction

This is an inventory, not `instruction.md`. Nothing here is prose for the agent-facing brief. Stages 4–8 have not been built. Paths under `/app` and `tests/` are the contract a later bundle would have to implement. They do not exist yet.

## Outcome the instruction would have to require

From the public files only, produce a material card and one schedule `d = (t1, t2, t3, h)` such that every public scenario, replayed by an independent solver the agent does not supply, satisfies all of the following. Any legal schedule passes. Matching a reference schedule is not required.

## Public inputs the instruction would name

- Geometry, elastic constants, known thermal constants, parameter bounds, scenarios, schedule family, sampling operators and the numerical limits in `specification.md` section 9, which the pilot validated.
- Calibration CSVs with the columns and noise model in specification section 8. The pilot file is `results/calibration_public.csv` (SHA-256 `d43146d85f1789f1db43065ef4e9bcc5d30ee685298700f0f59b4e26b6baf12b`, seed 20261108). A later bundle would copy a regenerated file into `environment/data/` and quote its hash.
- Units: mm, s, N, MPa, K. Full 360° force. Gauge opening is `2 u_z` of the outer-surface material point that starts at `Z = 4 mm`. Temperatures are reference-volume averages on `0 ≤ Z ≤ 1.5 mm` (neck, capped) and `6 ≤ Z ≤ 8.5 mm` (shoulder, reported, not capped).

## Schedule definition the instruction would state exactly

- Increments `Δ1 = 0.70 mm`, `Δ2 = 0.75 mm`, `Δ3 = 0.55 mm` on the half-model end face.
- Shape `s(ξ) = ξ - sin(2πξ)/(2π)` on each ramp. `v_max = 2δ/t`, `a_max = 2πδ/t²`.
- Hold at the stage-2 / stage-3 boundary only. `h = 0` means stage 3 starts immediately. No extra smoothing.
- `τ = t1 + t2 + t3 + h`.
- Inclusive bounds: t1 in [0.34, 12] s, t2 in [0.35, 12] s, t3 in [0.30, 12] s, h in [0, 16] s, each ramp positive.
- Illegal if `v_max > 8 mm/s` or `a_max > 40 mm/s²`.

## Constraints the instruction would state as requirements

Evaluated on the replay solver’s accepted steps, including steps inside ramps and inside the hold.

- Neck-band average ≤ 345 K at every accepted step, every scenario.
- `τ ≤ 8 s`.
- Final gauge opening in [2.40, 2.95] mm, every scenario.
- Minimum current outer radius ≥ 2.20 mm at every accepted step, every scenario.
- Axial force ≤ 14000 N at every accepted step, every scenario.

Comfortable numerical margins observed in the pilot, not extra hidden tests: 5 K, 1 s, 0.05 mm on the gauge window, 0.05 mm on radius, 500 N. A later verifier must publish its own refinement rule and its acceptance margins. It must not use a tighter secret test.

## Agent artifacts a later instruction would name

Proposed, not implemented:

| Path | Contents |
|---|---|
| `/app/output/material_card.json` | Eight parameters inside the inclusive bounds, plus `sigma_y = K ε0^n` and the analytical flow checks |
| `/app/output/calibration_predictions.csv` | Predictions of the public calibration rows |
| `/app/output/schedule.json` | Only `t1, t2, t3, h` and a schema version |
| `/app/output/design_predictions.csv` | Force, gauge, minimum radius, neck and shoulder temperatures at public normalized times, per scenario |

Predicted histories would be compared with replay, not accepted as proof that a constraint was met.

## Verifier behavior a later bundle would need

- Reconstruct the schedule from the public definition. Do not execute submitted code.
- Run a solver that does not import the reference constitutive update, assembly or time integrator.
- Use the sealed generating card only after the public data have been shown, as in this pilot, not to invent a pass/fail the public cards disagree on.
- Two fitted cards and the generating card agreed on M2, FAST and SLOW. That agreement is the fairness evidence so far. It is not a finished verifier.
- Binary reward, pytest 8.4.1, pytest-json-ctrf 0.3.5, offline verifier image, agent timeout exactly 28800 s, verifier timeout measured rather than left at 900 s.
- Instruction ending, when it is eventually written by the human author, must be the guide’s exact sentence: `You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.`
- Plain prose, no Markdown headings, in that future `instruction.md`.

## Physics statements the instruction would have to keep

- Multiplicative J2, Hencky elasticity, Swift–Voce plus the multiplicative rate and temperature factors. Eight unknowns. Conductivity, heat capacity, `chi`, grip temperature and the scenarios are known.
- Heating and yield use Kirchhoff `q` per reference volume. Conduction uses `k J C^{-1}` pulled back from the true `F`.
- No thermal expansion and no reversible thermoelastic heat.
- No damage, fracture, contact, anisotropy, kinematic hardening, or inertia. Quasi-static because the actuator caps produce a few pascals of inertial stress.
- Neck temperature is a material-coordinate average, not a pointwise maximum.

## Known solver gaps a later implementation must not ignore

- Hold steps at `u_end = 1.45 mm` on the deeper imperfection, and some `k = 0` replays, failed Newton and were not forced through by retuning. A verifier that dies on a legal schedule is a wrong fail.
- The independent adiabatic temperature on a 3×6 mesh differs by up to 6.1 K from the primary adiabatic temperature. Force agrees to 0.06%. A verifier tolerance tighter than that gap, on that mesh, would punish a correct independent code.
- Radius and force limits did not bind. Do not describe them as the source of difficulty.
- Numba, JAX, surrogates and extrapolation stay allowed for agents. This pilot did not use them.
