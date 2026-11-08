# Amendment to the scientific specification

This amends section 11 of `specification.md`. The equations, units, geometry, scenarios, schedule shape, and the numerical limits are unchanged. What changes is the status of those limits as a design problem.

## Limits stay where they are

Neck-band average ≤ 345 K, cycle time ≤ 8 s, final gauge in [2.40, 2.95] mm, minimum radius ≥ 2.20 mm, force ≤ 14000 N. These were not tightened after the proportional family was found to pass. The 8 s cap was not lowered.

## What the corrected runs assign to each limit

- Temperature and conduction are active. On the warm start, `(1.4, 1.5, 1.1, 0)` s is 340.2 K at `k = 45` and 376.1 K at `k = 0` in the same coupled code. The zero-conductivity run fails 345 K. The end face remains a Dirichlet boundary at 293 K for both values of `k`.
- Cycle time does not separate 4, 5, 6, and 7 s. All four are under 8 s. It would reject only a much slower schedule, such as the 14.5 s case in the first sweep.
- Gauge and radius move with the mesh and the time step by several hundredths of a millimetre on the primary solver. On the warm 4 s family they stay inside the limits on every completed primary grid from 3×6 through 5×12 and from 4 to 16 substeps, including the combined 5×12 / 16 case. They are not the constraint that forces a search over four durations.
- Force does not bind.

## Schedule family that already works

`h = 0` and `(t1, t2, t3) = (0.35, 0.375, 0.275) * tau` for `tau` in {4, 5, 6, 7} s. This is one duration, not an optimized `(t1, t2, t3, h)`. Two joint fits that match the regenerated noisy records (standardized RMS about 0.74, with the mesh discrepancy in the denominator) give the same warm 4 s neck temperature as the generating card, to 0.01 K.

## Solver rule that is now part of the numerical specification

A time step is accepted only if the mechanical free residual meets its tolerance and, when conduction is on, the staggered temperature change meets the stated tolerance (0.25 K in the decision runs; 0.05 K and 0.50 K were repeated on the warm 4 s case). A failed Newton or a stagger that hits its iteration cap without meeting the tolerance does not update displacement, plastic strain, the elastic Finger tensor, the deformation gradient, or nodal temperature. The t = 0 state is part of every extremum. Histories of accepted steps are the npz files under `results/corrected/histories/`.

## What is still not specified as validated

A continuous-time extremum between accepted steps is not certified beyond the refinement table. A profile likelihood or bootstrap of the eight parameters was not computed. Two of four fit starts never returned a real residual. The hold at 1.45 mm still fails on some coarse and zero-conductivity runs. None of that is fixed by changing the public limits.
