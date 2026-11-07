# Provenance

The pulley and spindle are synthetic. No measured machine is represented, and nothing is drawn from a published example with this geometry.

`environment/data/problem.json` is the full problem statement in numbers. `authoring/generate_truth.py` solves it on a 16-by-16 by 32 contact mesh and writes `tests/truth.json`. The reference solution solves the same equations on an 8-by-8 by 16 mesh. There is no random seed; both solves are deterministic.

The long equal-material sleeve check in `authoring/evidence/convergence.json` compares the mid-face pressure with an independent plane-strain Lamé formula. Agreement there is the check that the axisymmetric operator, the hoop strain, and the contact jump are the right ones. The short aluminum-on-steel hub is then kept because that same uniform formula misses the station pressures.
