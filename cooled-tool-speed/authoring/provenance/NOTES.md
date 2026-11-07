# Provenance

The tool section is synthetic. It is not a commercial insert, and the dimensions are not taken from a DEFORM tutorial or a published turning example.

`environment/data/problem.json` is the full problem statement in numbers. `authoring/generate_truth.py` solves it on the fine fitted mesh and writes `tests/truth.json`. The reference solution uses the medium fitted mesh. There is no random seed. Both solves are deterministic.

The annulus check in `authoring/evidence/convergence.json` compares the same conduction and convection assembly with the exact radial solution of a constant-conductivity annular wall. The carbide tip then keeps the hole and the softening conductivity, because a solid tool and a constant-conductivity tool both miss the rake stations.
