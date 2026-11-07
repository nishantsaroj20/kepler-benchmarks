# Provenance

The bus is synthetic. The length, the capacitance values, and the inverter coefficients are not taken from a process design kit or a published interconnect paper.

`environment/data/problem.json` is the full problem in numbers. `authoring/generate_truth.py` integrates a 48-segment ladder at 0.1 ps and writes `tests/truth.json`. The reference solution uses 24 segments and 0.2 ps. There is no random seed.

The uncoupled run in `authoring/evidence/convergence.json` repeats the fine transient with the neighbor capacitance set to zero. That run is the check that coupling, not the ladder spacing, is what sets the crossing times.
