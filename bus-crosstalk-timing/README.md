# Coupled on-chip bus timing

## Difficulty

Five inverter-driven wires share nearest-neighbor capacitance. The far-end half-supply times move when that coupling is removed, by much more than the acceptance band. The driver law is a level-1 MOSFET, so a linear resistor in place of the inverter is a different circuit. The length, parasitics, and switching pattern are synthetic.

## Reference solution

A medium RC ladder is integrated with backward Euler. Each step solves the inverter currents with Newton. The reported times are the first half-supply crossings, and the limit clock follows from the latest crossing plus the setup time. The reference ladder is coarser than the one that produced the sealed truth.

## Verification

The truth file is the fine-ladder result from `authoring/generate_truth.py`. Tests compare the five crossing times and the limit clock with the bands stated in `instruction.md`. Those bands sit above the recorded medium-to-fine spread in `authoring/evidence/convergence.json` and still reject the uncoupled transient. The verifier only sees `/app/output/bus_timing.json` plus the files baked into `tests/`.
