# Concept: one outer controller for two grid-following LCL converters

Working name: `afterquery/converter-robust-control`. Field: Electrical & Computer Engineering. This note is authoring rationale, not an agent instruction.

## Equipment and the decision

The synthetic installation is two three-phase, two-level voltage-source converters on a 400 V, 50 Hz feeder. Each unit is rated 100 kVA and has its own LCL filter and a short cable to a common point of connection. Behind that point is a utility Thevenin source whose short-circuit strength and X/R are known for each published case. The DC buses are stiff at 750 V. No laboratory, site, or measured waveform is claimed.

The engineering decision is whether one set of outer-loop coefficients, commissioned once, can keep both inverters at their scheduled active and reactive power when the feeder strength, the power schedule, or the cable asymmetry changes among four published operating points. Inner current regulators, phase-locked loops, sensor filters, and the PWM lag stay at factory values. They are part of the plant.

## Why grid-following

The feeder already has a voltage source at a fixed nominal frequency. The converters are scheduled to exchange power with that feeder. They are not the sole sources of an island. Grid-following operation is the matching architecture: each inverter locks a local dq frame to its own filtered capacitor voltage and regulates current in that frame. Grid-forming droop would be the architecture for an islanded microgrid, which is not the installation described here.

Huang, Xin and Dörfler (arXiv:1906.11331) use weighting functions to make one converter look grid-following or grid-forming in small signal, and they give a decentralized small-gain certificate for a network. That paper supports the relevance of structured H-infinity design across grid conditions. It is not the model used here. Their numerical example is a 50 kVA, 380 V converter with 0.05 pu LCL elements and a static 3-by-7 gain. None of those values or that controller are used.

## Why one controller, and why it is local

The four cases are operating conditions of the same two inverters, not four different products. A commissioned firmware image should not be rewritten when the scheduled power changes or when one lateral is longer. There is no control network between the cabinets, so each block may use only its own power errors. The same numerical blocks are applied at every case. Case-specific current-reference biases are the scheduled equilibrium, not extra controller freedom.

## Measurements, actuators, and order

Each inverter measures its own filtered capacitor voltage and filtered converter-side and grid-side currents in its PLL frame. The outer controller sees only the two normalized power errors formed from the filtered voltage and the filtered grid-side current. It commands the d- and q-axis current-reference perturbations, in per-unit of the rated peak current, about the published bias.

Two controller states per inverter are the classical choice: one integrator for active power and one for reactive power. A general 2-by-2 dynamic block includes that PI pair and also cross terms. The order is kept because it is the industrial outer loop, not because a larger order would be harder to search. `Ak` is not required to be Hurwitz. A zero eigenvalue is an integrator and is legitimate if the closed loop is internally stable.

## What the requirements mean

Internal stability is a decay-rate requirement on the physical plant plus the controller, not on the analysis weights. The performance requirement is one stacked H-infinity norm of a weighted map from normalized grid-voltage, power-reference, and measurement-noise channels to weighted power errors and weighted current-reference effort. The pass level is 1. A witness used for numerical comfort is tighter, but it is not a second physical specification. See `REQUIREMENTS_RATIONALE.md`.

## Limits

The four linearizations do not certify behavior between the cases, under switching from one case to another, under PWM, or under current or voltage saturation. Anti-windup is not modeled. No continuous parameter set is covered. Equilibria sit inside sine-triangle modulation headroom; that is an audit of the chosen operating points, not a nonlinear robustness claim. No real-device certification is implied.

## Simplest plausible solution

A tied diagonal PI, identical on both inverters, is the cheapest competent experiment. The executed pilot found a wide set of such gains that meet the frozen specification. The pilot verdict is therefore STOP. Details and the raw runs are in `PILOT_REPORT.md`.
