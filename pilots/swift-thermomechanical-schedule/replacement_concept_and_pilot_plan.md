# Replacement concept and cheapest decisive pilot

One candidate only. Nothing here is implemented, and nothing here is run. No numerical acceptance tolerance is frozen. This is not an instruction.

## The engineering question

A mechanical analyst is given three synthetic tensile records and the written specimen definition, and is asked for a certificate rather than a point forecast. The records are the ones already in the original task: `tensile_qs.csv`, `tensile_v30.csv`, and `tensile_v3000.csv`. Each file has exactly two channels, `extensometer_opening_mm` and `axial_load_N`. No temperature, radius, or notch load is in those files. Geometry, the Swift–Voce rate/temperature family, the known elastic and thermal constants, and the parameter bounds are in `specs.json`. The observations are synthetic. No laboratory campaign and no personal design project are available, and none is invented here.

The question is: using only those public channels, what interval of notch peak load is implied by every parameter vector that agrees with the three curves, and is that interval small enough to be a certificate or so wide that the records do not support one?

A single fitted card followed by one forward notch simulation does not answer it. That is the job the earlier task already asked, and capable solvers already did it. Reporting one feasible card, or searching a loading schedule with a known law, is the same kind of job again. The decision here is whether the public curves determine the structural prediction tightly enough that every adequate explanation lies in a useful band.

## What would count as success, later

Not specified as a number now. The shape of the criterion is fixed so it cannot be tightened after the pilot to manufacture difficulty.

- The output is an interval of notch peak load for each of the two notch geometries already defined in `specs.json`, plus the finite set of cards used to build it.
- Every card in the set must agree with all three public load–opening curves under a discrepancy rule written down before the notch results are inspected.
- The interval must cover the notch peak loads of those cards.
- An independent calculation, using a card that meets the same public discrepancy rule and is not one of the submitted cards, must land inside the interval.
- If the only way to pass is to emit an interval so wide that it contains every peak load allowed by the published parameter bounds, the records are not sufficient and the attempt fails for fairness, not for lack of a solver.

The v1 relative tolerances on analytical flow stress and on notch peak load already exist in that task. The pilot may use them as a pre-existing yardstick for whether a band is trivial. They are not retuned here.

## Cheapest pilot, before any implementation

**Strongest shortcut first.** Do not start with a finite-element mesh. Convert the pre-peak portion of each public curve to approximate true stress and true strain with the usual incompressible round-bar formulas, using the gauge length and radius already stated in `specs.json`. Fit the eight-parameter Swift–Voce rate/temperature expression by bounded least squares, in closed form, with no spatial discretization. Keep every trial card whose pre-peak residual is within a small multiple of the stated noise, 15 N and 0.0015 mm. That multiple is chosen before looking at notch loads; three and five are enough to see whether the conclusion depends on the multiple. Evaluate only the analytical flow stress at the six plastic strains already listed in `specs.json`, and the two dynamic flow states already listed there. Those quantities do not require a notch mesh.

**What this tests.** One generating card is already stored with the original task, together with its analytical flow stresses in `tests/ground_truth.json`. The pilot uses that card only as an offline check that a certificate would have covered the known synthetic truth. It is not a fitting target and it is not a new hidden case. The public curves are the three existing CSVs. The hypotheses are: the pre-peak formulas already force every noise-compatible card into the existing flow-stress tolerances; or they admit a spread larger than those tolerances; or they admit no card at all except by violating the noise model.

**Proposed budget, not a measurement.** A few hundred to a few thousand algebraic evaluations, one NumPy script, well under an hour on a laptop. No Abaqus run and no axisymmetric mesh in this pilot.

**Go.** Continue toward a later implementation only if all three of the following hold. At least two cards with materially different `K` and `α` both fit all three pre-peak curves inside the predeclared noise multiple. Their analytical flow stresses at one or more of the published report strains differ by more than the existing v1 relative tolerance. Adding the post-peak points of the same three curves, still with a one-dimensional Considère estimate and still without a mesh, does not remove that difference. The essential content is then “which features of these curves certify the prediction,” and a pure pre-peak fit is not enough. A later finite-element check would be justified only as a test of that claim, not as the starting point.

**No-go, abandon.** Stop if the algebraic ensemble already lies inside the existing flow-stress tolerances. The certificate is then a low-dimensional fit, the structural mesh is unnecessary, and this candidate is the same job as “fit the known law.” Also stop if the ensemble is empty unless the noise multiple is made absurd, or if the band is essentially the whole parameter box. In those cases the three channels are not a fair basis for a certificate, and the missing measurements should be named rather than simulated into existence. Temperatures, local strains, and notch loads are not in the public CSVs. They are not available.

**Why a failed shortcut would mean something.** The shortcut is the standard pre-peak reduction a tensile analyst would try before any mesh. If it fails the go criteria above, the failure is that distinct hardening mixtures agree with the same load–opening points and still disagree on a published flow stress. That is a property of the observations and the law, not an under-resolved mesh and not a missing column. If the shortcut fails because the script is too crude to fit even the generating card’s pre-peak curve within the noise, the pilot is inconclusive and must be repeated with a better one-dimensional residual before any claim. That repeat is still not a finite-element project.

**Not a post-hoc limit change.** The comparison yardstick is the tolerance set already published for those analytical flow stresses. The pilot does not choose a tighter band after seeing the spread. If the spread sits just outside or just inside that yardstick, the result is reported as borderline and the candidate is not “saved” by editing the yardstick.

A point prediction from one fitted card remains an acceptable scientific answer only when the pilot has shown that every other noise-compatible card gives the same prediction to within the existing tolerance. If that happens, the robust-certificate framing collapses and should be dropped. Legitimate reduced-order and analytical methods are in scope. Nothing here depends on forbidding a class of solver, and nothing here is a prediction that a particular model will fail.
