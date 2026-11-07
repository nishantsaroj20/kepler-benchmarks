# Short-pulley shrink fit

## Difficulty

The hub is short relative to its wall, so the contact pressure is higher near the ends than at mid-length. A uniform Lamé pressure, even with the right interference, does not meet the station requirements. The operating case then superimposes aluminum growth, a smaller steel growth, and centrifugal body force. The delivered speed is the root of that torque capacity, not the assembly pressure read off a formula. The geometry and properties are synthetic.

## Reference solution

Stuck axisymmetric linear quadrilaterals close half the stated diametral interference. Each part carries its own modulus, expansion coefficient, and density. Torque is the friction integral of pressure over the interface, and speed is bisected inside the compressive bracket from 1000 to 3000 rad/s. The reference mesh is coarser than the mesh that produced the sealed truth.

## Verification

The truth file is the fine-mesh result from `authoring/generate_truth.py`. Tests compare the five assembly pressures, the assembly torque, and the limit speed with the bands stated in `instruction.md`. Those bands sit above the recorded medium-to-fine spread in `authoring/evidence/convergence.json` and still reject a uniform long-cylinder pressure. The verifier only sees `/app/output/fit_result.json` plus the files baked into `tests/`.
