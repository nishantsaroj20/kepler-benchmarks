# Cooled turning-tool tip

## Difficulty

The carbide section is a rectangle with a circular coolant hole. Heat enters on a short rake patch, the hole convects to coolant, and conductivity falls as the tip heats. A solid tool misses the rake stations. A single room-temperature conductivity misses them again and sends the limit speed the wrong way, because the hotter solve is a poorer conductor. The geometry and properties are synthetic.

## Reference solution

A fitted triangular mesh resolves the hole. Conductivity is updated from the element temperature, the shank is held at the stated temperature and clamped, and plane-strain thermal stress gives the corner growth at the reference speed. The limit speed is the bisection of the corner temperature. The reference mesh is coarser than the mesh that produced the sealed truth.

## Verification

The truth file is the fine-mesh result from `authoring/generate_truth.py`. Tests compare the five rake temperatures, the edge growth, and the limit speed with the bands stated in `instruction.md`. Those bands sit above the recorded medium-to-fine spread in `authoring/evidence/convergence.json` and still reject a constant-conductivity rake. The verifier only sees `/app/output/tool_result.json` plus the files baked into `tests/`.
