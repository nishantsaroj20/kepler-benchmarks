# Provenance

Practice and held-out packs are synthetic. Generator: `generate_packs.py`.

- Practice seeds: 101–104 (fixed packs copied into `environment/data/practice/`).
- Held-out packs are not baked into the agent image. The verifier builds eight packs at grade time from `tests/hidden_family.py`, with the top-level seed derived from SHA-256 of the submitted `/app/solver/analyze.py` bytes (salt `contact-shift-v1`). Patch centers, radii, clearances, $E$, $\nu$, $\mu$, node counts, and load-path family are redrawn so cases are shifted relative to practice rather than reseeds.
- Continuum model: rigid punch on isotropic plane-strain half-space, uncoupled identical-material operators, piecewise-constant integrated log kernel with pack `rref`, Coulomb friction, force-controlled history.
