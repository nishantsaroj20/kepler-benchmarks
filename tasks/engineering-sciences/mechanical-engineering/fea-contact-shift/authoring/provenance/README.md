# Provenance

Practice and held-out packs are synthetic. Generator: `generate_packs.py`.

- Practice seeds: 201–205 (fixed packs in `environment/data/practice/` and mirrored to `tests/practice_data/` for the verifier image).
- Held-out packs are built at grade time by `tests/hidden_family.py`, with the top-level seed from SHA-256 of `/app/solver/analyze.py` (salt `contact-shift-v2-approach`). Cases are shifted vs practice (geometry, $E$, $\nu$, $\mu$, approach-path family), including shear-reversal and no-contact.
- Continuum model: rigid punch on isotropic plane-strain half-space, uncoupled identical-material operators, piecewise-constant log kernel with pack `rref`, Coulomb friction, **approach-controlled** $(\delta_n,\delta_t)$ history. Forces are solution outputs.
- Ground truth: `tests/truth_solver.py` uses Gauss quadrature for $K_{ij}$ and a separately authored active-set loop (does not import the reference analyzer).
