# Hidden campaign

The hidden campaign is not a new noise seed of the visible toluene week. `hidden_campaign()` in `authoring/provenance/generate_campaign.py` writes `tests/hidden_campaign/`. `write_truth()` in the same file evaluates the latent parameters on an 80-node, 200-step finite-volume grid and stores `tests/truth/expected.json`. That grid is finer than the reference program. The truth is not the reference fit.

`assert_shift()` checks that no visible charge is repeated. It also checks that the visible week already contains every condition that decides a hidden charge: uniform particles, a curved thermocouple, three chromatograph windows, an aborted row, a regenerated state, and a hold whose initial uptake is several times the hydrogen inventory. The solvents differ. The hidden sieve cuts are not the visible sieve cuts.

What changes, without introducing a new rule:

- Different fill. `solvent.json` is anisole, with its own Antoine constants. `h2_solubility.csv` is the anisole table. The visible solvent is toluene, and the vapor-pressure and solubility rules are the same.
- Stronger internal diffusion, same particle balance. Both campaigns set `metal_location` to `uniform`. The hidden diffusivity is lower and the forecast diameters are 230, 310, 390, and 470 µm. The visible diameters are 120, 280, and 420 µm. The visible crushed-versus-whole check is slower for the whole particle, and the visible 120 µm and 420 µm rates differ.
- Values outside the visible window. Bath temperatures run from 34 °C to 153 °C, stirrer speeds are 3, 8, 14, and 20 rps, and catalyst masses run from 1.6 g to 150 g. The visible speeds are 6, 9, 12, and 16 rps. Each hidden speed and diameter has its own row in the hidden coefficient tables.
- Same categories. `regenerated` is in the visible sites table and on visible charge `G1`. The hidden table has its own regenerated uptake.
- Same faults, different records. Both campaigns have aborted rows, three chromatograph windows, a blank tachometer row, a baffles-removed row, and a thermocouple table that is not a straight line. The hidden curve and the hidden calendar are not the visible ones.
- Same edge. Visible charge `S1` is a 5000 s hold in a 0.00012 m³ headspace. Hidden charge `E03` also empties a small headspace. `E01` is near boiling. The vapor-pressure subtraction is the visible rule.

The factorial map and the four edge charges are written explicitly in `hidden_campaign()`. The noise seeds in `authoring/provenance/seeds.json` apply only to the uptake noise. The rate equation, the temperature laws, the particle balance, and the closed-headspace balance are in `environment/data_note.md`.
