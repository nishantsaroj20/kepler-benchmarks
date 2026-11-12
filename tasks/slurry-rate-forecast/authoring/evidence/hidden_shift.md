# Hidden campaign

The hidden campaign is not a new noise seed of the visible toluene week. `hidden_campaign()` in `authoring/provenance/generate_campaign.py` writes `tests/hidden_campaign/`. `write_truth()` in the same file evaluates the latent parameters on an 80-node, 200-step finite-volume grid and stores `tests/truth/expected.json`. That grid is finer than the reference program. The truth is not the reference fit.

`assert_shift()` checks that no visible charge is repeated, that the solvent and metal location differ, that a regenerated state is present, and that aborted rows are present.

The structural shifts are:

- New regime. `catalyst.json` sets `metal_location` to `uniform` and the forecast diameters are 230, 310, 390, and 470 µm. Internal diffusion changes the rate. The visible campaign is eggshell Pd/C at 40, 80, and 150 µm.
- New material. `solvent.json` is anisole, with its own Antoine constants. `h2_solubility.csv` is the anisole table, not the toluene table.
- Values outside the visible window. Bath temperatures run from 34 °C to 153 °C, stirrer speeds are 3, 8, 14, and 20 rps, and catalyst masses run from 1.6 g to 150 g. The visible speeds are 6, 9, 12, and 16 rps.
- Unseen category. `catalyst_sites.csv` includes `regenerated`. The visible table has only `conditioned` and `fresh`.
- Injected faults. `rate_runs.csv` contains aborted rows with unusable areas. `gc_standards.csv` has three windows. `thermocouple_calibration.csv` is curved, not the visible straight line.
- Edges. `E01` is near boiling, with a hydrogen partial pressure of about 0.20 bar. `E03` uses a 0.00012 m³ headspace that the hold empties. Visible holds are 600 s in an 0.008 m³ headspace.

The factorial map and the four edge charges are written explicitly in `hidden_campaign()`. The noise seeds in `authoring/provenance/seeds.json` apply only to the uptake noise.
