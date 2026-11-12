Forecast the solute recovery and the organic-phase concentration for every charge in `/app/data/charges.csv`. Write the forecast to `/app/output/predictions.json`.

`/app/output/predictions.json` is one JSON object with the key `charges`. `charges` is an array with one object per row of `/app/data/charges.csv`, in that same order, and no other rows. Each object has exactly the keys `case_id`, `recovery`, and `organic_mol_m3`. `case_id` matches the table. The file is at most 2000000 bytes, it is not a symbolic link, and it contains no duplicate keys. Every numeric value is finite.

`recovery` is the fraction of feed solute that leaves the aqueous phase. Its absolute difference from the true recovery is at most the larger of 0.02 times the true recovery and 0.012. `organic_mol_m3` is the concentration in the combined solvent, in mol m^{-3}, and it does not include the organic assay blank. Its absolute difference from the true concentration is at most the larger of 0.02 times the true concentration and 0.60 mol m^{-3}.

The pilot log, the equilibrium shakes, and the meter record are in `/app/data`. `/app/data_note.md` is the campaign note for that record. Reconcile the log with that note and forecast every charge inside both limits.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
