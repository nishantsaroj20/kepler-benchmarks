Two structurally different identifications were scored against the latent forward model on all 52 charges.

The reference route is the joint nonlinear least squares in solution/forecast.py. Its maximum relative recovery error is 0.0030347837687722994 and its maximum relative organic error is 0.003035231390944039. It misses no charge.

The other route fits the equilibrium constants on shakes.csv only, with the rate law turned off, and then fits kLa_ref on the agitated pilot rows. Its maximum relative recovery error is 0.005556995228951198 and its maximum relative organic error is 0.005557028972866781. It misses no charge. The largest absolute gap between the two forecasts is 0.00118840947523835 in recovery and 0.19424654080047787 mol m^{-3} in the combined organic concentration.

A monomer-only distribution, fit to the same reconciled aqueous residuals, misses 37 of 52 charges. Fitting ph_sheet, tic_c, and recovery_sheet, and keeping aborted and unmixed rows, misses 46. Copying the mean recovery of the agitated pilot rows misses 49.
