K-17 is a monoprotic carboxylic acid recovered from aqueous broth into 4-methyl-2-pentanone. The campaign is synthetic. These notes are the plant record for this log, not a literature correlation to look up, and they do not write out the stage balance.

calibration.json is the meter record. thermocouple_celsius converts the controller reading tic_c to broth temperature in celsius as a0 + a1*tic_c + a2*tic_c^2. The extraction depends on the broth temperature. tic_c is only the controller reading.

ph_mv is the glass-electrode millivolt reading. Use the ph_windows entry whose start and end contain the row date: pH = intercept - (ph_mv - millivolt_offset) / slope_mv_per_ph. ph_sheet is the column the spreadsheet stored. Through 2024-05-18 that column was built with the slope sign reversed; sheet_slope_sign records the sign the sheet actually used, and the electrode was reset on 2024-05-19. ph_sheet is not the broth pH.

organic_assay_blank_mol_m3 sits inside c_org_assay_mol_m3. recovery_sheet was formed from that assay, so it still carries the blank. The aqueous residual does not. The removal measured on a historical contact is 1 - c_aq_mol_m3 / c_feed_mol_m3.

temperature_reference_k, stirrer_reference_rps, and stirrer_exponent are the plant agitation reference. A contact held for t seconds at impeller speed N revolutions per second approaches its equilibrium transfer by the fraction 1 - exp(-kLa * t), with kLa = kLa_ref * (N / stirrer_reference_rps)^stirrer_exponent. kLa_ref was not written down.

runs.csv is the pilot contact log. A row with status ok was an agitated contact. A row with status aborted was an export that was cut off. A row with status unmixed was a batch whose impeller never started. Those last two are not agitated contacts.

shakes.csv is the equilibrium study. Those bottles were held long enough that the transfer was complete. At the same broth pH and temperature, a richer feed loads the ketone by more than a partition coefficient that does not change with concentration. The process chemist’s conclusion from that study is that the extracted neutral acid dimerizes in the ketone, and that a monomer-only distribution does not describe the shakes. Across this temperature window the dimerization equilibrium was treated as constant. The acid dissociation and the monomer partition were not.

stage_count is how many contacts were made. When it is 2, the solvent is split into two equal portions and the raffinate from the first contact is contacted with the fresh second portion. Both contacts see the same broth, the same impeller speed, and the same contact time. S in the tables is the combined organic volume divided by the aqueous volume. The organic concentration that goes on a charge sheet is the solute that left the aqueous phase, divided by that combined solvent volume, and it does not include the assay blank.

charges.csv lists the contacts that have not been assayed. They use the same meters and the same split. Several of them lie outside the pilot ranges of controller temperature, solvent ratio, stirrer speed, and contact time.
