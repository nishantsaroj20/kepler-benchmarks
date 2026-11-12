K-17 is a monoprotic carboxylic acid recovered from aqueous broth into 4-methyl-2-pentanone. The campaign is synthetic. The numbers below are the plant conventions for this log, not a literature correlation to look up.

Only the neutral acid distributes. With broth temperature T in kelvin and T_ref = 298.15 K,

f = 1 / (1 + 10^(pH - pKa))

K = K_ref * f * exp(B * (1/T_ref - 1/T))

pKa, K_ref, B, and kLa_ref are constant across the log and are not written down. Identify them from the historical contacts.

A contact of t seconds at impeller speed N revolutions per second approaches equilibrium by

phi = 1 - exp(-kLa_ref * (N / 3.2)^0.85 * t)

The reference speed 3.2 rps and the exponent 0.85 are fixed.

When stage_count is 1, the solvent ratio S = V_org / V_aq is applied once. The fraction of feed solute that leaves the aqueous phase is

recovery = phi * (K*S) / (1 + K*S)

When stage_count is 2, the solvent is split into two equal portions. The raffinate from the first contact is contacted again with the fresh second portion. Each contact uses solvent ratio S/2, the same K, and the same phi. recovery is the fraction of feed solute removed after both contacts. In either case the organic concentration to report is the concentration in the combined solvent,

organic_mol_m3 = recovery * c_feed_mol_m3 / S

with no assay blank added.

The broth temperature in celsius is a0 + a1*tic_c + a2*tic_c^2 using thermocouple_celsius in calibration.json. tic_c is the controller reading. Do not fit K against tic_c.

ph_mv is the glass-electrode reading. Convert it with the ph_windows entry whose start and end contain the row date:

pH = intercept - (ph_mv - millivolt_offset) / slope_mv_per_ph

ph_sheet is the column the spreadsheet produced. Before 2024-05-19 that column was built with the slope sign reversed. Do not use ph_sheet.

The measured recovery of a historical contact is 1 - c_aq_mol_m3 / c_feed_mol_m3. c_org_assay_mol_m3 still includes organic_assay_blank_mol_m3, and recovery_sheet was computed from that assay. Neither is the fitting target, and the blank is not part of organic_mol_m3.

Drop every historical row whose status is aborted or unmixed. Those rows are not single agitated contacts. The charges to forecast all have status ok, and they still need the temperature and pH corrections above.
