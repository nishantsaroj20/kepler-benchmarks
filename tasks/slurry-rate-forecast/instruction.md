The campaign records in `/app/data` come from a baffled stirred autoclave hydrogenating 4-nitrotoluene in toluene over an eggshell palladium-on-carbon powder. The charges that still have to be run are listed in `/app/data/predict_cases.csv`. For each charge, determine the steady initial rate of hydrogen consumption and how the dissolved-hydrogen concentration falls from the gas-liquid interface to the outer surface of the catalyst.

Write `/app/output/predictions.json` as a regular file of at most 65536 bytes. It must not be a symbolic link. The file is one JSON object whose only key is `cases`. The value of `cases` is an array of objects, one for each `case_id` in `/app/data/predict_cases.csv` and no others. Each `case_id` appears once. No key is repeated inside any object. Each object has exactly the keys `case_id`, `rate_mol_m3_s`, `resistance_fraction_gas_liquid`, `resistance_fraction_liquid_solid`, and `resistance_fraction_kinetic`. Every numeric value is finite.

`rate_mol_m3_s` is moles of hydrogen consumed per second per cubic meter of slurry. With $ C^{*} $ the dissolved-hydrogen concentration at the gas-liquid interface, $ C_{b} $ the bulk-liquid concentration, and $ C_{s} $ the concentration at the outer surface of the catalyst, the gas-liquid fraction is $ (C^{*} - C_{b}) / C^{*} $, the liquid-solid fraction is $ (C_{b} - C_{s}) / C^{*} $, and the kinetic fraction is $ C_{s} / C^{*} $. The three fractions must sum to 1 within an inclusive absolute tolerance of $ 0.0001 $.

Each rate must lie within the larger of an inclusive relative tolerance of $ 0.040 $ and an inclusive absolute tolerance of $ 0.0006\,\mathrm{mol\,m^{-3}\,s^{-1}} $. Each fraction must lie within an inclusive absolute tolerance of $ 0.030 $.

Units, file layouts, and the measurement conventions that fix the temperature, the hydrogen partial pressure, the interface concentration, the transport coefficients, and the catalyst uptake are in `/app/data_note.md`.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
