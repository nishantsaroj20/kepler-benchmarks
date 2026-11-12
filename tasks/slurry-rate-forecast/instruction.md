Forecast the hydrogen rate, the hydrogen consumed over the requested hold, and the three interface-concentration shares for every charge in a campaign. Put the forecast program at `/app/forecast.py`. Running

`python3 /app/forecast.py --campaign /app/data --cases /app/data/predict_cases.csv --output /app/output/predictions.json`

reads the campaign in `/app/data` and writes `/app/output/predictions.json`. The same command is used with other campaign directories. Those directories use the filenames, columns, units, and conventions in `/app/data_note.md`.

`/app/output/predictions.json` is one JSON object with the key `cases`. `cases` is an array with one object per charge in the cases file, and no other charges. Each object has exactly the keys `case_id`, `rate_mol_m3_s`, `hydrogen_consumed_mol_m3`, `resistance_fraction_gas_liquid`, `resistance_fraction_liquid_solid`, and `resistance_fraction_kinetic`. `case_id` matches the cases file. The file is at most 2000000 bytes, it is not a symbolic link, and it contains no duplicate keys. Every numeric value is finite.

`rate_mol_m3_s` is the initial slurry-volume-specific hydrogen rate in mol m$ ^{-3} $ s$ ^{-1} $. Its absolute difference from the true rate is at most the larger of 0.090 times the true rate and 0.00015 mol m$ ^{-3} $ s$ ^{-1} $. `hydrogen_consumed_mol_m3` is the slurry-volume-specific hydrogen consumed over the hold, in mol m$ ^{-3} $. Its absolute difference from the true consumption is at most the larger of 0.060 times the true consumption and 0.20 mol m$ ^{-3} $. Each fraction's absolute difference from the true share is at most 0.030. The three fractions sum to 1 within 0.0001. `resistance_fraction_gas_liquid` is $ (c^* - c_b) / c^* $, `resistance_fraction_liquid_solid` is $ (c_b - c_s) / c^* $, and `resistance_fraction_kinetic` is $ c_s / c^* $. $ c^* $ is the interface concentration, $ c_b $ the bulk liquid, and $ c_s $ the outer particle surface. `/app/data_note.md` states the rate equation, the particle balance, and the closed-headspace balance.

Each run of `python3 /app/forecast.py` must finish within 20 seconds on 1 CPU.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
