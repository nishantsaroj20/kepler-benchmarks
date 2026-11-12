Files are UTF-8. Each CSV file has one header row and no index column. Each JSON file is one object. Dates are calendar dates in the form YYYY-MM-DD. Stirrer speeds are revolutions per second. Particle diameters are micrometers. A pressure written as gauge is bar gauge. The site atmospheric pressure in bar, the slurry volume in cubic meters, and the baffles flag are in `/app/data/vessel.json`. Another campaign directory follows this same layout and these same column definitions.

`/app/data/rate_runs.csv` is the hydrogenation campaign. A row is a measurement only when `status` is `ok`. Any other status is an aborted export and is not a rate. `tic201_C` is thermowell TIC-201 in degrees Celsius. `P_total_barg` is the headspace gauge pressure. `m_cat_g` is the dry catalyst charge in grams. `h2_peak_area` is the chromatograph area for the hydrogen consumed during `window_s` seconds. `sheet_rate_mol_m3_s` is the number exported by the acquisition spreadsheet.

Moles of hydrogen equal `h2_peak_area` times `response_factor_mol_per_count` from `/app/data/gc_standards.csv`. Use the standards row whose `start_date` and `end_date` both include `run_date`. The slurry-volume-specific rate, in mol m$ ^{-3} $ s$ ^{-1} $, is those moles divided by `window_s` and by the slurry volume. The substrate is in large excess and each window stays at low conversion, so this rate is the initial hydrogen rate. Catalyst mass concentration w, in kg m$ ^{-3} $, is the charge in kilograms divided by the slurry volume. Every historical rate was measured after the catalyst had been conditioned for 4 h at that run's stirrer speed.

Absolute headspace pressure in bar is `P_total_barg` plus the site atmospheric pressure. `/app/data/solvent.json` names the solvent and gives its vapor pressure through log10($ P_{sat} $ / 1 bar) = A - B / ($ T_K $ + C), with $ T_K $ the reference-probe temperature in kelvin and $ P_{sat} $ in bar. Hydrogen partial pressure in bar is the absolute headspace pressure minus that vapor pressure.

`/app/data/h2_solubility.csv` lists the dissolved-hydrogen concentration, in mol m$ ^{-3} $, at a hydrogen partial pressure of 1 bar, against reference-probe temperature in degrees Celsius. Between adjacent rows, interpolate the natural logarithm of that concentration linearly in reciprocal thermodynamic temperature. Outside the table, continue that line from the nearest end interval. The interface concentration is the interpolated 1 bar value multiplied by the hydrogen partial pressure in bar.

`/app/data/thermocouple_calibration.csv` pairs `tic201_C` with `reference_C`. Interpolate `reference_C` linearly in `tic201_C`. Outside the table, continue the nearest end interval. Temperatures in `/app/data/predict_cases.csv` are already the reference probe, in degrees Celsius.

`/app/data/kla_measurements.csv` gives the baffled-vessel gas-liquid coefficient `kla_per_s`, in s$ ^{-1} $. `/app/data/kc_measurements.csv` gives the liquid-solid coefficient `kc_m_per_s`, in m s$ ^{-1} $. A row with a blank stirrer speed is not usable. A row whose note contains `baffles removed` is not the baffled campaign vessel. Use the usable gas-liquid row at the charge's stirrer speed, and the usable liquid-solid row at that speed and particle diameter. External catalyst area per slurry volume is $ 6 w / (\rho_p d_p) $, where $ \rho_p $ is the particle density in `/app/data/catalyst.json` and $ d_p $ is the particle diameter in meters. The liquid-solid volumetric coefficient is `kc_m_per_s` times that area.

The catalytic step scales in proportion to CO chemisorption uptake, in micromoles per gram, divided by the fresh uptake in `/app/data/catalyst.json`. `/app/data/catalyst_sites.csv` lists that uptake by `catalyst_state`. Use the row for the charge's state and stirrer speed. A blank stirrer speed on a state applies at every speed.

Particles are spheres of diameter `dp_um`. `metal_location` in `/app/data/catalyst.json` is `uniform` or `eggshell`, and `effective_diffusivity_m2_per_s` is the diffusivity $ D $ of dissolved hydrogen inside the particle, in m$ ^{2} $ s$ ^{-1} $. The hydrogenation is Langmuir-Hinshelwood in dissolved hydrogen. One campaign has one set of four constants. $ k_{75} $ is the saturated surface coefficient at 348.15 K, in mol kg$ ^{-1} $ s$ ^{-1} $. $ K_{75} $ is the adsorption coefficient at 348.15 K, in m$ ^{3} $ mol$ ^{-1} $. $ E $ and $ \Delta H $ are in J mol$ ^{-1} $. With $ R = 8.314462618 $ J mol$ ^{-1} $ K$ ^{-1} $ and $ T_K $ in kelvin,

$ k(T) = k_{75} \exp[-E/R (1/T_K - 1/348.15)] $

$ K(T) = K_{75} \exp[-\Delta H/R (1/T_K - 1/348.15)] $

Let $ U $ be the CO uptake selected for the charge and $ U_{fresh} $ the fresh uptake. At local concentration $ c $, the uniform-particle volumetric rate, in mol m$ ^{-3} $ s$ ^{-1} $, is

$ r_v = \rho_p (U/U_{fresh}) k(T) \dfrac{K(T) c}{1 + K(T) c} $

For `eggshell`, multiply that rate by $ R_p^3 / (R_p^3 - (R_p - \delta)^3) $, where $ \delta $ is `eggshell_depth_um` in meters. The core inside the shell is inert and impermeable. Steady diffusion in the metal-bearing region is

$ \dfrac{1}{r^2} \dfrac{d}{dr} \left( r^2 D \dfrac{dc}{dr} \right) = r_v(c) $

For `uniform` the region is $ 0 < r < R_p $ and no flux crosses the center. For `eggshell` the region is $ R_p - \delta < r < R_p $ and no flux crosses the inner face. The outer surface is at concentration $ c_s $. The slurry-volume-specific rate is

$ r = \dfrac{w}{\rho_p V_p} 4 \pi R_p^2 D \left( \dfrac{dc}{dr} \right)_{R_p} $

with $ V_p $ the particle volume. The two films close the surface concentration: $ c^* - c_b = r / (k_L a) $ and $ c_b - c_s = r / k_s $, where $ k_L a $ is `kla_per_s` and $ k_s $ is the liquid-solid volumetric coefficient. The four constants are fixed by the corrected historical rates.

`/app/data/predict_cases.csv` gives the charges to forecast, including `hold_s` and `headspace_m3`. The reported rate and the three fractions are the initial steady values at the starting hydrogen partial pressure. `resistance_fraction_gas_liquid` is $ (c^* - c_b) / c^* $. `resistance_fraction_liquid_solid` is $ (c_b - c_s) / c^* $. `resistance_fraction_kinetic` is $ c_s / c^* $. `hydrogen_consumed_mol_m3` is the slurry-volume-specific hydrogen consumed during `hold_s` seconds, the time integral of the instantaneous rate. The headspace is closed. Its initial hydrogen inventory, in moles, is the starting hydrogen partial pressure in pascal times `headspace_m3`, divided by $ R T $, with $ R = 8.314462618 $ J mol$ ^{-1} $ K$ ^{-1} $ and $ T $ the reference-probe temperature in kelvin. Solvent vapor stays at the solvent vapor pressure. While the hydrogen partial pressure is positive, its change in bar per second is $ -r V_{slurry} R T / (V_{head} \times 10^{5}) $. The partial pressure never goes below zero, and both the rate and further consumption are zero while it is zero.

`/app/data/effectiveness_check.csv` is a crushed-versus-whole side test. It is not a row of `/app/data/rate_runs.csv`. `/app/data/operator_note.txt` is the operator's weekly note. It is not the specification of the rate or the forecast.
