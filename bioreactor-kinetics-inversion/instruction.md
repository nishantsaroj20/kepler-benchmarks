Two historical fed-batch telemetry files are provided at `/workspace/data/sensor_run_01.csv` and `/workspace/data/sensor_run_02.csv`. Each CSV contains irregularly sampled columns `time_h`, `X`, `S`, `P`, `C_L`, `V`, and `F` (biomass, substrate, inhibitory product, dissolved oxygen, broth volume, and instantaneous feed). Known operating constants, yields, aeration parameters, pump limits, the product-inhibition exponent `gamma`, horizon settings, and the Run 03 initial state are provided in `/workspace/data/reactor_specs.json`.

The culture obeys Haldane–Andrews growth with product inhibition and one metabolic phase shift. With state $ x = [X, S, P, C_L]^{\mathsf T} $ and volume $ V $,

$ \mathrm{d}V/\mathrm{d}t = F(t) $,
$ \mathrm{d}X/\mathrm{d}t = \mu(S,P)\,X - (F/V)\,X $,
$ \mathrm{d}S/\mathrm{d}t = -(1/Y_{X/S})\,\mu\,X - m_s\,X + (F/V)(S_{\mathrm{in}} - S) $,
$ \mathrm{d}P/\mathrm{d}t = q_p\,X - (F/V)\,P $,
$ \mathrm{d}C_L/\mathrm{d}t = k_L a\,(C^{\star} - C_L) - q_{O_2}\,X - (F/V)\,C_L $.

Specific growth is
$ \mu(S,P) = \mu_{\max}\,S / (K_S + S + S^{2}/K_I)\cdot(1 - P/P_{\mathrm{crit}})^{\gamma} $,
using the known `gamma` from `/workspace/data/reactor_specs.json`. At an unknown switch time $ t_{\mathrm{switch}} $, $ \mu_{\max} $ is multiplied by a factor $ \alpha \in [0.4, 0.7] $. Before the switch, product formation is growth-associated with $ q_p = Y_{P/X}\,\mu $ using `Y_px` from `/workspace/data/reactor_specs.json`; after the switch, $ q_p = \beta $. Oxygen uptake is $ q_{O_2} = Y_{X/O}\,\mu + m_o $ with `Y_xo` and `m_o` from the same file. The tabulated `F` column is mildly noisy and must not be treated as an exact control archive.

Estimate $ \mu_{\max} $, $ K_S $, $ K_I $, $ P_{\mathrm{crit}} $, $ t_{\mathrm{switch}} $, $ \alpha $, and $ \beta $ from the two telemetry files. Write `/workspace/output/estimated_parameters.json` as a JSON object that contains exactly the keys `mu_max`, `K_S`, `K_I`, `P_crit`, `t_switch`, `alpha`, and `beta`. Every value must be a finite number and must be non-negative, with rates in 1/h or g/(g·h), concentrations in g/L, and time in h. Values of `mu_max`, `P_crit`, `alpha`, and `beta` must lie within 10% relative error of the process truth implied by the data. Values of `K_S` and `K_I` must lie within 20% relative error of that truth. The value of `t_switch` must be within 0.5 h of truth.

Using those kinetics, design a piecewise-constant feed for Run 03 on $ t \in [0, t_{\mathrm{final}}] $, where `t_final` and `n_feed_intervals` come from `/workspace/data/reactor_specs.json` and `n_feed_intervals` equals 10 equal intervals. Initialize from `run03_initial_state` in `/workspace/data/reactor_specs.json`, using keys `X0`, `S0`, `P0`, `CL0`, and `V0`. Every feed rate must be non-negative and must not exceed `max_pump_rate` from `/workspace/data/reactor_specs.json`. Broth volume must remain at or below `V_max` from `/workspace/data/reactor_specs.json`. Throughout the Run 03 trajectory, substrate must satisfy $ S(t) \le 15 $ g/L and dissolved oxygen must satisfy $ C_L(t) \ge 0.15\,C^{\star} $, where $ C^{\star} $ is `C_star` in `/workspace/data/reactor_specs.json`. Final product mass $ V(t_f)\,P(t_f) $ must be at least 249.33 g. Write `/workspace/output/optimal_feed.json` containing the key `feed_rates_L_per_h` as an array of exactly ten finite, non-negative floats (L/h) in chronological interval order. The file may also include an `intervals` array of objects with keys `index`, `t_start`, `t_end`, and `F`.

Both `/workspace/output/estimated_parameters.json` and `/workspace/output/optimal_feed.json` must exist and must be valid JSON. All required numeric fields in those files must be finite.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
