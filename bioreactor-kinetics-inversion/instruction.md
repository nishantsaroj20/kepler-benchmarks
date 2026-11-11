A 10 L-class fed-batch bioreactor produced two historical telemetry files at `/workspace/data/sensor_run_01.csv` and `/workspace/data/sensor_run_02.csv`. Each CSV provides irregularly sampled columns `time_h`, `X`, `S`, `P`, `C_L`, `V`, and `F` (biomass, substrate, inhibitory product, dissolved oxygen, broth volume, and instantaneous feed). Known operating constants, yields, aeration parameters, pump limits, the product-inhibition exponent `gamma`, and the Run 03 initial state are given in `/workspace/data/reactor_specs.json`.

The culture follows Haldane–Andrews growth with product inhibition and a single metabolic phase shift. With state $ x = [X, S, P, C_L]^{\mathsf T} $ and volume $ V $,

$ \mathrm{d}V/\mathrm{d}t = F(t) $,
$ \mathrm{d}X/\mathrm{d}t = \mu(S,P)\,X - (F/V)\,X $,
$ \mathrm{d}S/\mathrm{d}t = -(1/Y_{X/S})\,\mu\,X - m_s\,X + (F/V)(S_{\mathrm{in}} - S) $,
$ \mathrm{d}P/\mathrm{d}t = q_p\,X - (F/V)\,P $,
$ \mathrm{d}C_L/\mathrm{d}t = k_L a\,(C^{\star} - C_L) - q_{O_2}\,X - (F/V)\,C_L $.

Specific growth is
$ \mu(S,P) = \mu_{\max}\,S / (K_S + S + S^{2}/K_I)\cdot(1 - P/P_{\mathrm{crit}})^{\gamma} $
with the known `gamma` from `/workspace/data/reactor_specs.json`. At an unknown switch time $ t_{\mathrm{switch}} $, $ \mu_{\max} $ is multiplied by a factor $ \alpha \in [0.4, 0.7] $. Before the switch, product formation is growth-associated with $ q_p = Y_{P/X}\,\mu $ using $ Y_{P/X} $ from the specs file; after the switch, $ q_p = \beta $ (non-growth-associated). Oxygen uptake follows $ q_{O_2} = Y_{X/O}\,\mu + m_o $ with the coefficients in the same file. Tabulated feed $ F $ is mildly noisy and must not be treated as an exact control archive.

Estimate $ \mu_{\max} $, $ K_S $, $ K_I $, $ P_{\mathrm{crit}} $, $ t_{\mathrm{switch}} $, $ \alpha $, and $ \beta $ from the two telemetry runs. Write `/workspace/output/estimated_parameters.json` as a JSON object containing exactly the keys `mu_max`, `K_S`, `K_I`, `P_crit`, `t_switch`, `alpha`, and `beta` with finite numeric values (rates in 1/h or g/(g·h), concentrations in g/L, time in h). Values of `mu_max`, `P_crit`, `alpha`, and `beta` must lie within 10% relative error of the process truth implied by the data; `K_S` and `K_I` must lie within 20% relative error; and `t_switch` must be within 0.5 h of truth.

Using the estimated kinetics, design a piecewise-constant feed for Run 03 over $ t \in [0, t_{\mathrm{final}}] $ with `t_final` and `n_feed_intervals` from `/workspace/data/reactor_specs.json` (ten equal intervals). Initialize from `run03_initial_state` in that file. Feed rates must be non-negative and must not exceed `max_pump_rate`. Broth volume must remain at or below `V_max`. Along the resulting trajectory under the true process, substrate must satisfy $ S(t) \le 15 $ g/L and dissolved oxygen must satisfy $ C_L(t) \ge 0.15\,C^{\star} $. The schedule must achieve final product mass $ V(t_f)\,P(t_f) $ of at least 249.33 g. Write `/workspace/output/optimal_feed.json` containing key `feed_rates_L_per_h` as an array of exactly ten non-negative finite floats (L/h) in chronological interval order. An `intervals` array of objects with `index`, `t_start`, `t_end`, and `F` may also be included.

Both output files must exist, use valid JSON, and contain only finite numbers for all graded fields.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
