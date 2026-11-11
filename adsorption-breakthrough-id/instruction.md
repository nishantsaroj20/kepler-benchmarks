Two noisy fixed-bed breakthrough records are provided at `/workspace/data/breakthrough_exp_A.csv` and `/workspace/data/breakthrough_exp_B.csv`. Each CSV contains columns `time_s`, `c_outlet_mol_m3`, and `c_midbed_mol_m3` (mobile-phase solute concentration at the bed outlet and at mid-bed). Known packing constants, the two experiment operating conditions, and a held-out design duty are given in `/workspace/data/bed_specs.json`, including keys `bed_voidage`, `particle_density_kg_m3`, `lab_bed_length_m`, `experiment_A`, `experiment_B`, and `design_duty`.

The packed bed obeys one-dimensional convection–dispersion with a Langmuir isotherm and linear-driving-force (LDF) solid uptake. With interstitial superficial velocity $ u $, bed voidage $ \varepsilon $, particle density $ \rho_p $, axial dispersion $ D_{\mathrm{ax}} $, and solid loading $ q $,

$ \partial c/\partial t + u\,\partial c/\partial z = D_{\mathrm{ax}}\,\partial^{2}c/\partial z^{2} - ((1-\varepsilon)/\varepsilon)\,\rho_p\,\partial q/\partial t $,
$ \partial q/\partial t = k_{\mathrm{LDF}}\bigl(q^{\star}(c) - q\bigr) $,
$ q^{\star}(c) = q_m K c / (1 + K c) $.

Use a Dirichlet inlet $ c(0,t)=c_{\mathrm{in}} $ and a zero-Neumann outlet. Both historical experiments start from a clean bed. Values of `bed_voidage`, `particle_density_kg_m3`, and `lab_bed_length_m` in `/workspace/data/bed_specs.json` are known. Inlet concentration and superficial velocity for each experiment are given under `experiment_A` and `experiment_B` as `c_in_mol_m3` and `superficial_velocity_m_s`.

Estimate the hidden parameters $ q_m $, $ K $, $ k_{\mathrm{LDF}} $, and $ D_{\mathrm{ax}} $ using both breakthrough files. Using those estimates, compute the minimum packed length $ L_{\min} $ for the held-out duty in `design_duty` of `/workspace/data/bed_specs.json`. That duty specifies `c_in_mol_m3`, `superficial_velocity_m_s`, `service_time_s`, and `outlet_mole_fraction_limit`. The design bed must keep the outlet concentration at or below `outlet_mole_fraction_limit` times the design `c_in_mol_m3` for all times in $ [0, \texttt{service\_time\_s}] $, starting from a clean bed under the same transport model.

Write `/workspace/output/adsorption_identification.json` as a JSON object that contains exactly the keys `qm`, `K`, `k_LDF`, `D_ax`, and `L_min_m`. Every value must be a finite positive number, with `qm` in mol/kg, `K` in m^3/mol, `k_LDF` in 1/s, `D_ax` in m^2/s, and `L_min_m` in m. Relative to the process truth implied by the data, `qm` must lie within 10%, `K` within 12%, `k_LDF` within 20%, `D_ax` within 30%, and `L_min_m` within 10%. The reported `L_min_m` must not be undersized relative to that truth and must keep the true-model outlet within the design purity limit over the service window.

`/workspace/output/adsorption_identification.json` must exist and must be valid JSON. All required numeric fields in that file must be finite.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
