Put your solver at `/workspace/output/column_solver.py`. It will be called repeatedly on case folders; each call must finish in 120 seconds on a single CPU and must work without network access (held-out runs are offline even though the agent environment may have public networking).

Sample cases live under `/workspace/data/sample_cases/` and are indexed in `/workspace/data/run_index.json`. Look through those first: some beds are densely packed, some runs are fast with strong axial dispersion, some mid-bed traces have dropouts or spikes, some feeds are concentrated, and one isotherm is nearly linear. Held-out runs use other cases built the same way, just with different numbers.

A case folder has `case.json`, `breakthrough_exp_A.csv`, and `breakthrough_exp_B.csv`. Run like this:

`python3 /workspace/output/column_solver.py --case /absolute/path/to/case.json --out /absolute/path/to/result.json`

Paths listed under the `csv` keys are relative to that folder. Create the `--out` file before you exit.

In `case.json` the known packing values are `bed_voidage`, `particle_density_kg_m3`, and `lab_bed_length_m`. Blocks `experiment_A` and `experiment_B` each give `csv`, `c_in_mol_m3`, `interstitial_velocity_m_s`, and `t_end_s`. The design block `design_duty` has `c_in_mol_m3`, `interstitial_velocity_m_s`, `service_time_s`, and `outlet_mole_fraction_limit`. The CSVs have `time_s`, `c_outlet_mol_m3`, and `c_midbed_mol_m3`. If the mid-bed column is junk in places, lean on the outlet.

Model: 1-D column, Langmuir equilibrium, LDF kinetics, axial dispersion. `interstitial_velocity_m_s` is the $ u $ in the balances below. Voidage $ \varepsilon $, particle density $ \rho_p $, dispersion $ D_{\mathrm{ax}} $, solid loading $ q $:

$$
\frac{\partial c}{\partial t} + u\frac{\partial c}{\partial z} = D_{\mathrm{ax}}\frac{\partial^{2}c}{\partial z^{2}} - \frac{1-\varepsilon}{\varepsilon}\rho_p\frac{\partial q}{\partial t}
$$

$$
\frac{\partial q}{\partial t} = k_{\mathrm{LDF}}\bigl(q^{\star}(c) - q\bigr)
$$

$$
q^{\star}(c) = \frac{q_m K c}{1 + K c}
$$

Inlet is Dirichlet $ c(0,t)=c_{\mathrm{in}} $. Outlet is zero Neumann. Start from a clean bed.

From the two breakthroughs, get $ q_m $, $ K $, $ k_{\mathrm{LDF}} $, and $ D_{\mathrm{ax}} $. Then find the shortest bed length $ L_{\min} $ for `design_duty` such that the outlet never exceeds `outlet_mole_fraction_limit` times that duty's `c_in_mol_m3` on $ [0,\mathrm{service\_time\_s}] $.

`--out` must be JSON with keys `qm`, `K`, `k_LDF`, `D_ax`, and `L_min_m` (you may also write `case_id`). All five required values must be finite and positive. Units: `qm` mol/kg, `K` m^3/mol, `k_LDF` 1/s, `D_ax` m^2/s, `L_min_m` m. Stay within 12% of truth on `qm`, 25% on `K`, 50% on `k_LDF`, and 45% on `D_ax`. The length you report has to actually meet the purity duty on the true process and should not be an obviously padded giant bed.

BLAS/OpenMP are forced to one thread. The solver file must be a real, non-empty `.py` that `python3` can run.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
