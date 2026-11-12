Write a Python program at `/workspace/output/column_solver.py` that identifies Langmuir–LDF–dispersion parameters from two fixed-bed breakthrough experiments and computes a minimum design packed length. The verifier re-runs that program on held-out cases. Each run of `/workspace/output/column_solver.py` must finish within 120 seconds on 1 CPU.

Visible sample cases are listed in `/workspace/data/run_index.json` under `/workspace/data/sample_cases/`. Those samples already include every deciding condition class used later: a clean baseline, dense packing (low voidage / high particle density), a high-Peclet dispersion regime, mid-bed sensor dropouts and spikes, elevated feed concentrations, and a near-linear (small-$ K $) isotherm. Held-out evaluation reuses those same condition classes with different numeric values.

Every case directory contains `case.json`, `breakthrough_exp_A.csv`, and `breakthrough_exp_B.csv`. Invoke the program as
`python3 /workspace/output/column_solver.py --case /absolute/path/to/case.json --out /absolute/path/to/result.json`.
CSV paths inside `case.json` are relative to that case directory. The program must create the `--out` JSON file before exiting.

Each `case.json` provides known packing fields `bed_voidage`, `particle_density_kg_m3`, and `lab_bed_length_m`, plus `experiment_A` and `experiment_B` with keys `csv`, `c_in_mol_m3`, `interstitial_velocity_m_s`, and `t_end_s`, and a `design_duty` object with keys `c_in_mol_m3`, `interstitial_velocity_m_s`, `service_time_s`, and `outlet_mole_fraction_limit`. Breakthrough CSVs contain columns `time_s`, `c_outlet_mol_m3`, and `c_midbed_mol_m3`. When mid-bed values are missing or spiked, treat the outlet series as the primary identification signal.

The packed bed obeys one-dimensional convection–dispersion with a Langmuir isotherm and linear-driving-force uptake. The field `interstitial_velocity_m_s` is the advection speed $ u $ in the PDE below (interstitial velocity). With bed voidage $ \varepsilon $, particle density $ \rho_p $, axial dispersion $ D_{\mathrm{ax}} $, and solid loading $ q $,

$ \partial c/\partial t + u\,\partial c/\partial z = D_{\mathrm{ax}}\,\partial^{2}c/\partial z^{2} - ((1-\varepsilon)/\varepsilon)\,\rho_p\,\partial q/\partial t $,
$ \partial q/\partial t = k_{\mathrm{LDF}}\bigl(q^{\star}(c) - q\bigr) $,
$ q^{\star}(c) = q_m K c / (1 + K c) $.

Use a Dirichlet inlet $ c(0,t)=c_{\mathrm{in}} $ and a zero-Neumann outlet. Experiments start from a clean bed.

For every case the program must estimate $ q_m $, $ K $, $ k_{\mathrm{LDF}} $, and $ D_{\mathrm{ax}} $, then compute a minimum packed length $ L_{\min} $ for that case's `design_duty` so the outlet stays at or below `outlet_mole_fraction_limit` times the design `c_in_mol_m3` throughout $ [0, \mathrm{service\_time\_s}] $. Any standard numerical solution of this model that recovers the parameters and a feasible $ L_{\min} $ is acceptable; the task is not about matching one proprietary discretization.

Write a JSON object to `--out` containing exactly the keys `qm`, `K`, `k_LDF`, `D_ax`, and `L_min_m` (and optionally `case_id`). Every required value must be a finite positive number, with `qm` in mol/kg, `K` in m^3/mol, `k_LDF` in 1/s, `D_ax` in m^2/s, and `L_min_m` in m. Relative to the process truth for that case, `qm` must lie within 12%, `K` within 25%, `k_LDF` within 50%, and `D_ax` within 45%. The reported `L_min_m` must keep the true-model outlet within the design purity limit over the service window and must not be trivially oversized.

Threading environment variables are pinned to one thread during evaluation. A run that exceeds the 120-second limit fails the task.

`/workspace/output/column_solver.py` must exist, must be a regular non-empty Python file, and must be runnable via `python3`.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
