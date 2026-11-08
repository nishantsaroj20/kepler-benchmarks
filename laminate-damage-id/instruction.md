Your goal is to identify two unordered, compact elliptical zones of a breathing midplane split and two boundary mounting stiffness multipliers for a laminated glass-fibre reinforced polymer (GFRP) plate, determine the preparation branch sign for each of four modal observation records, compute the negative-branch snap-through limit load, and predict the lowest four positive-preparation natural frequencies at two prescribed forecast forces, writing your final deliverable to `/app/output/damage.json`. How you formulate the forward model, perform the modal inversion, and execute the nonlinear continuation is entirely up to you.

All model specifications, initial conditions, material properties, parameter bounds, hold forces, and record identifiers are provided in `/app/data/plate_specs.json`. The governing equations, nonlinear breathing contact formulation, branch preparation and equilibrium continuation definitions, sensor coordinates, and observation conventions are detailed in `/app/data/README.md`. Refer directly to these files for all authoritative equations, baseline dimensions, and numerical conventions.

The deliverable must be a regular file (not a symlink) of at most 1 MiB in size. It must be a valid JSON object with unique keys at every level, whose root key `plates` maps to an array containing exactly one object with `plate_id` equal to `"plate_01"`. All physical quantities must be finite JSON numbers, not booleans or strings. The plate object must contain the following keys:
- `zones`: an array of exactly two damage zone objects, each containing finite numeric values for `xc`, `yc`, `a`, `b`, `theta_deg`, and `d`.
- `fixture`: an object containing finite numeric values for `x0_multiplier` and `yL_multiplier`.
- `branches`: an object mapping each of the four record IDs defined in `plate_specs.json` to its numeric preparation branch sign, `-1` or `1`.
- `snap_load_N`: a finite positive number representing the snap-through stability-loss force in newtons on the negative branch.
- `open_branch_forecast_Hz`: an array of exactly two objects ordered to match `forecast.forces_N`, where each object contains `force_N` matching the prescribed load within 10^-12 and `frequencies_Hz`, a list of the lowest four positive natural frequencies in ascending order.

Reported parameters must lie in the specified bounds: `xc` in [0.06, 0.44] m, `yc` in [0.055, 0.345] m, `a` in [0.026, 0.055] m, `b` in [0.010, 0.024] m, `theta_deg` in [-90.0, 90.0] degrees, `d` in [0.40, 0.95], and both fixture multipliers in [0.6, 1.6]. Each zone must have a > b.

All four branch signs must match the physical preparation states. The relative error in `snap_load_N` must be at most 0.0005. Each of the eight forecast frequencies must be within a relative error of 0.0005 of the corresponding physical open-branch frequency. Every continuous bound and error limit includes its endpoint, with a roundoff allowance of 1e-12.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
