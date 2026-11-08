Your goal is to identify two unordered compact elliptical zones of a breathing midplane split and two boundary mounting multipliers for a single symmetric eight-ply laminated GFRP plate, determine the static preparation branch sign for each of four modal observation records, predict the critical positive snap-through stability-loss force on the negative branch, and forecast the lowest four natural frequencies of the positive preparation at two forces that were not held during the test. How you formulate the nonlinear equilibrium path, solve the modal inverse problem, and compute the limit point is entirely up to you.

All model specifications, load values, search bounds, and baseline properties are provided under `/app/data/`. The file `/app/data/plate_specs.json` defines laminate material constants, geometry, stacking sequence, initial compressive stress resultants, known imperfection amplitudes, transverse hold forces, the two forecast forces, record identifiers, parameter intervals, and synthetic noise scales. Full documentation of the variational potential, the breathing split, the function spaces, the boundary springs, and the observation conventions is in `/app/data/README.md`.

`snap_load_N` is the first positive force at which the negative preparation loses stability. `open_branch_forecast_Hz` contains the lowest four ascending natural frequencies of the positive preparation at each force listed under `forecast` in `/app/data/plate_specs.json`, in that order. The breathing law, the function spaces, the observation conventions, and both of those definitions are specified in `/app/data/README.md`.

Write your results to `/app/output/damage.json` as a regular file (not a symlink) of at most 1 MiB. It must be a valid JSON object with unique keys at every level, whose root key `plates` contains a list with exactly one object having `plate_id` equal to `"plate_01"`. All physical values must be finite JSON numbers, not booleans or strings. This object must contain:
- `zones`: an array of exactly two zone objects, each with `xc`, `yc`, `a`, `b`, `theta_deg`, and `d`.
- `fixture`: an object with `x0_multiplier` and `yL_multiplier`.
- `branches`: an object mapping exactly the four record IDs in the specification to numeric branch signs `-1` or `1`.
- `snap_load_N`: a positive number, the closed-branch snap-through limit load in newtons.
- `open_branch_forecast_Hz`: a list of two objects in the same order as `forecast.forces_N`, each with `force_N` equal to that requested force and `frequencies_Hz`, a list of four positive frequencies in ascending order.

Reported parameters must lie in the specified bounds: `xc` in [0.06, 0.44] m, `yc` in [0.055, 0.345] m, `a` in [0.026, 0.055] m, `b` in [0.010, 0.024] m, `theta_deg` in [-90.0, 90.0] degrees, `d` in [0.40, 0.95], and both fixture multipliers in [0.6, 1.6]. Each zone must have $ a > b $.

All four branch signs must match the physical preparation states. The relative error in `snap_load_N` must be at most 0.0005. Each of the eight forecast frequencies must be within a relative error of 0.0005 of the corresponding physical open-branch frequency. Every continuous bound and error limit includes its endpoint, with a roundoff allowance of $ 10^{-12} $.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
