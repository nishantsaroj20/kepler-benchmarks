A compressed eight-ply GFRP panel has a breathing split on the midplane between ply 4 and ply 5. Separated area sits in three compact elliptical patches. The nominal rotational springs are given in the specification; the multipliers on the edges x = 0 and y = Ly are the unknown ones. Under `/app/data/` there are four vibration holds. Each hold kept four modes out of the lowest seven, and the preparation label was never written down.

`/app/data/plate_specs.json` lists the layup, springs, imperfection, precompression, the four hold forces, the two forecast forces, the noise scales and the search intervals. `/app/data/README.md` states the potential, the ellipse profile, the nine-mode space and the branch sign `q11`.

Recover all three patches, both spring multipliers and the preparation sign of each hold. That sign is the sign of the hold's q11. +1 is the upward open preparation and -1 is the downward closed preparation. Also report one closed-branch snap: the first positive force at which the downward preparation loses stability. Report the lowest four natural frequencies of the upward preparation at every force named under `forecast`. Those two forecast forces were not held. Any equilibrium method is acceptable.

Save `/app/output/damage.json`. Keep it at or below 1 MiB, and do not use a symbolic link. The text is one JSON object, and each key appears once. `plates` is a JSON array of exactly one object, and that object's `plate_id` is `"plate_01"`. That object carries:
- `zones`: three objects, in any order, with `xc`, `yc`, `a`, `b`, `theta_deg` and `d`.
- `fixture`: `x0_multiplier` and `yL_multiplier`.
- `branches`: the four record ids from the specification, each set to the number `-1` or `1`.
- `snap_load_N`: the closed-branch limit load in newtons, and it is positive.
- `open_branch_forecast_Hz`: two objects in `forecast.forces_N` order. Each has `force_N` equal to the requested force and `frequencies_Hz`, four positive frequencies sorted from low to high.

Physical entries are JSON numbers. A boolean or a quoted string in those positions is rejected.

Keep every reported parameter inside its interval, including the endpoints: `xc` 0.06 to 0.44 m, `yc` 0.055 to 0.345 m, `a` 0.014 to 0.040 m, `b` 0.008 to 0.020 m, `theta_deg` -90.0 to 90.0, `d` 0.40 to 0.95, and both multipliers 0.6 to 1.6. For each patch, `a` is greater than `b`.

Match the three patches to the physical patches in any order. For each matched patch, `xc` must lie within 0.013 m and `yc` must lie within 0.013 m, semi-axis `a` within 0.012 m, semi-axis `b` within 0.008 m, and `d` within 0.16. `theta_deg` must lie within 90 degrees once a 180-degree flip of the major axis is treated as the same axis. Each spring multiplier must lie within 0.07 of its physical value. The four signs must be the preparations that were held. `snap_load_N`, and each of the eight forecast frequencies, must come within 0.0008 of the true value when the difference is divided by that true value. Every limit includes its endpoint. Another 1e-12 is allowed for roundoff.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
