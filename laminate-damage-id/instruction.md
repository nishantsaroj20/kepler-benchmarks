A compressed eight-ply GFRP panel has a breathing split on the midplane between ply 4 and ply 5. Separated area sits in three compact elliptical patches. The nominal rotational springs are given in the specification; the multipliers on the edges x = 0 and y = Ly are the unknown ones. Under `/app/data/` there are four vibration holds. Each hold kept four modes out of the lowest seven, and the preparation label was never written down.

`/app/data/plate_specs.json` lists the layup, springs, imperfection, precompression, the four hold forces, the two forecast forces, the noise scales and the search intervals. `/app/data/README.md` states the potential, the ellipse profile, the nine-mode space and the branch sign `q11`.

Recover all three patches, both spring multipliers and the preparation sign of each hold. That sign is the sign of the hold's q11. +1 is the upward open preparation and -1 is the downward closed preparation. Also report one closed-branch snap: the first positive force at which the downward preparation loses stability. Report the lowest four natural frequencies of the upward preparation at every force named under `forecast`. Those two forecast forces were not held. Any equilibrium method is acceptable.

Save `/app/output/damage.json`. Keep it at or below 1 MiB, and do not use a symbolic link. The text is one JSON object, and each key appears once inside its object. `plates` is a JSON array of exactly one object, and that object's `plate_id` is `"plate_01"`. That object carries:
- `zones`: a JSON array of exactly three objects, in any order. Each object has `xc`, `yc`, `a`, `b`, `theta_deg` and `d`.
- `fixture`: a JSON object with `x0_multiplier` and `yL_multiplier`. `x0_multiplier` scales the rotational spring on the edge x = 0, and `yL_multiplier` scales the rotational spring on the edge y = Ly.
- `branches`: a JSON object. Its keys are the four record ids from the specification, and each value is the number `-1` or `1`.
- `snap_load_N`: the closed-branch limit load in newtons, and it is positive.
- `open_branch_forecast_Hz`: a JSON array of exactly two objects, in `forecast.forces_N` order. Each object has `force_N` equal to the requested force and `frequencies_Hz`, a JSON array of four positive frequencies sorted from low to high.

Physical entries are JSON numbers. Booleans and quoted strings are not numbers.

Keep every reported parameter inside its interval, including the endpoints: `xc` 0.06 to 0.44 m, `yc` 0.055 to 0.345 m, `a` 0.014 to 0.040 m, `b` 0.008 to 0.020 m, `theta_deg` -90.0 to 90.0, `d` 0.40 to 0.95, and both multipliers 0.6 to 1.6. For each patch, `a` is greater than `b`.

Match the three patches to the physical patches in any order. For each matched patch, `xc` is within 0.013 m of that patch's `xc`, and `yc` is within 0.013 m of that patch's `yc`, each on its own. Semi-axis `a` is within 0.012 m of that patch's `a`, semi-axis `b` is within 0.008 m of that patch's `b`, and `d` is within 0.55 of that patch's `d`. `theta_deg` is within 90 degrees of that patch's angle once a 180-degree flip of the major axis is the same axis. Each spring multiplier is within 0.07 of its physical value. The four signs are the preparations that were held. `snap_load_N` and each of the eight forecast frequencies are within 0.0008 of the true value, with the difference divided by that true value. Every limit includes its endpoint. Another 1e-12 is allowed for roundoff.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
