# GFRP plate with a breathing midplane split

## Difficulty

The plate is a compressed eight-ply GFRP laminate with one midplane split. On the upward preparation the split is open: bending moves toward the two four-ply sublaminates, and the membrane is reduced on the separated area. On the downward preparation the membrane stays intact, but bending keeps that full split change. The separated fraction is the product of four compact elliptical zones. A scan that stops after three patches leaves the fourth centre outside the cap.

The two unknown edge springs are shared by every hold, and closed-branch bending still carries the split, so a nominal-spring scan prefers a compromise blob. Four of seven modes are missing, and on one hold the fundamental is missing too. Ellipse angle still trades with area and severity, so those caps stay wide. The centre cap does not. A polish that keeps the angle, area and severity it first accepted misses the snap limit and the forecast limit. The forecast forces sit outside the held set, and the snap is the first stability loss of the closed branch. Twenty-six continuous unknowns and four preparation signs remain.

The four holds are synthetic. They are generated from the published breathing-split law, then reduced to four modes out of seven with the preparation label omitted. That is the incompleteness a laboratory modal survey of a compressed panel actually leaves, so the identification is at a realistic measurement difficulty rather than a complete noise-free spectrum. A composites structural-dynamics engineer in nondestructive evaluation does this work as part of the job: recover the midplane splits and the unknown edge springs from those holds, then forecast the closed-branch snap and the open-branch frequencies at forces the test stand did not apply.

## Reference solution

The reference builds both operators from the published breathing law and keeps the lower-residual preparation at each hold. Centres are scanned on a 12 mm grid, with patches not yet placed parked at a corner. After each accepted centre, a bounded least-squares step updates that patch together with both spring multipliers before the next scan. A 4 mm local refinement follows. Two full polishes, one started with every angle shifted by 30 degrees, keep the lower residual. Each patch is then tried on a grid of major axes, severities and angles, and that lower residual is polished again. No authored parameter is a start.

The closed-branch snap is the last positive force at which a bisection still finds a stable negative equilibrium. The open-branch forecast is the lowest four frequencies of the positive equilibrium at the two requested forces.

```sh
python3 authoring/provenance/generate_data.py
LDI_DATA_DIR=environment/data LDI_OUT_FILE=/tmp/damage.json python3 solution/solve_damage.py
```

On the seeded records this search finishes at cost 4.24e-6, with the authored signs. The largest centre misses are 15.8 mm and 12.5 mm. The snap relative error is 2.57e-4 and the largest forecast relative error is 2.00e-4.

## Verification

The isolated checks read only the JSON artifact and `tests/ground_truth.json`. Zone geometry, severities, fixture multipliers and preparation signs in that file are the authored generator inputs. The snap and the eight forecast frequencies come from `authoring/provenance/independent.py` on the published 30 by 26 rule. `generate_data.py` writes them into the truth file.

On the authored parameters the direct model and the reference condensation differ by at most 6.3e-10 relative in frequency and by 3.4e-7 relative in the snap. A polish of the same unweighted residual, started at the authored parameters, stays inside 1 mm on the centres and 0.01 on severity, with snap relative error 2.1e-4 and forecast relative error 1.4e-4. The geometry caps stay wide enough for both of the inverses below, not only for a polish started at the authored parameters: centres 0.026 m, `a` 0.022 m, `b` 0.010 m, severity 0.40, angle 90 degrees, and each multiplier 0.07. The snap cap is 0.00115 and the forecast cap is 0.00080. The angle cap stays at the half-period because the reference still misses one angle by 55 degrees. The severity cap stays wide because the direct inverse misses one severity by 0.195, and a 0.16 cap is what rejected noise-level inversions when only three patches were scored. A continuous polish that does not retry major axis, severity and angle finishes at snap relative error 1.46e-3 and forecast relative error 9.85e-4. `pytest tests/test_outputs.py` rejects that artifact on the snap.

A second inverse, `authoring/provenance/independent_inverse.py`, does not import the reference condensation, that condensation's residual, or its least-squares fit. It integrates strains on the 30 by 26 Gauss rule, continues equilibria by Newton, and matches each hold on that direct model's frequencies and shapes. The search starts from corner-parked patches, not from the authored parameters. The real verifier scored the JSON it wrote: `LDI_ARTIFACT` pointed at that file and `pytest tests/test_outputs.py -q` reported 3 passed. The misses below sit inside the caps. The reference misses a centre by 15.8 mm and an angle by 55 degrees; this direct inverse misses a centre by 13.8 mm, an angle by 49 degrees and a severity by 0.195.

| quantity | direct Newton inverse | acceptance cap |
|---|---|---|
| centre x | 0.0099 m | 0.026 m |
| centre y | 0.0138 m | 0.026 m |
| semi-axis a | 0.0145 m | 0.022 m |
| semi-axis b | 0.0060 m | 0.010 m |
| angle | 49.1 degrees | 90 degrees |
| severity | 0.195 | 0.40 |
| x = 0 spring | 0.0089 | 0.07 |
| y = Ly spring | 0.0018 | 0.07 |
| snap, relative | 0.000873 | 0.00115 |
| forecast, relative | 0.000622 | 0.00080 |

The midpoint of every search interval, centre (0.25 m, 0.20 m), matches none of the four patches. Against the nearest patch the y error is 0.036 m, outside the 0.026 m centre cap. The same verifier rejects that plate: `test_preparations_snap_and_open_branch_forecast` fails and the other two tests pass. Moving any one recovered patch to a corner fails that same centre cap. The reference-search counts are in `authoring/evidence/zone-calibration.json`, and the direct-inverse artifact plus both pytest summaries are in `authoring/evidence/independent-inverse.json`.

The artifact must be a regular file of at most 1 MiB, with unique keys. `plates` is a JSON array of one object. That object has four in-bound zones (`a>b`), two in-bound multipliers, numeric branch signs and a positive snap load. The forecast must be two rows at the requested forces with four positive frequencies each. Branch signs must match. Each zone is matched, in any order, to the authored patch inside the caps above. The snap relative error must be at most 0.00115, and each forecast-frequency relative error must be at most 0.00080. A 1e-12 allowance is included.

```sh
LDI_ARTIFACT=/tmp/damage.json pytest tests/test_outputs.py -q
```
