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

On the authored parameters the direct model and the reference condensation differ by at most 6.3e-10 relative in frequency and by 3.4e-7 relative in the snap. A polish of the same unweighted residual, started at the authored parameters, stays inside 1 mm on the centres and 0.01 on severity, with snap relative error 2.1e-4 and forecast relative error 1.4e-4.

Every cap is set from two inverses that do not share a residual, and from a case that must fail when one was measured. The reference search retries major axis, severity and angle. The direct Newton inverse in `authoring/provenance/independent_inverse.py` integrates strains on the 30 by 26 rule, does not import the reference residual, and does not start from the authored parameters. `pytest tests/test_outputs.py` passed on that JSON. Angle acceptance is 90 degrees because a 180-degree flip is the same major axis; both inverses are inside that half-period. Severity acceptance stays at 0.40 because the direct inverse misses one severity by 0.195, and a 0.16 cap rejected noise-level inversions when only three patches were scored. Semi-axis caps sit above both inverses because area still trades with severity. A multiplier of 0.07 rejects a nominal value of 1, which is 0.24 away from 1.24 and from 0.76, while both inverses are within 0.009. Snap acceptance of 0.00115 and forecast acceptance of 0.00080 sit above both inverses and below the continuous polish that does not retry the ellipse. That polish is at snap relative error 1.46e-3 and forecast relative error 9.85e-4, and the same tests reject it on the snap. The 1e-12 allowance is only for inclusive endpoints and for equality of a reported forecast force with the requested force.

| quantity | reference | direct Newton | cap | outside the cap |
|---|---|---|---|---|
| centre x | 0.0158 m | 0.0099 m | 0.026 m | a patch moved to a corner |
| centre y | 0.0125 m | 0.0138 m | 0.026 m | interval midpoint, 0.036 m off |
| semi-axis a | 0.0165 m | 0.0145 m | 0.022 m | above both inverses |
| semi-axis b | 0.0058 m | 0.0060 m | 0.010 m | above both inverses |
| angle | 55.4 degrees | 49.1 degrees | 90 degrees | half-period identity |
| severity | 0.056 | 0.195 | 0.40 | a 0.16 cap rejected noise-level inversions |
| x = 0 spring | 0.0083 | 0.0089 | 0.07 | nominal 1 misses 1.24 by 0.24 |
| y = Ly spring | 0.0006 | 0.0018 | 0.07 | nominal 1 misses 0.76 by 0.24 |
| snap, relative | 0.000257 | 0.000873 | 0.00115 | continuous polish 0.00146 |
| forecast, relative | 0.000200 | 0.000622 | 0.00080 | continuous polish 0.000985 |

The midpoint of every search interval, centre (0.25 m, 0.20 m), matches none of the four patches. The same tests reject that plate on the zone match, and they reject a corner-moved patch on the centre cap. The reference-search counts are in `authoring/evidence/zone-calibration.json`, and the direct-inverse artifact plus both pytest summaries are in `authoring/evidence/independent-inverse.json`.

The artifact must be a regular file of at most 1 MiB, with unique keys. `plates` is a JSON array of one object. That object has four in-bound zones (`a>b`), two in-bound multipliers, numeric branch signs and a positive snap load. The forecast must be two rows at the requested forces with four positive frequencies each. Branch signs must match. Each zone is matched, in any order, to the authored patch inside the caps above. The snap relative error must be at most 0.00115, and each forecast-frequency relative error must be at most 0.00080. A 1e-12 allowance is included.

```sh
LDI_ARTIFACT=/tmp/damage.json pytest tests/test_outputs.py -q
```
