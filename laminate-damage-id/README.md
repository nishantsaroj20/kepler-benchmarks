# GFRP plate with a breathing midplane split

## Difficulty

The plate is a compressed eight-ply GFRP laminate with one midplane split. On the upward preparation the split is open: bending moves toward the two four-ply sublaminates, and the membrane is reduced on the separated area. On the downward preparation the membrane stays intact, but bending keeps that full split change. The separated fraction is the product of four compact elliptical zones. A scan that stops after three patches leaves the fourth centre outside the cap.

The two unknown edge springs are shared by every hold, and closed-branch bending still carries the split, so a nominal-spring scan prefers a compromise blob. Four of seven modes are missing, and on one hold the fundamental is missing too. Ellipse angle still trades with area and severity, so those caps stay wide. The centre cap does not. The forecast forces sit outside the held set, and the snap is the first stability loss of the closed branch. Twenty-six continuous unknowns and four preparation signs remain.

## Reference solution

The reference builds both operators from the published breathing law and keeps the lower-residual preparation at each hold. Centres are scanned on a 12 mm grid, with patches not yet placed parked at a corner. After each accepted centre, a bounded least-squares step updates that patch together with both spring multipliers before the next scan. A 4 mm local refinement follows. Two full polishes, one started with every angle shifted by 30 degrees, keep the lower residual. No authored parameter is a start.

The closed-branch snap is the last positive force at which a bisection still finds a stable negative equilibrium. The open-branch forecast is the lowest four frequencies of the positive equilibrium at the two requested forces.

```sh
python3 authoring/provenance/generate_data.py
LDI_DATA_DIR=environment/data LDI_OUT_FILE=/tmp/damage.json python3 solution/solve_damage.py
```

On the seeded records this search finishes at cost 1.40e-5, with the authored signs. The largest centre misses are 18.6 mm and 14.5 mm. The snap relative error is 1.64e-3 and the largest forecast relative error is 1.18e-3.

## Verification

The isolated checks read only the JSON artifact and `tests/ground_truth.json`. Zone geometry, severities, fixture multipliers and preparation signs in that file are the authored generator inputs. The snap and the eight forecast frequencies come from `authoring/provenance/independent.py` on the published 30 by 26 rule. `generate_data.py` writes them into the truth file.

On the authored parameters the direct model and the reference condensation differ by at most 6.3e-10 relative in frequency and by 3.4e-7 relative in the snap. A polish of the same unweighted residual, started at the authored parameters, stays inside 1 mm on the centres and 0.01 on severity, with snap relative error 2.1e-4 and forecast relative error 1.4e-4. The sequential search above is farther out, and the caps sit above that search rather than above only the authored start: centres 0.026 m, `a` 0.022 m, `b` 0.010 m, severity 0.40, angle 90 degrees, and each multiplier 0.07. Snap and forecast caps are 0.0024. The angle cap stays at the half-period because the sequential search already misses one angle by 75.6 degrees. The severity cap stays wide because area still trades with severity; a 0.16 cap is what rejected noise-level inversions when only three patches were scored.

The midpoint of every search interval, centre (0.25 m, 0.20 m), matches none of the four patches. Against the nearest patch the y error is 0.036 m, outside the 0.026 m centre cap. Moving any one recovered patch to a corner fails that same cap. The counts are in `authoring/evidence/zone-calibration.json`.

The artifact must be a regular file of at most 1 MiB, with unique keys. `plates` is a JSON array of one object. That object has four in-bound zones (`a>b`), two in-bound multipliers, numeric branch signs and a positive snap load. The forecast must be two rows at the requested forces with four positive frequencies each. Branch signs must match. Each zone is matched, in any order, to the authored patch inside the caps above. The snap relative error and each forecast-frequency relative error must be at most 0.0024. A 1e-12 allowance is included.

```sh
LDI_ARTIFACT=/tmp/damage.json pytest tests/test_outputs.py -q
```
