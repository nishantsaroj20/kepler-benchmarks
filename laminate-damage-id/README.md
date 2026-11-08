# GFRP plate with a breathing midplane split

## Difficulty

The plate is a compressed eight-ply GFRP laminate with one midplane split. On the upward preparation the split is open: bending stiffness moves from the intact laminate toward the two four-ply sublaminates, and extensional stiffness is reduced by the open-membrane factor on the separated area. On the downward preparation the membrane stays intact, but bending keeps that full split change instead of returning to the intact plate. The separated fraction is the product of three compact elliptical zones. A trial zone that misses a real ellipse does not change the operators, and one ellipse cannot stand in for the other two.

The two unknown edge springs are shared by every hold. The closed preparation still carries bending damage, so those springs are not a damage-free reading of the downward records. A centre scan that freezes the springs at a nominal value drops the later ellipses on a compromise blob, and the modal residual stays well above the noise. Updating the springs on the patches already found is what moves the next scan onto a real patch. Four of seven modes are missing, the preparation label was not stored, and ellipse orientation trades against severity. Noise-level fits can span almost a half-turn in angle while the centres, the open-branch forecast and the closed-branch snap remain informative.

The forecast forces sit outside the held set. The snap load is the first stability loss of the closed branch, not an equal-energy jump or a linear buckling load. The records are seeded synthetic output of the nine-term shallow-plate model documented with the data. They are not measurements from a physical panel. Twenty continuous unknowns plus four preparation signs remain. A composite structural-health-monitoring engineer, or a structural-dynamics researcher, performs this work in practice when a compressed laminate panel must be assessed from incomplete vibration records.

## Reference solution

The reference builds the open and closed operators from the published breathing law, including the spatial membrane reduction on the upward preparation, and keeps the lower-residual preparation at each hold. Centres are scanned on a 12 mm grid, with patches not yet placed parked at a corner. After each accepted centre, a short bounded least-squares step updates that patch's axes, angle and severity together with both spring multipliers. Those updated springs are what make the next scan select a real patch rather than the compromise blob. A 4 mm local refinement follows. Two full polishes, one started with every angle shifted by 30 degrees, keep the lower residual. No authored parameter is a start.

The closed-branch snap is the last positive force at which a bisection still finds a stable negative equilibrium. The open-branch forecast is the lowest four frequencies of the positive equilibrium at the two requested forces.

```sh
python3 authoring/provenance/generate_data.py
LDI_DATA_DIR=environment/data LDI_OUT_FILE=/tmp/damage.json python3 solution/solve_damage.py
```

On the seeded records this search finishes at cost 4.37e-6, with the authored signs, a snap relative error of 1.67e-4 and a largest forecast relative error of 3.22e-4. The largest centre miss is 7.8 mm.

## Verification

The isolated checks read only the JSON artifact and `tests/ground_truth.json`. They do not import a solver. Zone centres, ellipse axes, severities, fixture multipliers and preparation signs in that file are the authored generator inputs. The snap load and the eight forecast frequencies are not taken from the reference condensation. `authoring/provenance/independent.py` integrates the membrane strains, including the open-branch membrane reduction, on the published 30 by 26 rule, eliminates the in-plane coefficients by a dense factorisation, and forms the tangent by central differences. A Newton continuation then brackets the closed-branch eigenvalue crossing. `generate_data.py` writes those independent values into the truth file.

On the authored parameters the direct model and the reference condensation differ by at most 7.8e-10 relative in the compared frequencies and by 5.3e-7 relative in the snap load. The stored reference inversion, checked against the direct model, is 3.22e-4 relative on the forecast and 1.67e-4 relative on the snap. A 40 by 36 quadrature in the same direct code moves the open-branch frequencies at -90 N by 5.1e-3 relative. That is a different rule from the published 30 by 26 model, so it is outside the contract. The environment image contains the model statement and the four records, not the generator seed or the truth file.

The discrimination bound is a second inverse, not a restart of the reference estimator. `authoring/provenance/independent_inverse.py` builds frequencies and sensor shapes from the direct Newton's tangent, assigns them to the four records, and polishes with its own bounded least squares. It does not import the reference condensation. Its JSON was scored by `tests/test_outputs.py` and passed. That run finished at cost 4.71e-6 with the authored signs, snap relative error 2.03e-4 and forecast relative error 2.04e-4. Its largest misses are 0.0078 m in centre, 0.0025 m in `a`, 0.0031 m in `b`, 66 degrees in orientation, 0.109 in severity, and 0.0007 and 0.0005 in the two multipliers.

Five further polishes of the reference estimator, from that fit, from the authored parameters, and from shifted angles, centres and severities, stay below a modal cost of 4.8e-6. Their largest forecast relative error is 3.72e-4 and their largest misses are 0.0099 m in centre, 0.0064 m in `a`, 0.0038 m in `b`, 89.9 degrees in orientation, 0.133 in severity, 0.033 in `x0_multiplier` and 0.0034 in `yL_multiplier`. The orientation span is the half-period, so a tighter angle cap would reject noise-level fits. The reporting caps sit above both clouds: centres 0.013 m, `a` 0.012 m, `b` 0.008 m, orientation 90 degrees modulo 180, severity 0.16, and each multiplier 0.07. Snap and forecast caps are 0.0008.

The midpoint of every search interval, centre (0.25 m, 0.20 m), axes 0.027 m and 0.014 m, angle 0 and severity 0.675, matches none of the three patches. Against the middle patch the y error is 0.014 m and the severity error is 0.165, outside 0.013 m and 0.16. The other two patches miss the centre cap. A polished compromise-blob placement has forecast relative error 1.70e-3 and centres more than 0.06 m off. An intact closed branch, or a membrane left whole, leaves forecast relative errors 0.00284 and 0.00319. Parking one real patch at a corner misses the centre cap. The counts are in `authoring/evidence/zone-calibration.json`.

The artifact must be a regular file of at most 1 MiB, with unique keys. `plates` is a JSON array of one object. That object has three in-bound zones (`a>b`), two in-bound multipliers, numeric branch signs and a positive snap load. The forecast must be two rows at the requested forces with four positive frequencies each. Branch signs must match. Each zone is matched, in any order, to the authored patch inside the caps above. The snap relative error and each forecast-frequency relative error must be at most 0.0008. A 1e-12 allowance is included.

```sh
LDI_ARTIFACT=/tmp/damage.json pytest tests/test_outputs.py -q
```
