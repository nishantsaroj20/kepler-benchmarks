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

Severity is the flat direction of these holds. With the reference centres and axes held fixed, the weak patch at angle 60 degrees and severity 0.45 has noise-weighted chi-square about 127, against about 119 at the reference fit. The weights are the published scales, a frequency standard deviation of 8e-5 and a shape-component standard deviation of 0.0003, so this residual is not the unweighted reference residual. That alternative was scored by `tests/test_outputs.py` and passed. Its snap relative error is 3.4e-5 and its forecast relative error is 2.3e-4. Its severity is 0.46 from the authored value 0.91. The severity window is 0.55, which accepts every severity inside the search interval. A window of 0.16 would reject this noise-level plate.

`authoring/provenance/independent_inverse.py` builds the same records from the direct Newton tangent instead of the reference condensation, and its JSON passed the same tests: cost 4.71e-6, snap relative error 2.03e-4, forecast relative error 2.04e-4, centre miss 0.0078 m, severity miss 0.109. Five restarts of the reference residual stay within forecast 3.72e-4, centre 0.0099 m and severity 0.133, and they already span 89.9 degrees, so the angle cap stays at the 90-degree half-period. The other caps sit above that cloud: centres 0.013 m, `a` 0.012 m, `b` 0.008 m, and each multiplier 0.07. Snap and forecast caps are 0.0008.

The midpoint of every search interval, centre (0.25 m, 0.20 m), axes 0.027 m and 0.014 m, angle 0 and severity 0.675, matches none of the three patches. Against the middle patch the y error is 0.014 m, outside the 0.013 m centre cap. The other two patches miss that cap as well. A polished compromise-blob placement has forecast relative error 1.70e-3 and centres more than 0.06 m off. An intact closed branch, or a membrane left whole, leaves forecast relative errors 0.00284 and 0.00319. Parking one real patch at a corner misses the centre cap. The counts are in `authoring/evidence/zone-calibration.json`.

The artifact must be a regular file of at most 1 MiB, with unique keys. `plates` is a JSON array of one object. That object has three in-bound zones (`a>b`), two in-bound multipliers, numeric branch signs and a positive snap load. The forecast must be two rows at the requested forces with four positive frequencies each. Branch signs must match. Each zone is matched, in any order, to the authored patch inside the caps above. The snap relative error and each forecast-frequency relative error must be at most 0.0008. A 1e-12 allowance is included.

```sh
LDI_ARTIFACT=/tmp/damage.json pytest tests/test_outputs.py -q
```
