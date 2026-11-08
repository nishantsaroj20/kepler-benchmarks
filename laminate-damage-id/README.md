# GFRP plate with a breathing midplane split

## Difficulty

The plate is a compressed eight-ply GFRP laminate with one midplane split that is open on the upward preparation and closed on the downward preparation. Membrane stiffness stays intact on both preparations. Bending stiffness switches between the intact laminate and the sum of the two four-ply sublaminate rigidities, blended by a separated-area fraction. That fraction is carried by two compact elliptical zones, so a trial zone that misses the real ellipse does not change the operators. Edge restraint on two sides is uncertain, four of seven modes are missing, and the preparation label was not stored with the record.

A shared ply-group knockdown, or one bending field used for every hold, cannot reproduce both preparations. Random local optimization from generic starts settles in decoy basins whose modal residual is several times the known noise; those basins miss the open-branch frequency forecast. The forecast forces are outside the held set, and the snap load is the first stability loss of the closed branch rather than an equal-energy or linear-buckling load.

The records are seeded synthetic output of the nine-term shallow-plate model documented with the data. They are not measurements from a physical panel. Geometric stiffening, the known imperfection, support uncertainty, incomplete shapes and the breathing interface are the modelling content. Sixteen continuous unknowns plus four preparation signs remain; the compact support is what makes a blind local search fail. A composite structural-health-monitoring engineer, or a structural-dynamics researcher, performs this work in practice when a compressed laminate panel must be assessed from incomplete vibration records.

## Reference solution

The reference assembles the intact membrane energy with the existing in-plane condensation, then replaces the bending and spring operator by the open or closed split field. For each held force it compares both preparations, assigns four measured columns to distinct predicted modes, and keeps the lower residual. Centres are scanned on a plate grid with one zone active, because the elliptical support has no gradient until it overlaps the split. The two best separated centres seed a bounded least-squares polish over both zones, the two multipliers and a few orientation pairs. No authored parameter is used as a start.

The closed-branch snap is the last positive force at which a bisection still finds a stable negative equilibrium. The open-branch forecast evaluates the lowest four frequencies of the positive equilibrium at the two requested forces.

```sh
python3 authoring/provenance/generate_data.py
LDI_DATA_DIR=environment/data LDI_OUT_FILE=/tmp/damage.json python3 solution/solve_damage.py
```

On the seeded records this search finishes at a noise-level residual, with the authored preparation signs, a snap error of about 4.5e-5 relative, and a largest forecast error of about 4.7e-5 relative.

## Verification

The isolated checks read only the JSON artifact and `tests/ground_truth.json`. They do not import a solver. Zone centres, ellipse axes, severities, fixture multipliers and preparation signs in that file are the authored generator inputs. The snap load and the eight forecast frequencies are not taken from the reference condensation. `authoring/provenance/independent.py` integrates the membrane strains on the published 30 by 26 rule, eliminates the in-plane coefficients by a dense factorisation, and forms the tangent by central differences of that residual. A Newton continuation then brackets the closed-branch eigenvalue crossing. `generate_data.py` writes those independent values into the truth file.

On the authored parameters that direct model and the reference polynomial condensation differ by at most 6.4e-10 relative in the compared frequencies and by 5.4e-7 relative in the snap load. The stored reference inversion, checked against the direct model, is 4.7e-5 relative on the forecast and 4.6e-5 relative on the snap. The 0.0005 limits are about ten times those inversion errors and about a thousand times the direct-versus-reference difference. Four random bounded starts that share the reference residual, and do not scan centres, miss the forecast limit; their errors are 9.9e-4, 1.5e-3, 1.8e-3 and 2.8e-3.

A 40 by 36 quadrature in the same direct code moves the open-branch frequencies at -90 N by 1.17e-3 relative. That is a different rule from the published 30 by 26 model, so it is outside the 0.0005 contract. The environment image contains the model statement and the four records, not the generator seed or the truth file.

The artifact must be a regular file of at most 1 MiB, with unique keys, one plate, two in-bound zones (`a>b`), two in-bound multipliers, numeric branch signs and a positive snap load. The forecast must be two rows at the requested forces with four positive frequencies each. Branch signs must match. Each zone is matched, in either order, to the authored patch: centres within 0.020 m, semi-axis `a` within 0.010 m, semi-axis `b` within 0.008 m, orientation within 30 degrees modulo 180, and severity within 0.12. Each spring multiplier must be within 0.05 of the authored value. The snap relative error and each forecast-frequency relative error must be at most 0.0005. A 1e-12 allowance is included. The reference search above sits inside every one of those caps.

```sh
LDI_ARTIFACT=/tmp/damage.json pytest tests/test_outputs.py -q
```
