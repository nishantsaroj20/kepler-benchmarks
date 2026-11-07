# Breathing midplane-split identification

## Difficulty

The plate is a compressed eight-ply GFRP laminate with one midplane split that is open on the upward preparation and closed on the downward preparation. Membrane stiffness stays intact on both preparations. Bending stiffness switches between the intact laminate and the sum of the two four-ply sublaminate rigidities, blended by a separated-area fraction. That fraction is carried by two compact elliptical zones, so a trial zone that misses the real ellipse does not change the operators. Edge restraint on two sides is uncertain, four of seven modes are missing, and the preparation label was not stored with the record.

A shared ply-group knockdown, or one bending field used for every hold, cannot reproduce both preparations. Random local optimization from generic starts settles in decoy basins whose modal residual is several times the known noise; those basins miss the open-branch frequency forecast. The forecast forces are outside the held set, and the snap load is the first stability loss of the closed branch rather than an equal-energy or linear-buckling load.

The records are seeded synthetic output of the nine-term shallow-plate model documented with the data. They are not measurements from a physical panel. Geometric stiffening, the known imperfection, support uncertainty, incomplete shapes and the breathing interface are the modelling content. Sixteen continuous unknowns plus four preparation signs remain; the compact support is what makes a blind local search fail.

## Reference solution

The reference assembles the intact membrane energy with the existing in-plane condensation, then replaces the bending and spring operator by the open or closed split field. For each held force it compares both preparations, assigns four measured columns to distinct predicted modes, and keeps the lower residual. Centres are scanned on a plate grid with one zone active, because the elliptical support has no gradient until it overlaps the split. The two best separated centres seed a bounded least-squares polish over both zones, the two multipliers and a few orientation pairs. No authored parameter is used as a start.

The closed-branch snap is the last positive force at which a bisection still finds a stable negative equilibrium. The open-branch forecast evaluates the lowest four frequencies of the positive equilibrium at the two requested forces.

```sh
python3 authoring/provenance/generate_data.py
LDI_DATA_DIR=environment/data LDI_OUT_FILE=/tmp/damage.json python3 solution/solve_damage.py
```

On the seeded records this search finishes at a noise-level residual, with the authored preparation signs, a snap error of about 4.5e-5 relative, and a largest forecast error of about 4.7e-5 relative.

## Verification

The isolated checks read only the JSON artifact and the authored truth. They do not import a solver. Zone and fixture values are the generator inputs. The snap load and the eight forecast frequencies are forward predictions of that same input, recomputed by the generator before any inversion. The environment image contains the model statement and the four records, not the generator seed or the truth file.

The artifact must be a regular file of at most 1 MiB, with unique keys, one plate, two in-bound zones (`a>b`), two in-bound multipliers, numeric branch signs and a positive snap load. The forecast must be two rows at the requested forces with four positive frequencies each. Branch signs must match. The snap relative error must be at most 0.01. Each forecast frequency relative error must be at most 0.0005. A 1e-12 allowance is included.

Four random bounded starts on this realization, using the same residual and no centre scan, produced forecast errors of 9.9e-4, 1.5e-3, 1.8e-3 and 2.8e-3. All four miss 0.0005. The centre scan above meets it. The limit is several times the reference forecast error and below those decoy forecasts. It is not a claim that every noise-level description of the ellipses is unique in angle or severity; those details move more than the open-branch frequencies do, so the frequencies, the signs and the snap are the required physical results.

```sh
LDI_ARTIFACT=/tmp/damage.json pytest tests/test_outputs.py -q
```
