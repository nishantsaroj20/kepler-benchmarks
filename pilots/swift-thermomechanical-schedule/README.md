# Swift thermomechanical schedule pilot

Stage 1–3 only: diagnosis, specification, and an executed feasibility pilot. This directory is not a Harbor task and not a submission ZIP. `instruction.md` is intentionally absent.

## Reproduce

From this directory:

```bash
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 src/run_pilot.py
```

Or `scripts/run_pilot.sh`. Existing `results/*.json` checkpoints are reused. `PILOT_RERUN=1` recomputes them. The run needs Python 3.12, NumPy 2.1.3 and SciPy 1.14.1. No network, no Numba, no JAX.

The logged pilot (`logs/pilot.log`) took about 15 minutes of checkpointed solves on 4 CPUs. A later 28 s supplement, `results/thermal_warm_supplement.json`, reran FAST, H4, M2 and SLOW on the warm start in both adiabatic and coupled modes. That file is how the M2 feasibility flip (334.45 K coupled, 373.82 K adiabatic) was measured.

## Provenance

| Item | Value |
|---|---|
| Python | 3.12.3 |
| NumPy | 2.1.3 |
| SciPy | 1.14.1 |
| CPUs | 4, one thread |
| Calibration seed | 20261108 |
| Generating card | K=680, eps0=0.015, n=0.18, alpha=0.5, Q=220, beta=14, C=0.012, m=0.8 |
| Original v1 ZIP | SHA-256 `0c154a40d1f6d594dab71b35147d6abf3fb66a69589d6593c2d491d4b95fc43d` (not modified, not in this tree) |

SHA-256 of the files that define the executed study:

| File | SHA-256 |
|---|---|
| `src/coupled_current.py` | `e60110c8ee9bddbc59e7eaefff700e695e27c2f93345d8d62d9ddf8e7d411f18` |
| `src/reference_pk1.py` | `155b7b70607ecce9c1b841abe53b09403848ae3a3948a0f0f2d3b11b8c4752f6` |
| `src/material_checks.py` | `e51c6c42d4f538580823ad13d2200591e9db120f90c43e72e1ce058e1d0d8bb1` |
| `src/problem_spec.py` | `dd08456ea555f864b39965aef8a54345e4a7a118ce8e8fb4963649e8e505c6ad` |
| `src/run_pilot.py` | `fd0963a7737e56a3c019b2dd8ab8195347dbc98ab1f5f3bbddb49fb0ebf8fe29` |
| `results/calibration_public.csv` | `d43146d85f1789f1db43065ef4e9bcc5d30ee685298700f0f59b4e26b6baf12b` |
| `results/pilot_results.csv` | `af9861d976ee2a9f14b680cd5511bb1b9c2526595559448b2c2a41f02d2d2b4c` |
| `results/gate_verdicts.json` | `e3d8b78ce891f91be2ee774429ed3038730197ac23814ec8eedff80ca4acbd6f` |

## Read first

`correction.md` supersedes the PASS table in `individualreport.md`. The old report and `results/gate_verdicts.json` are kept. Corrected runs are in `results/corrected/`. The proportional family was also run on the pre-correction solver and saved at `results/precorrection_immutable/proportional_family.json`.

- `SPECIFICATION_AMENDMENT.md` — what the corrected runs change about the specification's status
- `specification.md` — equations and limits, with section 11 marked as amended
- `REQUIREMENTS_INVENTORY.md` — still not an instruction

