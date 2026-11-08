# Two-converter robust-control pilot

Research code for Stages 0–4 of an original grid-following converter problem.
This directory is not a Harbor bundle, does not contain `instruction.md`, and
does not contain a sealed verifier.

Fresh entry, from this directory, with a new output path:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python scripts/run_fresh.py --output runs/<new-run-id>
```

The script refuses an existing output directory. Dependencies are pinned in
`requirements.txt` and were installed in a virtual environment; `python-control`
could not be installed because PyPI lists that project as quarantined.
