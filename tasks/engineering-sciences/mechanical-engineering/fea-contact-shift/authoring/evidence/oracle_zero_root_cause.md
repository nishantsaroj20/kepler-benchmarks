# Oracle reward-0 root cause

## Symptom

Portal reference verification: oracle scored `[0,0,0]` in ~2–3 s; nop scored 0 (expected).

## Root causes (two layers)

1. **Output path (earlier):** setuid `runner` wrote `--out` into a root-owned pytest `TemporaryDirectory` and could not create the file.
2. **Pack path (current failure):** held-out packs are materialised under pytest temp dirs with mode `0700` (root-only). Chmodding the leaf pack to `0555` is not enough — `runner` still cannot traverse the parent, so `Path.is_dir()` raises `PermissionError` and the analyzer exits 1. Practice packs under `/app/data/practice` (world-traversable) passed, which hid the bug in partial local checks.

## Fix

`tests/harness.py` copies each pack into a runner-owned `/tmp/run_*/pack` tree and writes `out.json` there, then the root verifier copies the JSON back. Full local suite under sealed `/tests`: 5 passed (~3 s).

## Contrast

`swift-notch-plasticity` grades a static `/app/output` artifact and never setuid-execs during pytest.
