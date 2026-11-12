# Oracle reward-0 root cause (2026-11-12)

## Symptom

AfterQuery Reference verification: oracle scored `[0,0,0]` in ~2 s; nop scored 0 (expected).

## Comparison with `swift-notch-plasticity`

Swift grades a static `/app/output/...` artifact and does **not** setuid-exec the agent during pytest. Our task must re-run `/app/solver/analyze.py`, so isolation is required, but the first harness put the output JSON inside pytest’s root-owned `TemporaryDirectory`. The unprivileged `runner` child could not traverse that directory → nonzero exit → reward 0.

## Fix

`tests/harness.py` now writes `out.json` inside a runner-owned `/tmp/run_*` workdir, then the root verifier copies it to the pytest path. Also raised `RLIMIT_AS` to 8 GiB so NumPy/OpenBLAS VAS reservation does not kill the child, and aligned `test.sh` reward writing with the swift atomic write pattern.

Local confirmation (sudo harness + Gauss truth compare) passes on contacting and no-contact practice packs.
