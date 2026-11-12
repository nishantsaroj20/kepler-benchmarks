# Harbor / Docker status

This environment has no Docker daemon and no `harbor` CLI.

Local stand-in (2026-11-12):

```bash
python3 authoring/evidence/run_local_validate.py
```

Result summary in `local_validate.log`:

- Oracle analyzer matches independent upsampled truth on 4 practice + 8 held-out packs → treated as reward **1**
- Missing program / zero / constant baselines fail → treated as nop / defective reward **0**
- Max practice oracle wall ~0.095 s; instruction per-run cap 8 s (~80× oracle, above the 2–3× guideline for CI margin on small packs)

When Docker+Harbor are available on a submit machine:

```bash
harbor run -p . -a oracle   # expect reward 1
harbor run -p . -a nop      # expect reward 0
```
