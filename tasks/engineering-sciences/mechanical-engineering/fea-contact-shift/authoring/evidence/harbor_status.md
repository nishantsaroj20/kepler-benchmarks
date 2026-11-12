# Harbor / Docker status

This authoring VM may lack Docker/`harbor`. Local stand-in:

```bash
python3 authoring/evidence/run_local_validate.py
```

Quality-review remediations applied (2026-11-12):

- Approach-controlled $(\delta_n,\delta_t)$ so $F_t$ cannot be copied from the load file
- Graded frictional-field metrics (`q_peak`, stick/slip, `cum_slip_l1`)
- Independent Gauss-quadrature truth (not oracle import / not force-file echo)
- Unprivileged runner harness, sealed `/tests`, `O_NOFOLLOW` outputs
- Practice packs + no-contact cases exercised; exact JSON key set; integer `schema_version`

When Docker+Harbor are available:

```bash
harbor run -p . -a oracle   # expect reward 1
harbor run -p . -a nop      # expect reward 0
```
