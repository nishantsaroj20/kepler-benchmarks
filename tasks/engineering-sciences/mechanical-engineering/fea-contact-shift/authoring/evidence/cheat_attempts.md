# Cheat / defective-method attempts

1. **All-zero JSON** — fails force and frictional-field checks on contacting packs.
2. **Constant / unit outputs** — fails relative pressure and force checks.
3. **Frictionless / $\delta_t=0$ during solve** — matches normals on some packs but fails `ft_total`, `q_peak`, and `cum_slip_l1` whenever the true history carries shear.
4. **Isolated Hertz per patch** — fails multi-patch packs on `p_peak` / `contact_length` once clearances engage neighbors.
5. **Import `/tests/truth_solver.py`** — blocked: analyzer runs as unprivileged `runner` with `/tests` mode `go-rwx`.
6. **Symlink output to sealed truth** — blocked: harness opens outputs with `O_NOFOLLOW` and refuses symlinks.
7. **Nop / missing `/app/solver/analyze.py`** — presence test fails; reward 0.
