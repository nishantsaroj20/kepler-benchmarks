# Cheat / defective-method attempts

Attempted against practice + held-out families:

1. **All-zero JSON** — fails `fn_total`, `p_peak`, `contact_length` on every pack.
2. **Constant pressure / unit reactions** — fails relative pressure and force checks.
3. **Frictionless normal-only IC solve** — can match `fn_total` / `p_peak` on some packs but fails `ft_total` and `mobilization` whenever $ |F_t| > 0 $.
4. **Isolated Hertz per patch** (ignore elastic coupling) — fails multi-patch held-out packs on `p_peak` / `contact_length` once clearances engage neighbors.
5. **Final-step Coulomb projection without slip history** — drifts on `shear_reverse` packs (`ft_total` sign/mobilization).
6. **Nop / missing `/app/solver/analyze.py`** — verifier presence test fails; reward 0.
