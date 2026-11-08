# Requirements rationale

All numerical thresholds below were fixed in `src/crc/params.py` and copied into `FROZEN_REQUIREMENTS.json` at the start of `runs/2026-11-08-stage0-4`, before the baseline search. They were not edited after the PI grid succeeded.

## Channels

The controller input `y` is dimensionless:

- `e_p = (P_schedule - P_filtered) / 100 kVA`
- `e_q = (Q_schedule - Q_filtered) / 100 kVA`

Positive `e_p` means the inverter is short of its scheduled active power. Positive reactive power is the leading-current convention: the inverter injects vars when it behaves like a capacitor. The schedule used in the linear model is the equilibrium schedule plus the reference disturbance.

The controller output `u` is a perturbation of the current reference in per-unit of `I_base = 204.124 A` peak. The total reference seen by the inner loop is the published bias plus `u * I_base`.

Disturbance channels, each with unit infinity-norm weight in the stacked problem:

| Channel | Physical amplitude when the channel equals 1 |
|---|---|
| `v_gd_5pct`, `v_gq_5pct` | 5% of the 326.599 V peak phase base, on the grid Thevenin voltage |
| `p_ref_*`, `q_ref_*` | 10% of 100 kVA on that power reference |
| `n_p_*`, `n_q_*` | 0.5% of 100 kVA added to the measured power before the error is formed |

These are design amplitudes for a feeder voltage disturbance, a dispatch change, and ordinary sensor noise. They are not recorded events.

## One stacked norm

The requirement is `||T_zw||_∞ ≤ 1` on the weighted generalized plant, the same controller, every case. Submatrix norms are smaller than or equal to the stacked norm, so separate channel limits are not additional pass tests. They are reported only as diagnostics. A frequency-grid peak is not this norm.

## Weights

Error weight, each of the four power errors:

`W_e(s) = (1/3) (s + 37.699) / (s + 0.12566)`

This is `(s/M + ω_b) / (s + ω_b ε)` with `ω_b = 2π·2` rad/s, `M = 3`, `ε = 0.01`. The 2 Hz corner is an outer-loop bandwidth below the 12 Hz PLL and far below the current loop. DC gain is 100, so a sustained 0.01 pu power error already produces a weighted output of 1. High-frequency gain is 1/3. The weight asks for integral action without demanding that the outer loop outrun the PLL.

Effort weight, each current-reference channel:

`W_u(s) = 2 (s + 188.50) / (s + 1885.0)`

DC gain is 0.2, so a sustained 5 pu current command would be required to saturate the weight. That is intentionally loose at DC: delivering the scheduled power needs a current of order 1 pu. High-frequency gain is 2, so a 0.5 pu high-frequency current command saturates the weight. Together with the 0.5% noise channel, this limits chatter without forbidding a normal proportional gain.

## Stability

On the scaled physical state plus the controller state, every eigenvalue must satisfy `Re(λ) + δ ≤ -0.5` rad/s, where `δ` is a Bauer-Fike style budget `κ_eig · 1e-8 · ||A||_F`. If the uncertainty interval straddles -0.5, the case is unresolved and does not pass. The 0.5 rad/s floor is a clear left-half-plane margin, about a 2 s dominant time constant, not a bandwidth specification. Bandwidth sits in `W_e`.

Analysis-weight poles are excluded from that margin. `W_e` has a pole at -0.126 rad/s by construction. Counting it would make the specification fail for every stabilizing controller. Those poles do not feed back into the equipment. The weighted loop used for the H-infinity norm must still be clearly Hurwitz (`Re + δ ≤ -0.001`), because a finite L-infinity number is not an H-infinity number on an unstable system. Python Control's `linfnorm` documentation states that distinction. The package itself was not imported: PyPI serves no files for `python-control` while the project is quarantined.

## Gain bound

`σ_max(D)` of each 2-by-2 block is at most 15. `D` is the high-frequency gain and is realization-invariant. Fifteen per-unit of current per per-unit of power error is far above a normal outer-loop proportional gain and is an implementation ceiling, not an oracle match. No entrywise bound is placed on `A`, `B`, or `C`, because a change of state basis can make those entries arbitrary while leaving the input-output behavior unchanged.

## What was not required

No continuous-parameter μ certificate. No anti-windup. No PWM spectrum. No requirement that `Ak` itself be stable. No match to a hidden coefficient vector. The witness levels `γ ≤ 0.95` and `Re + δ ≤ -1` are reporting margins for a comfortable certificate. The pass levels remain 1 and -0.5.
