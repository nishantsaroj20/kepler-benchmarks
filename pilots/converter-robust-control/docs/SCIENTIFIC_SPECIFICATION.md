# Scientific specification

Plant set `crc-two-gfl-lcl-v1`. Controller schema `1.0.0`. Authoritative matrices for the executed run are in `runs/2026-11-08-stage0-4/plants/`. Source commit `595ff19f1b24955dfbf15423ab5f889342f49fa6`.

## 1. Circuit

Two converters, indexed 1 and 2, share a PCC. From each converter terminal:

- converter inductor `Lc`, `Rc`
- series damping resistor `Rd` into capacitor `Cf` whose internal voltage is `v_c`
- the capacitor terminal voltage is `v_f = v_c + Rd (i_c - i_o)`
- output branch `Lo = Lg_filter + L_line`, `Ro = Rg_filter + R_line`, current `i_o` into the PCC
- both output currents sum to the grid current `i_g` through `Lg`, `Rg` to the Thevenin voltage `v_g`

There is no PCC capacitor. Kirchhoff's current law `i_g = i_o1 + i_o2` is enforced by eliminating `v_pcc` and `di_o/dt` together. The 6-by-6 real system, with unknowns `(di_o1d, di_o1q, di_o2d, di_o2q, v_pccd, v_pccq)`, is

```
Lo1 di_o1d/dt + v_pccd = v_f1d - Ro1 i_o1d + ω0 Lo1 i_o1q
Lo1 di_o1q/dt + v_pccq = v_f1q - Ro1 i_o1q - ω0 Lo1 i_o1d
```

and the same pattern for converter 2, plus

```
Lg (di_o1d/dt + di_o2d/dt) - v_pccd = -v_gd - Rg i_gd + ω0 Lg i_gq
Lg (di_o1q/dt + di_o2q/dt) - v_pccq = -v_gq - Rg i_gq - ω0 Lg i_gd
```

The complex transcription uses the same laws as `Lo di_o/dt = v_f - v_pcc - Ro i_o - j ω0 Lo i_o` and solves a 3-by-3 complex matrix. It does not call the real stamp.

Inductor and capacitor equations in the grid frame, amplitude-invariant Park transform with `θ = ω0 t` and

```
v_d = (2/3) Σ v_k cos θ_k
v_q = (2/3) Σ -v_k sin θ_k
```

are

```
Lc di_cd/dt = v_td - v_fd - Rc i_cd + ω0 Lc i_cq
Lc di_cq/dt = v_tq - v_fq - Rc i_cq - ω0 Lc i_cd
Cf dv_cd/dt = (i_cd - i_od) + ω0 Cf v_cq
Cf dv_cq/dt = (i_cq - i_oq) - ω0 Cf v_cd
```

The cross-coupling signs follow by differentiating `i_a = i_d cos θ - i_q sin θ` and matching `v = Ri + L di/dt`. Positive `ω0 = 100π`.

Power, using filtered PLL-frame voltage and filtered grid-side current:

```
P = 1.5 (v_d i_od + v_q i_oq)
Q = 1.5 (v_d i_oq - v_q i_od)
```

`Q > 0` when current leads a d-axis voltage (`i_q > 0`, `v_q = 0`). That is var injection.

## 2. Frames, gauge, and fixed loops

The grid frame is fixed by placing the Thevenin voltage on the positive d-axis at every equilibrium: `v_gq = 0`, `v_gd = V_base`. There is no free common angle left. Each PLL angle `δ` is the angle of that inverter's frame ahead of the grid frame:

```
x_pll = R(δ) x_grid,  R = [[cos δ, sin δ], [-sin δ, cos δ]]
```

The PLL is a second-order loop on the filtered q-axis voltage:

```
d ξ_pll / dt = v_qf
d δ / dt = Kp_pll v_qf + Ki_pll ξ_pll
```

with `Kp_pll = 0.32643 rad/(V·s)` and `Ki_pll = 17.406 rad/(V·s²)`, from `ω_n = 2π·12`, damping 0.707, and division by `V_base`. At equilibrium `v_qf = 0` and `ξ_pll = 0`, so the frame sits at grid frequency.

The inner current PI is in the PLL frame and is not a decision variable:

```
d ξ_d / dt = i_d_ref - i_cdf
v_d_cmd = Kp_i (i_d_ref - i_cdf) + Ki_i ξ_d - ω0 Lc i_cqf + v_df
```

and the q-axis equation with the opposite decoupling sign `+ ω0 Lc i_cdf`. `Kp_i = 0.54325 V/A`, `Ki_i = 4.2667 V/(A·s)`, from `Kp = Lc / (3 τ_a)` and cancellation of the converter-inductor ESR pole only. Decoupling uses fixed `ω0`, not the instantaneous PLL speed. The command is rotated back to the grid frame and applied through the lag `τ_a = 1.5/8000 s`:

```
τ_a dv_t/dt = v_cmd_grid - v_t
```

Sensor filters are first-order at 800 Hz on `v_f`, `i_c`, and `i_o` after the rotation into the PLL frame.

`D_yu = 0` because `y` is a function of those filtered states and of the reference and noise channels only. The current reference affects `ξ` and `v_t` through derivatives and the actuator lag. There is no algebraic loop. The published `Dyu` matrices are exactly zero; they were not overwritten after the fact.

## 3. Parameters

Synthetic, not measured. Bases use `P = 1.5 V_peak I_peak`, so `V_base = 326.5986 V`, `I_base = 204.1241 A`, `Z_base = 1.6 Ω`.

| Quantity | Value | Rule |
|---|---|---|
| `Lc` | 0.30558 mH | 6% of `L_base`, inside the 1–10% Liserre band |
| `Rc` | 2.40 mΩ | `ω0 Lc / 40` |
| `Lg_filter` | 0.10695 mH | `0.35 Lc` |
| `Cf` | 79.577 µF | 4% of rated reactive power, under the 5% limit |
| `Rd` | 0.9978 Ω | `1/(ω_res Cf)` |
| `f_res` | 2004.5 Hz | inside `10·50 Hz` and `8000/2 Hz` |
| cable R | 0.3207 Ω/km | IEC 60228 Class 2 copper, 70 mm², 0.268 Ω/km at 20 °C, scaled to 70 °C with α = 0.00393 |
| cable L | 0.25 mH/km | synthetic; IEC 60228 does not specify inductance |

Grid impedance is `Z = V_ll² / (SCR · 200 kVA)` with the stated X/R, defined at the PCC and not including the inverter cables.

| Case | SCR | X/R | Cable lengths | P pu | Q pu |
|---|---|---|---|---|---|
| nominal | 10 | 4 | 150 m, 150 m | 0.55, 0.55 | 0.08, 0.08 |
| weak_grid | 3 | 1.5 | 150 m, 150 m | 0.55, 0.55 | 0.08, 0.08 |
| high_power | 10 | 4 | 150 m, 150 m | 0.85, 0.85 | -0.15, -0.15 |
| asymmetric | 6 | 3 | 80 m, 800 m | 0.70, 0.25 | 0.12, 0.00 |

## 4. Equilibria

Residuals are the infinity norm of the SI right-hand side after a nonlinear least-squares solve. Power imbalance is converter terminal power minus grid-source absorption minus ohmic loss, in watts.

| Case | RHS inf | Imbalance W | δ rad | v_t peak V | u bias pu (id, iq) |
|---|---|---|---|---|---|
| nominal | 3.94e-10 | -9.1e-12 | 0.07217, 0.07217 | 330.31, 330.31 | 0.5417, 0.1193 |
| weak_grid | 3.70e-8 | 4.4e-11 | 0.1737, 0.1737 | 348.42, 348.42 | 0.5138, 0.1175 |
| high_power | 6.24e-10 | -3.3e-11 | 0.09313, 0.09313 | 348.09, 348.09 | 0.8048, -0.0996 |
| asymmetric | 5.71e-10 | 3.6e-12 | 0.09544, 0.09078 | 330.33, 341.94 | (0.6879, 0.1585), (0.2394, 0.0419) |

Sine-triangle peak phase voltage is `750/2 = 375 V`. Every equilibrium is below that. Space-vector headroom is larger. This does not model saturation.

Symmetric cases produce symmetric states. Active-power export produces a small positive angle, as expected on an inductive feeder. Weak-grid loss is 11.9 kW versus 3.6 kW nominal, which is the price of the larger Thevenin resistance at this schedule.

## 5. State vector and scaling

Eighteen states per converter, converter 1 then converter 2:

`i_cd, i_cq, v_cd, v_cq, i_od, i_oq, ξ_id, ξ_iq, v_td, v_tq, ξ_pll, δ, v_df, v_qf, ic_df, ic_qf, io_df, io_qf`

SI units are A, V, A·s, V·s, and rad. The contract plant uses `x_si = T x_scaled` with a fixed diagonal `T`: currents by `I_base`, voltages by `V_base`, current integrators by `I_base/ω0`, PLL integrator by `V_base/ω0`, angle by 1. The same `T` is used for every case. Open-loop condition numbers drop from about 5e6 in SI to about 2.7e3 after scaling. Eigenvalues are unchanged by the similarity.

Eight weight states are appended: four error-weight states, then four effort-weight states. Physical states are indices 0–35. Weight states are 36–43. They are analysis states.

Contract dimensions: `n_x = 44`, `n_w = 10`, `n_u = 4`, `n_y = 4`, `n_z = 8`.

Channel order:

- `w`: `v_gd_5pct, v_gq_5pct, p_ref_1, q_ref_1, p_ref_2, q_ref_2, n_p_1, n_q_1, n_p_2, n_q_2`
- `u`: `i_d_ref_1_pu, i_q_ref_1_pu, i_d_ref_2_pu, i_q_ref_2_pu`
- `y`: `e_p_1, e_q_1, e_p_2, e_q_2`
- `z`: the four weighted errors, then the four weighted efforts

Matrix SHA-256 of the contract JSON matrices:

| Case | SHA-256 prefix |
|---|---|
| nominal | `99121d1490123b16` |
| weak_grid | `250fa1e099f01c79` |
| high_power | `cb005460282d1283` |
| asymmetric | `991cf551542fbd17` |

Full hashes are the `matrix_sha256` fields in the plant files.

## 6. Two construction paths

Path A builds the Jacobian from the constant branch matrix and the analytic derivatives of the rotations, PI, PLL, and sensors (`jacobian_analytic.py`). Path B is a central-difference Jacobian of an independently coded real stamp (`rhs_real.py`). A third residual, the complex stamp (`rhs_complex.py`), is compared to the real residual and its finite-difference Jacobian is compared to path A.

At the four equilibria the real-versus-complex residual is at most `3.7e-8` (weak grid) and is 0 on the nominal point to solver tolerance. Relative Frobenius error of path A versus real finite differences is about `1.7e-11` at step `1e-5`, and larger steps and smaller steps remain below `2e-9`. Singular values of `Cy (jωI - A)^{-1} Bu` from path A and from the complex-stamp finite-difference matrix differ by at most `1.9e-11`.

The open physical loop, outer controller open, is stable in every case. The spectral abscissa is about `-7.86` rad/s. The LCL resonance is present near 1.9 kHz and is well damped by `Rd` and the current loop; it was not deleted.

Scaled DC gain from `u` in per-unit current to `y` is approximately `-1` on the matching P and Q diagonals and a few hundredths off-diagonal. Positive controller gain is therefore negative feedback.

## 7. Closed loop

For `D_yu = 0` the contract is

```
Acl = [[A + Bu Dk Cy, Bu Ck], [Bk Cy, Ak]]
Bcl = [[Bw + Bu Dk Dyw], [Bk Dyw]]
Ccl = [Cz + Dzu Dk Cy, Dzu Ck]
Dcl = Dzw + Dzu Dk Dyw
```

with `xc = [xc1, xc2]`, `y` and `u` in the orders above, and `Ak, Bk, Ck, Dk` the block-diagonal stacks. A second routine writes the same blocks from the series interconnection without calling that formula. They match on the replayed controllers.

Internal stability uses the same formula on the 36 physical states only. H-infinity uses the 44-state weighted plant. Weight poles do not appear in `y`, so they are not moved by the controller.

## 8. Controller file

`/app/output/controller.json` in a future bundle. The research checker is `schema.py`.

- Regular file, not a symlink, size 1 to 65536 bytes.
- UTF-8 JSON. Duplicate keys are rejected.
- Top-level keys exactly `schema_version`, `plant_set_id`, `blocks`.
- `schema_version` is `1.0.0`. `plant_set_id` is `crc-two-gfl-lcl-v1`.
- `blocks` has two objects with ids `converter_1` and `converter_2`, each exactly the keys `id, A, B, C, D`.
- Each array is 2 by 2, finite numbers. Booleans, strings, NaN, and infinities are rejected.
- Column order of `B` and `D` is `(e_p, e_q)`. Row order of `C` and `D` is `(i_d_ref_pu, i_q_ref_pu)`.
- Any real 2-state realization with the same input-output behavior is allowed. `σ_max(D) ≤ 15` on each block.

## 9. Evaluation used in the pilot

1. Physical stability with the margin in the rationale.
2. Weighted loop clearly Hurwitz, then SLICOT `AB13DD` (`slycot 0.7.1`, continuous, identity `E`, equilibrated). This is an L-infinity routine. It is treated as H-infinity only after the Hurwitz gate. On `1/(s-1)` it returns 1, and that result is recorded as L-infinity, not as a pass.
3. An independent refined frequency lower bound (grid plus bounded searches on log-frequency). Agreement with SLICOT within `max(0.03, 1% of the value)` is required for a pass.
4. A Hamiltonian bisection is reported when its imaginary-axis test both fires below the peak and clears above it. On the small catalog (`1/(s+1)`, feedthrough 2.5, second-order peak 10.01252, cancellation 0.5) the bracket matches SLICOT. On a resonance of damping `1e-4` the coarse probe misses the peak, the test is marked inactive, and SLICOT still returns 5000.000025. Inactive is not relabeled as a norm.
5. Bounded-real lemma in Clarabel, accepted only if the returned `P` reproduces an LMI whose largest eigenvalue is negative. This certifies `||G||_∞ < γ` at a chosen `γ`. It certified `γ = 1` for the passing PI loops (LMI eigenvalue about `-0.07`). It did not certify a tight level on the narrow catalog resonance. It is a threshold certificate, not a tight norm solver.

MATLAB `hinfstruct` was not run. No MATLAB installation is available. That baseline is UNVERIFIED, not failed.

## 10. Omitted effects

Switching harmonics and dead time, beyond the averaged lag. DC-bus dynamics. Negative-sequence unbalance. Frequency disturbances. Saturation, current limits, and anti-windup. Communication. Skin effect and core saturation. Cable capacitance. Any guarantee on a continuum of SCR or power between the four cases.

## 11. Citations read for this specification

- Huang, Xin, Dörfler, arXiv:1906.11331, HTML text. General H-infinity and grid-following discussion. Numerical example not used.
- Liserre, Blaabjerg, Hansen, IEEE Trans. Ind. Appl., 2005, abstract, plus the open step-by-step restatement at `https://busarello.prof.ufsc.br/LCL%20Filter%20Design%20Procedure.pdf` for the 5% capacitor limit, the resonance window, and `Rd = 1/(ω_res Cf)`.
- IEC 60228 Class 2 copper resistance, 70 mm², 0.268 Ω/km at 20 °C, as reproduced by Addison Cables, globalgrid2050, and the Lapp conductor table. Inductance is not taken from those tables.
- SLICOT `AB13DD` documentation, Bruinsma and Steinbuch, Systems & Control Letters, 1990, as cited there and by MathWorks `hinfstruct`.
- Python Control `linfnorm` documentation, read 8 Nov 2026. Package not installed.
- MathWorks `hinfstruct` documentation, read 8 Nov 2026. Not executed.
- Park sign and power were derived and checked in this repository against the circuit residual and the power balance, not copied from a solved benchmark.
