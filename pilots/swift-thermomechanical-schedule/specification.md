# Scientific specification: thermomechanical schedule design

Status: Stage-2 specification exercised by the Stage-3 pilot in this directory. Numerical limits below are the candidates the pilot actually ran. A limit is **validated** only when `results/gate_verdicts.json` and `results/pilot_results.csv` support it. Anything not backed by those files remains **proposed**. This document is not `instruction.md`. The human author writes that file.

## 1. Decision the task would ask for

Identify the eight Swift–Voce rate/temperature parameters from noisy public force, gauge, temperature and radius records, then choose a loading schedule

`d = (t1, t2, t3, h)`

that meets every published constraint on every public scenario. Any feasible schedule passes. There is no oracle schedule and no global-optimum requirement.

The v1 task (fixed notch forecasts after an inverse fit) was solved 8/8 by the frontier configuration that was shown to the author. This redesign changes the decision from “predict two fixed tests” to “choose durations under coupled constraints.” It does not claim a solve rate.

## 2. Units and conjugacy

| Quantity | Unit |
|---|---|
| Length, displacement | mm |
| Time | s |
| Force | N (full 360° axial reaction) |
| Stress | MPa = N/mm² |
| Temperature | K |
| Volumetric heat capacity `rho_c` | N/(mm²·K), per **reference** volume |
| Conductivity `k` | N/(s·K) |
| Diffusivity `k/rho_c` | mm²/s |
| Plastic work | N·mm per reference volume, because `q * d(ep)` with `q` in N/mm² |

`1 N/mm² = 10⁻³ J/mm³`. With `rho_c = 3.6 N/(mm²·K)` and `k = 45 N/(s·K)`,

`alpha = k / rho_c = 12.5 mm²/s = 1.25×10⁻⁵ m²/s`,

which is the usual mild-steel diffusivity (`k_SI ≈ 45 W/(m·K)`, `rho c ≈ 3.6 MJ/(m³·K)`).

**Kirchhoff versus Cauchy.** The yield stress and the plastic heating use the Kirchhoff von Mises stress

`q = sqrt(3/2) ||dev τ||`.

The dissipation per reference volume is `q * d(ep)`. It is not `σ : d` integrated on the current volume, and it is not multiplied by an extra Jacobian. The heating statement in the v1 specification (`rho_c dT = chi * sigma_flow * d(ep)` per reference volume) is the same statement, because `q = sigma_flow` on the yield surface.

Virtual work in the primary solver is the current-configuration form with Cauchy stress `σ = τ / J_bar` and the true current volume. The independent solver uses the equivalent reference form `P = (J / J_bar) τ F^{-T}`. Both include the F-bar factor `J/J_bar`.

## 3. Constitutive model (retained)

Elasticity is known and Hencky:

`ε_e = (1/2) ln b_e`, `τ = κ tr(ε_e) I + 2 μ dev(ε_e)`,

`μ = E / (2(1+ν))`, `κ = E / (3(1-2ν))`, `E = 210000 MPa`, `ν = 0.3`.

Multiplicative split `F = F_e F_p` with `det F_p = 1`. The plastic corrector is an exponential-map return on the principal logarithmic strains. The discrete consistency condition is

`q_trial - 3 μ Δγ - σ_flow(ep + Δγ, Δγ/Δt, T) = 0`.

Hardening, with `σ_y = K ε0^n`:

`σ_h(ep) = α K (ε0+ep)^n + (1-α) [σ_y + Q (1 - exp(-β ep))]`.

Flow stress, consistency viscoplasticity, no overstress, no creep below `edot0`:

`σ_flow = σ_h(ep) * (1 + C ln(max(ep_dot/edot0, 1))) * (1 - θ^m)`,

`θ = clip((T - T_ref)/(T_m - T_ref), 0, 1)`,

`edot0 = 0.1 /s`, `T_ref = 293 K`, `T_m = 1800 K`.

Unknowns, inclusive bounds (same numerical box as v1):

| Symbol | Bounds |
|---|---|
| K | 300–1500 MPa |
| ε0 | 0.001–0.1 |
| n | 0.05–0.45 |
| α | 0–1 |
| Q | 0–800 MPa |
| β | 0.5–80 |
| C | 0–0.1 |
| m | 0.3–3 |

Individual constants are not required to be unique. Cards are judged by the public observations and by whether they agree on the schedule constraints.

Synthetic generating card used only inside this pilot (not a hidden verifier key, because there is no verifier in this stage):

`K=680, ε0=0.015, n=0.18, α=0.5, Q=220, β=14, C=0.012, m=0.8`, so `σ_y = 319.305 MPa`.

Known thermal constants: `chi = 0.9`, `rho_c = 3.6`, `k = 45`.

### Explicit omissions

- No thermal expansion (`F_θ` is absent; `F = F_e F_p` only).
- No Gough–Joule / reversible thermoelastic heat. Only the plastic dissipation fraction `chi` heats the body.
- No damage, fracture, contact, anisotropy, kinematic hardening, or inertia.

## 4. Heat equation

Reference configuration, backward Euler, lumped capacity:

`∫ rho_c (T - T_n)/Δt  N_a dV0 + ∫ Grad N_a · (k J C^{-1} Grad T) dV0 = ∫ N_a chi q Δγ/Δt dV0`.

`C = F^T F` uses the **true** deformation gradient. The hoop temperature gradient is zero, so only the in-plane block of `C^{-1}` is stored. `J = det F_true`. The source uses reference volume and Kirchhoff `q`.

The mechanical step freezes the Gauss-point temperature. The heat step then updates the nodal temperature. The pair is repeated until the Gauss-point change is under 0.25 K or two passes have been done. That split is **proposed** as adequate; the pilot’s time-refinement block is what validates it.

Adiabatic mode (`k = 0` and a local fixed-point inside the return map) is a different update: `T = T_n + chi * σ_flow * Δγ / rho_c`, solved with the yield stress. It is the continuum adiabatic limit, not a renaming of the conduction matrix.

## 5. Geometry, boundary conditions, scenarios

Half specimen, midplane symmetry. One geometry family.

| Item | Value | Status |
|---|---|---|
| R0 | 3.0 mm | used in the pilot |
| L_half | 10.0 mm | used |
| Nominal imperfection dR | 0.015 mm | used |
| Profile | `r(Z) = R0 - dR cos²(π Z / (2 L))` | used |
| Gauge material point | outer surface, reference Z = 4 mm | used |
| Gauge opening | `2 u_z` of that material point | used |
| Neck monitor | reference-volume average, `0 ≤ Z ≤ 1.5 mm` | used |
| Shoulder monitor | reference-volume average, `6 ≤ Z ≤ 8.5 mm` | reported, not capped |
| End displacement | sum of the three fixed increments | used |

Axial Fourier time `L²/α = 8.0 s`. Radial time `R²/α = 0.72 s`. Holds of a few seconds can move heat to the grips; a 1 s ramp cannot.

Boundary conditions, all scenarios:

- `Z = 0`: `u_z = 0`, heat flux 0.
- `R = 0`: `u_r = 0`, heat flux 0.
- Lateral surface: traction free, heat flux 0.
- `Z = L`: `u_z = u_end(t)`, `u_r` free, `T = T_grip = 293 K`.

Public scenarios (no hidden condition):

| ID | T_init | dR |
|---|---|---|
| S0_nominal | 293 K | 0.015 mm |
| S1_warm | 333 K | 0.015 mm |
| S2_imperfection | 293 K | 0.030 mm |

Every scenario is required. A nominal-only schedule is not sufficient: the warm start raises the neck temperature by tens of kelvin.

## 6. Why the model is quasi-static

Actuator caps used for the duration floors: end speed `8 mm/s`, end acceleration `40 mm/s²`.

Elastic wave transit `L/c ≈ 0.01/5200 ≈ 2 μs`. The shortest ramp is about 0.3 s, so `10⁵` transits.

Inertia at the acceleration cap: `ρ L a ≈ 7800 kg/m³ * 0.01 m * 0.04 m/s² ≈ 3 Pa`, against a 300 MPa flow stress. `ρ v²` at 8 mm/s is under 1 Pa. The task is rate-dependent quasi-static thermomechanics.

## 7. Schedule family

Three fixed end-displacement increments of the half model:

`Δ1 = 0.70 mm`, `Δ2 = 0.75 mm`, `Δ3 = 0.55 mm` (total 2.00 mm).

These stop the bar past the load peak (end force about 1 kN under the peak on the nominal mesh) and before a deep neck whose radius is mesh-divergent. That choice is a **pilot result**, not an arbitrary extra unknown. Deeper extensions (total 2.6 mm, `ep ~ 1`) were run and rejected because the minimum radius still moved by several hundredths of a millimetre under one mesh level.

Inside a ramp of duration `t` and increment `δ`, the shape is the C2 smoothstep

`s(ξ) = ξ - sin(2πξ)/(2π)`, `ξ = (t_local)/t ∈ [0,1]`,

`u = u_start + δ s(ξ)`.

Then `v_max = 2 δ / t` and `a_max = 2π δ / t²`, attained inside the ramp, with `v = a = 0` at both ends. Concatenating ramps adds no extra transition and no hidden smoothing.

The optional hold `h ≥ 0` is at the stage-2 / stage-3 boundary, at fixed end displacement `Δ1+Δ2`. **Zero hold** means stage 3 starts at the instant stage 2 ends.

`τ = t1 + t2 + t3 + h`.

Inclusive box, lower edges equal to the actuator floors rounded up to 0.01 s:

| Variable | Inclusive bounds |
|---|---|
| t1 | 0.34–12 s |
| t2 | 0.35–12 s |
| t3 | 0.30–12 s |
| h | 0–16 s |

A schedule outside the box or with `v_max > 8` or `a_max > 40` is illegal. The fastest proportional schedule (durations proportional to the increments, scaled up until every floor is met) is legal and is one of the baselines.

## 8. Calibration observations

Three public tests on the nominal geometry, `T_init = T_grip = 293 K`, coupled conduction. The displacement history is the same smoothstep family.

| Test | (t1, t2, t3, h) s | Role |
|---|---|---|
| CAL_SLOW | (6, 6, 4, 0) | hardening, near-isothermal because `τ = 16 s` exceeds the axial diffusion time |
| CAL_FAST | (0.40, 0.42, 0.36, 0) | rate versus near-adiabatic heating |
| CAL_HOLD | (0.40, 0.42, 0.36, 6) | temperature decay at fixed displacement, then reload |

Sampling operator: linear interpolation of the accepted finite-element history onto the published times. Quantities:

- `t_s`, prescribed `u_end_mm`
- full-circumference axial force
- gauge opening
- neck-band and shoulder-band reference-volume average temperatures
- minimum current outer radius

The row at `t = 0` is exact: force 0, gauge 0, both temperatures 293 K, radius 2.985 mm. Every later row carries independent Gaussian noise with standard deviations 20 N, 0.002 mm, 0.4 K, 0.4 K and 0.003 mm. Seed `20261108`.

Generating mesh for the synthetic records: 5 radial × 12 axial, 8 substeps per ramp. The fit mesh is coarser (3×6, 4 substeps) so the identification is not an inverse crime on the generating grid.

## 9. Constraints

Evaluated on the solver’s accepted steps, including steps that fall inside a ramp and inside the hold, not on a sparse report. Cycle-wide extrema:

| Constraint | Candidate | Intended role |
|---|---|---|
| Neck-band average temperature | `≤ 345 K` at every accepted step | active |
| Cycle time | `τ ≤ 8 s` | active |
| Minimum outer radius | `≥ 2.20 mm` at every accepted step | safety floor; expected not to bind at this strain |
| Axial force | `≤ 14000 N` | safety floor; rate sensitivity is too weak to bind |
| Final gauge opening | in `[2.40, 2.95] mm` | active on the fast, lightly necked schedules |

Comfortable margin used when calling a schedule comfortably feasible: temperature at least 5 K under the cap, `τ` at least 1 s under the cap, radius at least 0.05 mm over the floor, force at least 500 N under the cap, gauge at least 0.05 mm inside each side of the window.

These numbers were placed in the gaps measured on a 4×8 mesh before the confirmatory refinement (fast warm neck ~368 K, simple 4 s hold on the warm start stuck at ~352 K, moderate schedule M2 near ~334 K, slow schedule `τ = 14.5 s`). They are not “the optimum minus a tiny epsilon.”

## 11. What the executed pilot validated

Authoritative tables: `results/pilot_results.csv` (SHA-256 `af9861d976ee2a9f14b680cd5511bb1b9c2526595559448b2c2a41f02d2d2b4c`) and `results/gate_verdicts.json`.

Validated on the primary solver, separate mesh refinement (4×8 to 5×12) and separate time refinement (8 to 16 substeps):

- M2 = `(1.5, 1.6, 1.3, 1.0)` s passes every scenario. Neck temperature changes by at most 0.16 K, force by at most 0.1%, and the pass/fail label does not flip.
- FAST, the fastest proportional schedule, and a 4 s hold at the actuator floor do not pass. On the warm start their neck temperatures are 368 K, 368 K and 352.5 K.
- SLOW = `(5, 5.5, 4, 0)` s has neck temperature 306–331 K and fails only because `τ = 14.5 s`.
- Ten one-at-a-time ±15% perturbations of M2, plus ±10% on all four durations, all pass every scenario.
- On S0, conduction lowers the neck temperature by 10.7 K (4 s hold), 21.1 K (M2) and 33.2 K (SLOW) relative to the adiabatic update. On S1, M2 is 334.45 K with conduction and 373.82 K adiabatic, so the same schedule passes the 345 K cap only when conduction is active (`results/thermal_warm_supplement.json`).
- Radius and force did not bind. The 2.20 mm and 14 kN floors stayed inactive across the tested family. They remain safety limits, not design drivers.
- A deep-neck trial at 2.6 mm end displacement was rejected before these limits were frozen: the minimum radius was still moving by several hundredths of a millimetre per mesh level.

Not validated: a fresh frontier solve rate, a verifier, and an `instruction.md`. The independent adiabatic cross-check on a coarse 3×6 mesh matches force to 0.06% and leaves a 4.4–6.1 K neck-temperature gap. That gap is larger than the primary solver’s own refinement scatter and is reported in the gate file rather than smoothed away.

Final gauge is not the prescribed end displacement. Uniform stretch would give `2 * (4/10) * 2.00 = 1.60 mm`. Computed gauges in the pilot sit near 2.3–2.7 mm because the neck is inside the gauge.

## 10. What this stage does not contain

No `instruction.md`, no `task.toml`, no verifier, no submission ZIP, no Harbor run, no frontier probe. Frontier difficulty is **UNVERIFIED**. Numba, JAX and surrogates are not banned for a future agent; this pilot simply did not need them. A one-dimensional surrogate was not used as evidence of feasibility.
