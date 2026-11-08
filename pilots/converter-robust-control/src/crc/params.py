"""Original synthetic parameters for two grid-following LCL converters.

Every derived number is computed from the design rules in this file.
These are not laboratory measurements and they are not the numerical
example of Huang, Xin and Dörfler, arXiv:1906.11331 (50 kVA, 380 V,
0.05 pu LCL, static 3-by-7 gain).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Ratings: synthetic commercial LV inverter class. Not a specific product.
S_RATED = 100_000.0  # VA, three-phase, per converter
V_LL_RMS = 400.0  # V, line-to-line rms
F0 = 50.0  # Hz
W0 = 2.0 * math.pi * F0
VDC = 750.0  # V, stiff DC bus used only for modulation-headroom audit
FSW = 8_000.0  # Hz, carrier used only to set the averaged PWM lag
TAU_A = 1.5 / FSW  # s, 1.5 carrier periods, first-order lag

# Amplitude-invariant dq base. P = 1.5*(vd*id + vq*iq).
V_BASE = V_LL_RMS * math.sqrt(2.0 / 3.0)  # peak phase voltage
I_BASE = S_RATED / (1.5 * V_BASE)
Z_BASE = V_BASE / I_BASE  # equals V_LL_RMS**2 / S_RATED
L_BASE = Z_BASE / W0
C_BASE = 1.0 / (W0 * Z_BASE)

# LCL rules, disclosed:
# converter reactance inside the 1–10% band used by Liserre-style design notes;
# capacitor reactive power below the 5% Liserre limit;
# grid-side filter inductor as a fixed fraction of Lc (ripple split);
# series resistance from a stated inductor quality factor at f0;
# series Rd equal to the capacitor reactance at the LCL resonance
# (design-note form of the Liserre passive-damping starting value).
LC_REACTANCE_PU = 0.06
CF_REACTIVE_PU = 0.04
LG_OVER_LC = 0.35
INDUCTOR_Q = 40.0

LC = LC_REACTANCE_PU * L_BASE
X_LC = W0 * LC
RC = X_LC / INDUCTOR_Q
LG_FILTER = LG_OVER_LC * LC
X_LG = W0 * LG_FILTER
RG_FILTER = X_LG / INDUCTOR_Q
Q_CAP = CF_REACTIVE_PU * S_RATED
CF = Q_CAP / (V_LL_RMS**2 * W0)
W_RES = math.sqrt((LC + LG_FILTER) / (LC * LG_FILTER * CF))
F_RES = W_RES / (2.0 * math.pi)
RD = 1.0 / (W_RES * CF)

# Cable class used for the inverter-to-PCC laterals.
# Resistance: IEC 60228 Class 2 copper, 70 mm^2, maximum DC resistance
# 0.268 ohm/km at 20 C, corrected to 70 C with alpha = 0.00393 /C.
# Inductance is NOT specified by IEC 60228. 0.25 mH/km is a synthetic
# value (0.0785 ohm/km at 50 Hz) placed inside the range of published
# 70 mm^2 cable reactances; it is not a copied catalog row.
COPPER_ALPHA = 0.00393
CABLE_R_PER_KM = 0.268 * (1.0 + COPPER_ALPHA * (70.0 - 20.0))
CABLE_L_PER_KM = 0.25e-3

# Factory-fixed loops. Not decision variables of the outer controller.
CURRENT_DELAY_FACTOR = 3.0  # kp = Lc / (factor * tau_a) => crossover ~ 1/(factor*tau)
KP_I = LC / (CURRENT_DELAY_FACTOR * TAU_A)
KI_I = KP_I * (RC / LC)  # cancels the converter-inductor ESR pole only
PLL_HZ = 12.0
PLL_XI = 0.707
PLL_WN = 2.0 * math.pi * PLL_HZ
KP_PLL = 2.0 * PLL_XI * PLL_WN / V_BASE
KI_PLL = PLL_WN**2 / V_BASE
SENSOR_HZ = 800.0
TAU_SENSOR = 1.0 / (2.0 * math.pi * SENSOR_HZ)

# Frozen performance normalizations and weights. Do not retune after baselines.
A_V_PU = 0.05  # 1.0 on a grid-voltage disturbance channel is a 5% voltage
A_REF_PU = 0.10  # 1.0 on a reference channel is a 10% power reference
A_NOISE_PU = 0.005  # 1.0 on a noise channel is 0.5% rated-power measurement noise
W_E_HZ = 2.0
W_E_M = 3.0
W_E_EPS = 0.01
W_E_WB = 2.0 * math.pi * W_E_HZ
W_E_KINF = 1.0 / W_E_M
W_E_ZERO = W_E_M * W_E_WB
W_E_POLE = W_E_WB * W_E_EPS
W_U_KINF = 2.0
W_U_ZERO = 2.0 * math.pi * 30.0
W_U_POLE = 2.0 * math.pi * 300.0
STAB_PASS = -0.5  # 1/s, spectral abscissa must lie at or left of this
STAB_WITNESS = -1.0
GAMMA_PASS = 1.0
GAMMA_WITNESS = 0.95
METHOD_AGREE_ABS = 0.03
D_SIGMA_MAX = 15.0
EIG_PERT_REL = 1.0e-8
PLANT_SET_ID = "crc-two-gfl-lcl-v1"
SCHEMA_VERSION = "1.0.0"
N_PER = 18
N_CONV = 2
N_X_PHYS = N_PER * N_CONV


def assert_design_bounds() -> None:
    """Liserre placement window: 10*f0 < fres < fsw/2, and Cf reactive < 5%."""
    if not (10.0 * F0 < F_RES < 0.5 * FSW):
        raise RuntimeError(f"LCL resonance {F_RES:.3f} Hz is outside the design window")
    if not (CF_REACTIVE_PU < 0.05):
        raise RuntimeError("capacitor reactive power exceeds the 5% design limit")
    if not (0.01 <= LC_REACTANCE_PU <= 0.10):
        raise RuntimeError("converter reactance is outside the stated 1-10% band")


assert_design_bounds()


@dataclass(frozen=True)
class Case:
    case_id: str
    role: str
    scr: float  # short-circuit ratio at the PCC, base 2*S_RATED
    x_over_r: float
    length_km: tuple[float, float]
    p_pu: tuple[float, float]
    q_pu: tuple[float, float]

    def grid_impedance(self) -> tuple[float, float]:
        """Phase Thevenin R and L. Z_phase = V_ll^2 / S_sc, S_sc = SCR * 2 * S_rated."""
        s_sc = self.scr * (2.0 * S_RATED)
        z_mag = V_LL_RMS**2 / s_sc
        r = z_mag / math.sqrt(1.0 + self.x_over_r**2)
        x = self.x_over_r * r
        return r, x / W0

    def line_rl(self) -> tuple[tuple[float, float], tuple[float, float]]:
        out = []
        for length in self.length_km:
            out.append((CABLE_R_PER_KM * length, CABLE_L_PER_KM * length))
        return out[0], out[1]


CASES: tuple[Case, ...] = (
    Case(
        "nominal",
        "Both inverters export the same moderate power on a strong PCC.",
        scr=10.0,
        x_over_r=4.0,
        length_km=(0.15, 0.15),
        p_pu=(0.55, 0.55),
        q_pu=(0.08, 0.08),
    ),
    Case(
        "weak_grid",
        "Same schedule on a weak, more resistive PCC.",
        scr=3.0,
        x_over_r=1.5,
        length_km=(0.15, 0.15),
        p_pu=(0.55, 0.55),
        q_pu=(0.08, 0.08),
    ),
    Case(
        "high_power",
        "Both inverters near rating and absorbing reactive power, strong PCC.",
        scr=10.0,
        x_over_r=4.0,
        length_km=(0.15, 0.15),
        p_pu=(0.85, 0.85),
        q_pu=(-0.15, -0.15),
    ),
    Case(
        "asymmetric",
        "Unequal schedules and unequal laterals on a moderately strong PCC.",
        scr=6.0,
        x_over_r=3.0,
        length_km=(0.08, 0.80),
        p_pu=(0.70, 0.25),
        q_pu=(0.12, 0.0),
    ),
)


def state_names_si() -> list[str]:
    names: list[str] = []
    for c in (1, 2):
        names.extend(
            [
                f"i_cd_{c}_A",
                f"i_cq_{c}_A",
                f"v_cd_{c}_V",
                f"v_cq_{c}_V",
                f"i_od_{c}_A",
                f"i_oq_{c}_A",
                f"xi_id_{c}_As",
                f"xi_iq_{c}_As",
                f"v_td_{c}_V",
                f"v_tq_{c}_V",
                f"xi_pll_{c}_Vs",
                f"delta_{c}_rad",
                f"v_df_{c}_V",
                f"v_qf_{c}_V",
                f"ic_df_{c}_A",
                f"ic_qf_{c}_A",
                f"io_df_{c}_A",
                f"io_qf_{c}_A",
            ]
        )
    return names


def scaling_diagonal() -> list[float]:
    """x_si = T x_scaled, T diagonal, identical for every case."""
    one = [
        I_BASE,
        I_BASE,
        V_BASE,
        V_BASE,
        I_BASE,
        I_BASE,
        I_BASE / W0,
        I_BASE / W0,
        V_BASE,
        V_BASE,
        V_BASE / W0,
        1.0,
        V_BASE,
        V_BASE,
        I_BASE,
        I_BASE,
        I_BASE,
        I_BASE,
    ]
    return one + one


W_NAMES = [
    "v_gd_5pct",
    "v_gq_5pct",
    "p_ref_1",
    "q_ref_1",
    "p_ref_2",
    "q_ref_2",
    "n_p_1",
    "n_q_1",
    "n_p_2",
    "n_q_2",
]
U_NAMES = ["i_d_ref_1_pu", "i_q_ref_1_pu", "i_d_ref_2_pu", "i_q_ref_2_pu"]
Y_NAMES = ["e_p_1", "e_q_1", "e_p_2", "e_q_2"]
Z_NAMES = [
    "z_ep_1",
    "z_eq_1",
    "z_ep_2",
    "z_eq_2",
    "z_ud_1",
    "z_uq_1",
    "z_ud_2",
    "z_uq_2",
]
WEIGHT_STATE_NAMES = [f"w_{name}" for name in Z_NAMES]
