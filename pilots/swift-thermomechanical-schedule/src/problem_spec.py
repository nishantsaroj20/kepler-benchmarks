"""Public problem constants and the four-variable schedule family.

This module is the shared specification (loading definition and published
numbers). Constitutive integrators do not live here: each solver implements
the return map itself.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# Synthetic truth used only to generate noisy public measurements.
# Steel-like Swift-Voce + Johnson-Cook-type factors from the v1 study.
TRUTH = dict(K=680.0, eps0=0.015, n=0.18, alpha=0.5, Q=220.0, beta=14.0, C=0.012, m=0.8)

E_MPA = 210000.0
NU = 0.3
EDOT0 = 0.1          # 1/s
CHI = 0.9            # Taylor-Quinney factor, known
RHO_C = 3.6          # N/(mm^2 K) = reference volumetric heat capacity
T_REF = 293.0        # K, thermal-softening origin and nominal room temperature
TM = 1800.0          # K
# k such that alpha = k / rho_c = 12.5 mm^2/s (mild-steel diffusivity).
# Units: rho_c [N/(mm^2 K)] * alpha [mm^2/s] = k [N/(s K)].
K_COND = 45.0        # N/(s K)
ALPHA_TH = K_COND / RHO_C  # mm^2/s

# Geometry (half specimen, midplane symmetry). Proposed until the pilot freezes it.
R0 = 3.0             # mm
DR_NOM = 0.015       # mm, midplane radius reduction
L_HALF = 10.0        # mm
Z_GAUGE = 4.0        # mm, material point defining half of the gauge opening
NECK_Z = 1.5         # mm, reference-coordinate neck band [0, NECK_Z]
SHOULDER = (6.0, 8.5)  # mm, reference-coordinate shoulder band

# Fixed end-displacement increments (half-model axial displacement), mm.
# Chosen so the cycle ends past the load peak but before a deep, mesh-fragile neck.
D1 = 0.70
D2 = 0.75
D3 = 0.55

# Actuator limits on the end motion. Inertia at these limits is negligible
# (see specification.md). The smoothstep shape converts them into duration floors.
V_ACT = 8.0          # mm/s
A_ACT = 40.0         # mm/s^2

# Inclusive search box. Lower ramp bounds equal the actuator floors (rounded up
# to 0.01 s). Upper bounds are a process window, not an optimality condition.
# Lower edges are the actuator floors rounded up to 0.01 s (see actuator_floors).
T1_BOX = (0.34, 12.0)
T2_BOX = (0.35, 12.0)
T3_BOX = (0.30, 12.0)
H_BOX = (0.0, 16.0)

# Proposed acceptance limits. The pilot runner treats these as candidates and
# records whether refinement still supports them. They are not a verifier.
T_NECK_CAP_K = 345.0
TAU_TARGET_S = 8.0
R_FLOOR_MM = 2.20
F_CAP_N = 14000.0
GAUGE_WINDOW_MM = (2.40, 2.95)

# Grip temperature, all scenarios. Lateral surface, axis and midplane are insulated.
T_GRIP = 293.0

# Fully public operating scenarios. No hidden conditions.
SCENARIOS = (
    {"id": "S0_nominal", "T_init": 293.0, "dR": 0.015},
    {"id": "S1_warm", "T_init": 333.0, "dR": 0.015},
    {"id": "S2_imperfection", "T_init": 293.0, "dR": 0.030},
)

# Calibration noise (synthetic Gaussian, seeded).
NOISE_FORCE_N = 20.0
NOISE_GAUGE_MM = 0.002
NOISE_TEMP_K = 0.4
NOISE_RADIUS_MM = 0.003
CAL_SEED = 20261108

PARAM_BOUNDS = {
    "K": (300.0, 1500.0),
    "eps0": (0.001, 0.1),
    "n": (0.05, 0.45),
    "alpha": (0.0, 1.0),
    "Q": (0.0, 800.0),
    "beta": (0.5, 80.0),
    "C": (0.0, 0.1),
    "m": (0.3, 3.0),
}


def actuator_floors(d1=D1, d2=D2, d3=D3, v_act=V_ACT, a_act=A_ACT):
    """Minimum ramp durations implied by v_max=2*delta/t and a_max=2*pi*delta/t^2."""
    def floor(delta):
        t_v = 2.0 * delta / v_act
        t_a = math.sqrt(2.0 * math.pi * delta / a_act)
        return t_v, t_a, max(t_v, t_a)
    return {"t1": floor(d1), "t2": floor(d2), "t3": floor(d3)}


def smoothstep(xi):
    """C2 monotone map s(0)=0, s(1)=1, s'=s''=0 at both ends.

    s = xi - sin(2 pi xi) / (2 pi)
    s' = 1 - cos(2 pi xi)          in [0, 2]
    s'' = 2 pi sin(2 pi xi)        in [-2 pi, 2 pi]
    """
    xi = np.asarray(xi, dtype=float)
    s = xi - np.sin(2.0 * math.pi * xi) / (2.0 * math.pi)
    sp = 1.0 - np.cos(2.0 * math.pi * xi)
    spp = 2.0 * math.pi * np.sin(2.0 * math.pi * xi)
    return s, sp, spp


@dataclass(frozen=True)
class Schedule:
    """d = (t1, t2, t3, h). Hold is at the stage-2 / stage-3 boundary.

    Zero hold is continuous: stage 3 starts at the instant stage 2 ends.
    No extra transition interval and no hidden smoothing are applied.
    Displacement increments are fixed; only durations change.
    """

    t1: float
    t2: float
    t3: float
    h: float
    d1: float = D1
    d2: float = D2
    d3: float = D3

    @property
    def tau(self) -> float:
        return float(self.t1 + self.t2 + self.t3 + self.h)

    def stages(self):
        """List of (kind, t_start, duration, u_start, delta)."""
        u1 = self.d1
        u2 = self.d1 + self.d2
        u3 = u2 + self.d3
        out = [
            ("ramp", 0.0, self.t1, 0.0, self.d1),
            ("ramp", self.t1, self.t2, u1, self.d2),
        ]
        t_hold = self.t1 + self.t2
        if self.h > 0.0:
            out.append(("hold", t_hold, self.h, u2, 0.0))
        out.append(("ramp", t_hold + self.h, self.t3, u2, self.d3))
        return out

    def u_end(self, t):
        """Prescribed half-model end displacement (mm) at time t (s)."""
        t = float(t)
        if t <= 0.0:
            return 0.0
        if t >= self.tau:
            return float(self.d1 + self.d2 + self.d3)
        for kind, t0, dur, u0, delta in self.stages():
            if t0 <= t <= t0 + dur + 1e-15:
                if kind == "hold" or dur <= 0.0:
                    return float(u0)
                xi = (t - t0) / dur
                s, _, _ = smoothstep(xi)
                return float(u0 + delta * s)
        return float(self.d1 + self.d2 + self.d3)

    def peaks(self):
        """Exact max |speed| and |acceleration| on each ramp, from the shape."""
        rows = []
        for name, delta, dur in (("t1", self.d1, self.t1), ("t2", self.d2, self.t2), ("t3", self.d3, self.t3)):
            vmax = 2.0 * delta / dur
            amax = 2.0 * math.pi * delta / dur ** 2
            rows.append({"stage": name, "vmax_mm_s": vmax, "amax_mm_s2": amax})
        return rows

    def violations(self):
        """Return reasons the schedule is outside the public box or actuator limits."""
        bad = []
        box = {"t1": (self.t1, T1_BOX), "t2": (self.t2, T2_BOX), "t3": (self.t3, T3_BOX), "h": (self.h, H_BOX)}
        for name, (val, (lo, hi)) in box.items():
            if not (lo - 1e-12 <= val <= hi + 1e-12):
                bad.append(f"{name}={val} outside [{lo}, {hi}]")
        if min(self.t1, self.t2, self.t3) <= 0.0:
            bad.append("ramp durations must be positive")
        if self.h < -1e-15:
            bad.append("hold must be >= 0")
        for row in self.peaks():
            if row["vmax_mm_s"] > V_ACT + 1e-9:
                bad.append(f"{row['stage']} speed {row['vmax_mm_s']:.4f} > {V_ACT}")
            if row["amax_mm_s2"] > A_ACT + 1e-8:
                bad.append(f"{row['stage']} accel {row['amax_mm_s2']:.4f} > {A_ACT}")
        return bad

    def as_dict(self):
        return {
            "t1_s": self.t1, "t2_s": self.t2, "t3_s": self.t3, "h_s": self.h,
            "tau_s": self.tau, "d1_mm": self.d1, "d2_mm": self.d2, "d3_mm": self.d3,
            "legal": not self.violations(), "violations": self.violations(), "peaks": self.peaks(),
        }


def constant_rate_durations(speed_target):
    """Durations proportional to increments so every ramp has the same v_max.

    v_max = 2 delta / t = speed_target => t = 2 delta / speed_target.
    Rejected later if acceleration exceeds the actuator limit.
    """
    v = float(speed_target)
    return Schedule(2.0 * D1 / v, 2.0 * D2 / v, 2.0 * D3 / v, 0.0)
