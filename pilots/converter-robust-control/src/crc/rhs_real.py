"""Path B: real-arithmetic component stamp of the averaged dq model.

The output-branch / PCC constraint is eliminated inside the right-hand
side by a 6-by-6 real linear solve. This file does not import the
complex transcription or the analytic Jacobian.
"""

from __future__ import annotations

import numpy as np

from .indexing import DELTA, IC, ICF, IO, IOF, N_PER, PLL, VC, VF, VT, XI, base
from .params import (
    CF,
    KI_I,
    KI_PLL,
    KP_I,
    KP_PLL,
    LC,
    LG_FILTER,
    RC,
    RD,
    RG_FILTER,
    TAU_A,
    TAU_SENSOR,
    W0,
    Case,
)


def branch_rl(case: Case) -> tuple[float, float, float, float, float, float]:
    rg, lg = case.grid_impedance()
    (r1, l1), (r2, l2) = case.line_rl()
    lo1 = LG_FILTER + l1
    ro1 = RG_FILTER + r1
    lo2 = LG_FILTER + l2
    ro2 = RG_FILTER + r2
    return lo1, ro1, lo2, ro2, lg, rg


def _solve_branches(
    vf1: np.ndarray,
    io1: np.ndarray,
    vf2: np.ndarray,
    io2: np.ndarray,
    vg: np.ndarray,
    lo1: float,
    ro1: float,
    lo2: float,
    ro2: float,
    lg: float,
    rg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unknowns: di_o1, di_o2, v_pcc. M is constant in the state."""
    mat = np.zeros((6, 6), dtype=float)
    mat[0, 0] = lo1
    mat[0, 4] = 1.0
    mat[1, 1] = lo1
    mat[1, 5] = 1.0
    mat[2, 2] = lo2
    mat[2, 4] = 1.0
    mat[3, 3] = lo2
    mat[3, 5] = 1.0
    mat[4, 0] = lg
    mat[4, 2] = lg
    mat[4, 4] = -1.0
    mat[5, 1] = lg
    mat[5, 3] = lg
    mat[5, 5] = -1.0
    ig = io1 + io2
    rhs = np.array(
        [
            vf1[0] - ro1 * io1[0] + W0 * lo1 * io1[1],
            vf1[1] - ro1 * io1[1] - W0 * lo1 * io1[0],
            vf2[0] - ro2 * io2[0] + W0 * lo2 * io2[1],
            vf2[1] - ro2 * io2[1] - W0 * lo2 * io2[0],
            -vg[0] - rg * ig[0] + W0 * lg * ig[1],
            -vg[1] - rg * ig[1] - W0 * lg * ig[0],
        ],
        dtype=float,
    )
    sol = np.linalg.solve(mat, rhs)
    return sol[0:2], sol[2:4], sol[4:6]


def _rot_into_pll(delta: float, vec: np.ndarray) -> np.ndarray:
    c = math_cos(delta)
    s = math_sin(delta)
    return np.array([c * vec[0] + s * vec[1], -s * vec[0] + c * vec[1]], dtype=float)


def _rot_into_grid(delta: float, vec: np.ndarray) -> np.ndarray:
    c = math_cos(delta)
    s = math_sin(delta)
    return np.array([c * vec[0] - s * vec[1], s * vec[0] + c * vec[1]], dtype=float)


def math_cos(delta: float) -> float:
    return float(np.cos(delta))


def math_sin(delta: float) -> float:
    return float(np.sin(delta))


def rhs_real(x: np.ndarray, u_amp: np.ndarray, vg: np.ndarray, case: Case) -> np.ndarray:
    """State derivative. u_amp is the total dq current reference in amperes."""
    lo1, ro1, lo2, ro2, lg, rg = branch_rl(case)
    f = np.zeros(2 * N_PER, dtype=float)
    blocks = []
    for conv in (0, 1):
        b = base(conv)
        ic = x[b + IC : b + IC + 2]
        vc = x[b + VC : b + VC + 2]
        io = x[b + IO : b + IO + 2]
        icap = ic - io
        vf = vc + RD * icap
        blocks.append((ic, vc, io, icap, vf))
    dio1, dio2, _vpcc = _solve_branches(
        blocks[0][4], blocks[0][2], blocks[1][4], blocks[1][2], vg, lo1, ro1, lo2, ro2, lg, rg
    )
    dios = (dio1, dio2)
    for conv in (0, 1):
        b = base(conv)
        ic, vc, io, icap, vf = blocks[conv]
        delta = float(x[b + DELTA])
        vdf = x[b + VF : b + VF + 2]
        icf = x[b + ICF : b + ICF + 2]
        iof = x[b + IOF : b + IOF + 2]
        xi = x[b + XI : b + XI + 2]
        vt = x[b + VT : b + VT + 2]
        xi_pll = float(x[b + PLL])
        uref = u_amp[2 * conv : 2 * conv + 2]

        f[b + IC] = (vt[0] - vf[0] - RC * ic[0] + W0 * LC * ic[1]) / LC
        f[b + IC + 1] = (vt[1] - vf[1] - RC * ic[1] - W0 * LC * ic[0]) / LC
        f[b + VC] = (icap[0] + W0 * CF * vc[1]) / CF
        f[b + VC + 1] = (icap[1] - W0 * CF * vc[0]) / CF
        f[b + IO : b + IO + 2] = dios[conv]

        vf_p = _rot_into_pll(delta, vf)
        ic_p = _rot_into_pll(delta, ic)
        io_p = _rot_into_pll(delta, io)
        f[b + VF : b + VF + 2] = (vf_p - vdf) / TAU_SENSOR
        f[b + ICF : b + ICF + 2] = (ic_p - icf) / TAU_SENSOR
        f[b + IOF : b + IOF + 2] = (io_p - iof) / TAU_SENSOR

        f[b + XI] = uref[0] - icf[0]
        f[b + XI + 1] = uref[1] - icf[1]
        vd_cmd = KP_I * (uref[0] - icf[0]) + KI_I * xi[0] - W0 * LC * icf[1] + vdf[0]
        vq_cmd = KP_I * (uref[1] - icf[1]) + KI_I * xi[1] + W0 * LC * icf[0] + vdf[1]
        vcmd_g = _rot_into_grid(delta, np.array([vd_cmd, vq_cmd], dtype=float))
        f[b + VT : b + VT + 2] = (vcmd_g - vt) / TAU_A

        f[b + PLL] = vdf[1]
        f[b + DELTA] = KP_PLL * vdf[1] + KI_PLL * xi_pll
    return f


def terminal_voltage(x: np.ndarray, conv: int) -> np.ndarray:
    b = base(conv)
    ic = x[b + IC : b + IC + 2]
    vc = x[b + VC : b + VC + 2]
    io = x[b + IO : b + IO + 2]
    return vc + RD * (ic - io)


def powers_from_filters(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Filtered PLL-frame power. Positive Q is leading-current (var injection)."""
    p = np.zeros(2)
    q = np.zeros(2)
    for conv in (0, 1):
        b = base(conv)
        vd = float(x[b + VF])
        vq = float(x[b + VF + 1])
        iod = float(x[b + IOF])
        ioq = float(x[b + IOF + 1])
        p[conv] = 1.5 * (vd * iod + vq * ioq)
        q[conv] = 1.5 * (vd * ioq - vq * iod)
    return p, q


def grid_current(x: np.ndarray) -> np.ndarray:
    i1 = x[base(0) + IO : base(0) + IO + 2]
    i2 = x[base(1) + IO : base(1) + IO + 2]
    return i1 + i2
