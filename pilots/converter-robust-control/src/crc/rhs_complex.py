"""Path A residual: complex dq transcription of the same averaged circuit.

Independent of rhs_real. The branch elimination uses a 3-by-3 complex
matrix rather than the 6-by-6 real stamp. Cross-coupling is the -j*w*L
term. Frame rotations use exp(-j*delta).
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


def _rl(case: Case) -> tuple[complex, complex, complex, complex, complex, complex]:
    rg, lg = case.grid_impedance()
    (r1, l1), (r2, l2) = case.line_rl()
    z_w = 1j * W0
    lo1 = LG_FILTER + l1
    lo2 = LG_FILTER + l2
    return lo1, RG_FILTER + r1, lo2, RG_FILTER + r2, lg, rg


def _as_c(vec: np.ndarray) -> complex:
    return complex(vec[0], vec[1])


def _as_r(val: complex) -> np.ndarray:
    return np.array([val.real, val.imag], dtype=float)


def rhs_complex(x: np.ndarray, u_amp: np.ndarray, vg: np.ndarray, case: Case) -> np.ndarray:
    lo1, ro1, lo2, ro2, lg, rg = _rl(case)
    f = np.zeros(2 * N_PER, dtype=float)
    vfs = []
    ios = []
    ics = []
    vcs = []
    for conv in (0, 1):
        b = base(conv)
        ic = _as_c(x[b + IC : b + IC + 2])
        vc = _as_c(x[b + VC : b + VC + 2])
        io = _as_c(x[b + IO : b + IO + 2])
        vf = vc + RD * (ic - io)
        vfs.append(vf)
        ios.append(io)
        ics.append(ic)
        vcs.append(vc)
    # Complex branch law:
    # Lo di_o/dt = vf - vpcc - Ro io - j w Lo io
    # Lg d(io1+io2)/dt = vpcc - vg - Rg ig - j w Lg ig
    mat = np.array(
        [
            [lo1, 0.0, 1.0],
            [0.0, lo2, 1.0],
            [lg, lg, -1.0],
        ],
        dtype=complex,
    )
    ig = ios[0] + ios[1]
    vg_c = _as_c(vg)
    rhs = np.array(
        [
            vfs[0] - ro1 * ios[0] - 1j * W0 * lo1 * ios[0],
            vfs[1] - ro2 * ios[1] - 1j * W0 * lo2 * ios[1],
            -vg_c - rg * ig - 1j * W0 * lg * ig,
        ],
        dtype=complex,
    )
    dio1, dio2, _vp = np.linalg.solve(mat, rhs)
    dios = (dio1, dio2)
    for conv in (0, 1):
        b = base(conv)
        delta = float(x[b + DELTA])
        rot = np.exp(-1j * delta)  # grid phasor -> PLL phasor
        vf_p = rot * vfs[conv]
        ic_p = rot * ics[conv]
        io_p = rot * ios[conv]
        vdf = _as_c(x[b + VF : b + VF + 2])
        icf = _as_c(x[b + ICF : b + ICF + 2])
        iof = _as_c(x[b + IOF : b + IOF + 2])
        xi = _as_c(x[b + XI : b + XI + 2])
        vt = _as_c(x[b + VT : b + VT + 2])
        uref = _as_c(u_amp[2 * conv : 2 * conv + 2])
        xi_pll = float(x[b + PLL])

        dic = (vt - vfs[conv] - RC * ics[conv] - 1j * W0 * LC * ics[conv]) / LC
        dvc = ((ics[conv] - ios[conv]) - 1j * W0 * CF * vcs[conv]) / CF
        f[b + IC : b + IC + 2] = _as_r(dic)
        f[b + VC : b + VC + 2] = _as_r(dvc)
        f[b + IO : b + IO + 2] = _as_r(dios[conv])
        f[b + VF : b + VF + 2] = _as_r((vf_p - vdf) / TAU_SENSOR)
        f[b + ICF : b + ICF + 2] = _as_r((ic_p - icf) / TAU_SENSOR)
        f[b + IOF : b + IOF + 2] = _as_r((io_p - iof) / TAU_SENSOR)
        f[b + XI : b + XI + 2] = _as_r(uref - icf)

        # Decoupling in the PLL frame uses fixed w0. Feedforward is filtered vf.
        # +j*w*Lc*i cancels the plant's -j*w*Lc*i cross coupling.
        v_cmd = KP_I * (uref - icf) + KI_I * xi + 1j * W0 * LC * icf + vdf
        v_cmd_g = np.exp(1j * delta) * v_cmd
        f[b + VT : b + VT + 2] = _as_r((v_cmd_g - vt) / TAU_A)
        f[b + PLL] = float(np.imag(vdf))
        f[b + DELTA] = KP_PLL * float(np.imag(vdf)) + KI_PLL * xi_pll
    return f
