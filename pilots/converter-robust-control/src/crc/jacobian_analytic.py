"""Analytic linearization of the real dq stamp.

The network map is the inverse of the constant 6-by-6 branch matrix.
Rotation, PI, PLL and sensor partial derivatives are written out.
This module does not call finite differences and does not import rhs_real
or rhs_complex.
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
    N_X_PHYS,
    RC,
    RD,
    RG_FILTER,
    TAU_A,
    TAU_SENSOR,
    W0,
    Case,
)


def _branch(case: Case):
    rg, lg = case.grid_impedance()
    (r1, l1), (r2, l2) = case.line_rl()
    return LG_FILTER + l1, RG_FILTER + r1, LG_FILTER + l2, RG_FILTER + r2, lg, rg


def _branch_matrix(lo1, lo2, lg) -> np.ndarray:
    mat = np.zeros((6, 6))
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
    return mat


def linearize_analytic(
    x: np.ndarray, u_amp: np.ndarray, vg: np.ndarray, case: Case
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lo1, ro1, lo2, ro2, lg, rg = _branch(case)
    mat = _branch_matrix(lo1, lo2, lg)
    # Columns of rhs sensitivity w.r.t.
    # vf1(2), io1(2), vf2(2), io2(2), vg(2). Direct, before vf(io) chain.
    rhs_s = np.zeros((6, 10))
    # vf1
    rhs_s[0, 0] = 1.0
    rhs_s[1, 1] = 1.0
    # io1: +w lo on the orthogonal axis, -ro on the same axis
    rhs_s[0, 2] = -ro1
    rhs_s[0, 3] = W0 * lo1
    rhs_s[1, 2] = -W0 * lo1
    rhs_s[1, 3] = -ro1
    # vf2
    rhs_s[2, 4] = 1.0
    rhs_s[3, 5] = 1.0
    rhs_s[2, 6] = -ro2
    rhs_s[2, 7] = W0 * lo2
    rhs_s[3, 6] = -W0 * lo2
    rhs_s[3, 7] = -ro2
    # vg and ig = io1+io2
    rhs_s[4, 8] = -1.0
    rhs_s[5, 9] = -1.0
    rhs_s[4, 2] += -rg
    rhs_s[4, 3] += W0 * lg
    rhs_s[5, 2] += -W0 * lg
    rhs_s[5, 3] += -rg
    rhs_s[4, 6] += -rg
    rhs_s[4, 7] += W0 * lg
    rhs_s[5, 6] += -W0 * lg
    rhs_s[5, 7] += -rg
    sol_s = np.linalg.solve(mat, rhs_s)  # 6 x 10
    # dio_c / d(packed)
    dio_s = (sol_s[0:2, :], sol_s[2:4, :])

    eye = np.eye(2)
    jw = np.array([[0.0, 1.0], [-1.0, 0.0]])
    a = np.zeros((N_X_PHYS, N_X_PHYS))
    bu = np.zeros((N_X_PHYS, 4))
    bv = np.zeros((N_X_PHYS, 2))

    vfs = []
    ics = []
    ios = []
    for conv in (0, 1):
        b = base(conv)
        ic = x[b + IC : b + IC + 2]
        vc = x[b + VC : b + VC + 2]
        io = x[b + IO : b + IO + 2]
        vf = vc + RD * (ic - io)
        vfs.append(vf)
        ics.append(ic)
        ios.append(io)

    pack_slice = {
        0: {"vf": slice(0, 2), "io": slice(2, 4)},
        1: {"vf": slice(4, 6), "io": slice(6, 8)},
    }

    for conv in (0, 1):
        b = base(conv)
        delta = float(x[b + DELTA])
        c = float(np.cos(delta))
        s = float(np.sin(delta))
        rot = np.array([[c, s], [-s, c]])
        drot = np.array([[-s, c], [-c, -s]])
        rot_back = np.array([[c, -s], [s, c]])
        drot_back = np.array([[-s, -c], [c, -s]])
        vf = vfs[conv]
        ic = ics[conv]
        io = ios[conv]
        vdf = x[b + VF : b + VF + 2]
        icf = x[b + ICF : b + ICF + 2]
        xi = x[b + XI : b + XI + 2]
        uref = u_amp[2 * conv : 2 * conv + 2]

        # Converter inductor and capacitor, through vf = vc + Rd (ic - io).
        # dic/dvt, dic/dic, dic/dvc, dic/dio
        a[b + IC : b + IC + 2, b + VT : b + VT + 2] = eye / LC
        a[b + IC : b + IC + 2, b + IC : b + IC + 2] = -(RC / LC) * eye + W0 * jw - (RD / LC) * eye
        a[b + IC : b + IC + 2, b + VC : b + VC + 2] = -eye / LC
        a[b + IC : b + IC + 2, b + IO : b + IO + 2] = (RD / LC) * eye
        a[b + VC : b + VC + 2, b + IC : b + IC + 2] = eye / CF
        a[b + VC : b + VC + 2, b + IO : b + IO + 2] = -eye / CF
        a[b + VC : b + VC + 2, b + VC : b + VC + 2] = W0 * jw

        # Output-current derivatives from the branch solve.
        sl = pack_slice[conv]
        # direct packed columns plus vf's dependence on vc, ic, io
        for src_conv in (0, 1):
            sb = base(src_conv)
            sls = pack_slice[src_conv]
            dio_vf = dio_s[conv][:, sls["vf"]]  # 2x2
            dio_io = dio_s[conv][:, sls["io"]]
            # vf = vc + Rd ic - Rd io
            a[b + IO : b + IO + 2, sb + VC : sb + VC + 2] += dio_vf
            a[b + IO : b + IO + 2, sb + IC : sb + IC + 2] += dio_vf * RD
            a[b + IO : b + IO + 2, sb + IO : sb + IO + 2] += dio_io - dio_vf * RD
        bv[b + IO : b + IO + 2, :] = dio_s[conv][:, 8:10]

        # Sensors.
        for offset, vec in ((VF, vf), (ICF, ic), (IOF, io)):
            a[b + offset : b + offset + 2, b + offset : b + offset + 2] = -eye / TAU_SENSOR
            a[b + offset : b + offset + 2, b + DELTA] = (drot @ vec) / TAU_SENSOR
            if offset == VF:
                gain = rot / TAU_SENSOR
                a[b + VF : b + VF + 2, b + VC : b + VC + 2] += gain
                a[b + VF : b + VF + 2, b + IC : b + IC + 2] += gain * RD
                a[b + VF : b + VF + 2, b + IO : b + IO + 2] += -gain * RD
            elif offset == ICF:
                a[b + ICF : b + ICF + 2, b + IC : b + IC + 2] += rot / TAU_SENSOR
            else:
                a[b + IOF : b + IOF + 2, b + IO : b + IO + 2] += rot / TAU_SENSOR

        # Current integrators and actuator.
        a[b + XI : b + XI + 2, b + ICF : b + ICF + 2] = -eye
        bu[b + XI : b + XI + 2, 2 * conv : 2 * conv + 2] = eye
        dvcmd_dicf = np.array([[-KP_I, -W0 * LC], [W0 * LC, -KP_I]])
        vd_cmd = KP_I * (uref[0] - icf[0]) + KI_I * xi[0] - W0 * LC * icf[1] + vdf[0]
        vq_cmd = KP_I * (uref[1] - icf[1]) + KI_I * xi[1] + W0 * LC * icf[0] + vdf[1]
        vcmd = np.array([vd_cmd, vq_cmd])
        # vt_dot depends on vcmd_g(delta, vcmd(states,u)) and vt
        # d vt_dot / d alpha = rot_back @ dvcmd/d alpha / tau
        a[b + VT : b + VT + 2, b + VT : b + VT + 2] = -eye / TAU_A
        a[b + VT : b + VT + 2, b + DELTA] = (drot_back @ vcmd) / TAU_A
        chain = rot_back / TAU_A
        a[b + VT : b + VT + 2, b + ICF : b + ICF + 2] += chain @ dvcmd_dicf
        a[b + VT : b + VT + 2, b + VF : b + VF + 2] += chain
        a[b + VT : b + VT + 2, b + XI : b + XI + 2] += chain * KI_I
        bu[b + VT : b + VT + 2, 2 * conv : 2 * conv + 2] += chain * KP_I

        # PLL. Uses the filtered q voltage, not the raw terminal voltage.
        a[b + PLL, b + VF + 1] = 1.0
        a[b + DELTA, b + VF + 1] = KP_PLL
        a[b + DELTA, b + PLL] = KI_PLL

    return a, bu, bv
