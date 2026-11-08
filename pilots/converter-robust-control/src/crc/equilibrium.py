"""Nonlinear equilibria of the averaged model.

The unknown vector is the physical state together with the four total
current references in amperes. Integral states are free; power schedules
are equality constraints. The grid-frame angle is fixed by putting the
Thevenin voltage on the positive d-axis.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares, root

from .indexing import DELTA, IC, ICF, IO, IOF, PLL, VC, VF, VT, XI, base
from .params import CF, I_BASE, KI_I, LC, RC, RD, S_RATED, V_BASE, W0, Case
from .rhs_real import grid_current, powers_from_filters, rhs_real, terminal_voltage


def _scales() -> np.ndarray:
    w = W0
    one = np.array(
        [
            I_BASE * w,
            I_BASE * w,
            V_BASE * w,
            V_BASE * w,
            I_BASE * w,
            I_BASE * w,
            I_BASE,
            I_BASE,
            V_BASE / 1.5e-4,
            V_BASE / 1.5e-4,
            V_BASE,
            1.0,
            V_BASE * w,
            V_BASE * w,
            I_BASE * w,
            I_BASE * w,
            I_BASE * w,
            I_BASE * w,
        ]
    )
    return np.concatenate([one, one])


def initial_guess(case: Case) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.zeros(36)
    u = np.zeros(4)
    vg = np.array([V_BASE, 0.0])
    for conv in (0, 1):
        b = base(conv)
        p = case.p_pu[conv] * S_RATED
        q = case.q_pu[conv] * S_RATED
        vd = V_BASE
        iod = p / (1.5 * vd)
        ioq = q / (1.5 * vd)
        vcd, vcq = vd, 0.0
        icap_d = -W0 * CF * vcq
        icap_q = W0 * CF * vcd
        icd = iod + icap_d
        icq = ioq + icap_q
        vf_d = vcd + RD * icap_d
        vf_q = vcq + RD * icap_q
        vtd = vf_d + RC * icd - W0 * LC * icq
        vtq = vf_q + RC * icq + W0 * LC * icd
        # Small positive angle for exporting active power into an inductive grid.
        delta = 0.03 * case.p_pu[conv]
        c, s = np.cos(delta), np.sin(delta)
        rot = np.array([[c, s], [-s, c]])
        x[b + IC : b + IC + 2] = (icd, icq)
        x[b + VC : b + VC + 2] = (vcd, vcq)
        x[b + IO : b + IO + 2] = (iod, ioq)
        x[b + VT : b + VT + 2] = (vtd, vtq)
        x[b + DELTA] = delta
        x[b + PLL] = 0.0
        vf_p = rot @ np.array([vf_d, vf_q])
        ic_p = rot @ np.array([icd, icq])
        io_p = rot @ np.array([iod, ioq])
        x[b + VF : b + VF + 2] = vf_p
        x[b + ICF : b + ICF + 2] = ic_p
        x[b + IOF : b + IOF + 2] = io_p
        # Grid-frame voltage guess rotated back is only approximate; start the
        # integrators from the PLL-frame resistive drop and let the solver finish.
        x[b + XI] = (RC * ic_p[0]) / KI_I
        x[b + XI + 1] = (RC * ic_p[1]) / KI_I
        u[2 * conv : 2 * conv + 2] = ic_p
    return x, u, vg


def residual_vector(z: np.ndarray, case: Case, vg: np.ndarray) -> np.ndarray:
    x = z[:36]
    u = z[36:]
    f = rhs_real(x, u, vg, case)
    p, q = powers_from_filters(x)
    pe = np.array(
        [
            (p[0] - case.p_pu[0] * S_RATED) / S_RATED,
            (q[0] - case.q_pu[0] * S_RATED) / S_RATED,
            (p[1] - case.p_pu[1] * S_RATED) / S_RATED,
            (q[1] - case.q_pu[1] * S_RATED) / S_RATED,
        ]
    )
    return np.concatenate([f / _scales(), pe])


def solve_equilibrium(case: Case) -> dict:
    x0, u0, vg = initial_guess(case)
    z0 = np.concatenate([x0, u0])

    def fun(z: np.ndarray) -> np.ndarray:
        return residual_vector(z, case, vg)

    fit = least_squares(fun, z0, method="trf", xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=400)
    if not fit.success or np.max(np.abs(fit.fun)) > 1e-8:
        alt = root(fun, fit.x if fit.x is not None else z0, method="lm", tol=1e-12)
        if alt.success and np.max(np.abs(alt.fun)) < np.max(np.abs(fit.fun)):
            z = alt.x
            message = str(alt.message)
            success = bool(alt.success)
            nfev = int(alt.nfev)
            final = alt.fun
        else:
            z = fit.x
            message = fit.message
            success = bool(fit.success)
            nfev = int(fit.nfev)
            final = fit.fun
    else:
        z = fit.x
        message = fit.message
        success = True
        nfev = int(fit.nfev)
        final = fit.fun
    x = z[:36]
    u = z[36:]
    f = rhs_real(x, u, vg, case)
    p, q = powers_from_filters(x)
    return {
        "x": x,
        "u_amp": u,
        "vg": vg,
        "success": bool(success) and float(np.max(np.abs(f))) < 1e-6,
        "message": message,
        "nfev": nfev,
        "residual_scaled_inf": float(np.max(np.abs(final))),
        "rhs_inf_si": float(np.max(np.abs(f))),
        "p_W": p,
        "q_var": q,
        "balance": power_balance(x, vg, case),
    }


def power_balance(x: np.ndarray, vg: np.ndarray, case: Case) -> dict:
    """Converter terminal power, grid-source absorption and ohmic loss."""
    from .params import RD, RC, RG_FILTER, LG_FILTER

    (r1, _l1), (r2, _l2) = case.line_rl()
    rg, _lg = case.grid_impedance()
    p_conv = 0.0
    loss = 0.0
    for conv, ro_extra in ((0, r1), (1, r2)):
        b = base(conv)
        vt = x[b + VT : b + VT + 2]
        ic = x[b + IC : b + IC + 2]
        io = x[b + IO : b + IO + 2]
        icap = ic - io
        p_conv += 1.5 * float(vt @ ic)
        loss += 1.5 * RC * float(ic @ ic)
        loss += 1.5 * RD * float(icap @ icap)
        loss += 1.5 * (RG_FILTER + ro_extra) * float(io @ io)
    ig = grid_current(x)
    loss += 1.5 * rg * float(ig @ ig)
    p_grid = 1.5 * float(vg @ ig)
    return {
        "p_converter_W": p_conv,
        "p_grid_absorbed_W": p_grid,
        "p_loss_W": loss,
        "imbalance_W": p_conv - p_grid - loss,
        "v_t_peak_V": [
            float(np.hypot(x[base(c) + VT], x[base(c) + VT + 1])) for c in (0, 1)
        ],
        "delta_rad": [float(x[base(c) + DELTA]) for c in (0, 1)],
        "v_f": [terminal_voltage(x, c).tolist() for c in (0, 1)],
    }
