"""Compare the two residuals and the analytic Jacobian against finite differences."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crc.jacobian_analytic import linearize_analytic
from crc.params import CASES, I_BASE, V_BASE
from crc.rhs_complex import rhs_complex
from crc.rhs_real import rhs_real


def main() -> None:
    case = CASES[0]
    rng = np.random.default_rng(7)
    x = np.zeros(36)
    x[0] = 80.0
    x[1] = -10.0
    x[2] = V_BASE
    x[4] = 70.0
    x[5] = -8.0
    x[8] = V_BASE
    x[11] = 0.04
    x[12] = V_BASE
    x[16] = 70.0
    x[18] = 40.0
    x[20] = V_BASE * 0.98
    x[22] = 35.0
    x[26] = V_BASE * 0.98
    x[29] = -0.03
    x[30] = V_BASE * 0.98
    x[34] = 35.0
    u = np.array([80.0, -10.0, 40.0, -5.0])
    vg = np.array([V_BASE, 0.0])
    fr = rhs_real(x, u, vg, case)
    fc = rhs_complex(x, u, vg, case)
    err = np.max(np.abs(fr - fc))
    scale = np.max(np.abs(fr)) + 1.0
    print(f"rhs max abs diff {err:.3e} rel {err / scale:.3e}")

    a, bu, bv = linearize_analytic(x, u, vg, case)

    def f_x(z):
        return rhs_real(z, u, vg, case)

    fd = np.zeros_like(a)
    for i in range(36):
        step = 1e-6 * max(1.0, abs(x[i]))
        xp = x.copy()
        xm = x.copy()
        xp[i] += step
        xm[i] -= step
        fd[:, i] = (f_x(xp) - f_x(xm)) / (2.0 * step)
    diff = a - fd
    print(f"A vs FD max {np.max(np.abs(diff)):.3e} rel {np.max(np.abs(diff)) / (np.max(np.abs(fd)) + 1):.3e}")
    worst = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
    print(f"worst A entry {worst} analytic {a[worst]:.6e} fd {fd[worst]:.6e}")

    fd_u = np.zeros_like(bu)
    for i in range(4):
        step = 1e-6 * I_BASE
        up = u.copy()
        um = u.copy()
        up[i] += step
        um[i] -= step
        fd_u[:, i] = (rhs_real(x, up, vg, case) - rhs_real(x, um, vg, case)) / (2.0 * step)
    print(f"Bu vs FD max {np.max(np.abs(bu - fd_u)):.3e}")

    fd_v = np.zeros_like(bv)
    for i in range(2):
        step = 1e-6 * V_BASE
        vp = vg.copy()
        vm = vg.copy()
        vp[i] += step
        vm[i] -= step
        fd_v[:, i] = (rhs_real(x, u, vp, case) - rhs_real(x, u, vm, case)) / (2.0 * step)
    print(f"Bv vs FD max {np.max(np.abs(bv - fd_v)):.3e}")
    # A few step sizes
    for eps in (1e-4, 1e-5, 1e-6, 1e-7):
        fd = np.zeros_like(a)
        for i in range(36):
            step = eps * max(1.0, abs(x[i]))
            xp = x.copy()
            xm = x.copy()
            xp[i] += step
            xm[i] -= step
            fd[:, i] = (f_x(xp) - f_x(xm)) / (2.0 * step)
        rel = np.max(np.abs(a - fd)) / (np.max(np.abs(a)) + 1.0)
        print(f"step {eps:.1e} rel {rel:.3e}")


if __name__ == "__main__":
    main()
