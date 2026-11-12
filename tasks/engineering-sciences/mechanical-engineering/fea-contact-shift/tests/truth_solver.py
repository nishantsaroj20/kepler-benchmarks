#!/usr/bin/env python3
"""Independent ground truth for approach-controlled frictional contact.

Physical independence from the reference analyzer:
- Influence operator assembled by 4-point Gauss quadrature of the log kernel
  (reference uses the closed-form antiderivative of ln|x-s|).
- Complementarity solved by a separately authored active-set loop in this file
  (no import of solution/analyze.py).

Continuum conventions match the instruction: plane-strain log kernel, prescribed
δ_n/δ_t, Coulomb stick-slip with accumulated plastic slip γ.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

CONTACT_FRAC = 1e-8


def _edges(x: np.ndarray) -> np.ndarray:
    n = len(x)
    e = np.empty(n + 1, dtype=float)
    if n == 1:
        e[0], e[1] = x[0] - 0.5, x[0] + 0.5
        return e
    e[1:-1] = 0.5 * (x[:-1] + x[1:])
    e[0] = x[0] - (e[1] - x[0])
    e[-1] = x[-1] + (x[-1] - e[-2])
    return e


def assemble_kernel(x: np.ndarray, rref: float) -> tuple[np.ndarray, np.ndarray]:
    n = len(x)
    edges = _edges(x)
    w = edges[1:] - edges[:-1]
    # Gauss–Legendre 4-point on each element
    gnodes = np.array([-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526])
    gwt = np.array([0.34785484513745385, 0.6521451548625461, 0.6521451548625461, 0.34785484513745385])
    K = np.zeros((n, n), dtype=float)
    log_r = float(np.log(rref))
    for i, xp in enumerate(x):
        xp = float(xp)
        for j in range(n):
            L, R = float(edges[j]), float(edges[j + 1])
            mid = 0.5 * (L + R)
            half = 0.5 * (R - L)
            acc = 0.0
            for g, wg in zip(gnodes, gwt):
                s = mid + half * g
                acc += wg * (np.log(max(abs(xp - s), 1e-45)) - log_r)
            K[i, j] = acc * half
    return K, w


def _load(pack_dir: Path):
    meta = json.loads((pack_dir / "meta.json").read_text())
    iface = np.loadtxt(pack_dir / "interface.csv", delimiter=",", skiprows=1)
    loads = np.loadtxt(pack_dir / "loads.csv", delimiter=",", skiprows=1)
    if iface.ndim == 1:
        iface = iface.reshape(1, -1)
    if loads.ndim == 1:
        loads = loads.reshape(1, -1)
    return meta, iface[:, 0].astype(float), iface[:, 1].astype(float), loads[:, 1].astype(float), loads[:, 2].astype(float)


def _solve_history(x, g0, A, mu, dn_hist, dt_hist, K, w, max_inner=160):
    n = len(x)
    p = np.zeros(n)
    q = np.zeros(n)
    gamma = np.zeros(n)
    active = np.zeros(n, dtype=bool)
    stick = np.zeros(n, dtype=bool)

    for delta_n, delta_t in zip(dn_hist, dt_hist):
        delta_n = float(delta_n)
        delta_t = float(delta_t)
        interf = g0 - delta_n
        if not np.any(active):
            cand = np.flatnonzero(interf < 0.0)
            if cand.size == 0:
                p[:] = 0.0
                q[:] = 0.0
                active[:] = False
                stick[:] = False
                continue
            active[:] = False
            active[cand] = True
            stick[active] = True

        for _ in range(max_inner):
            idx = np.flatnonzero(active)
            if idx.size == 0:
                break
            st = np.flatnonzero(active & stick)
            sl = np.flatnonzero(active & ~stick)
            m, ms = idx.size, st.size

            slip_sign = np.zeros(n)
            if sl.size:
                pred = A * (K @ q) + gamma
                s = np.sign(pred[sl] - delta_t)
                s[s == 0] = 1.0 if delta_t >= 0 else -1.0
                slip_sign[sl] = s

            nv = m + ms
            M = np.zeros((nv, nv))
            rhs = np.zeros(nv)
            st_map = {int(j): k for k, j in enumerate(st)}

            for a, i in enumerate(idx):
                for b, j in enumerate(idx):
                    M[a, b] = A * K[i, j]
                rhs[a] = g0[i] - delta_n

            for ka, i in enumerate(st):
                row = m + ka
                for b, j in enumerate(idx):
                    if stick[j]:
                        M[row, m + st_map[int(j)]] += A * K[i, j]
                    else:
                        M[row, b] += A * K[i, j] * slip_sign[j] * mu
                rhs[row] = delta_t - gamma[i]

            if nv == 0:
                break
            try:
                sol = np.linalg.solve(M, rhs)
            except np.linalg.LinAlgError:
                sol = np.linalg.lstsq(M, rhs, rcond=None)[0]

            pa = sol[:m]
            qs = sol[m:] if ms else np.zeros(0)
            p_new = np.zeros(n)
            q_new = np.zeros(n)
            p_new[idx] = pa
            if ms:
                q_new[st] = qs
            for i in sl:
                q_new[i] = slip_sign[i] * mu * max(p_new[i], 0.0)

            changed = False
            if np.any(pa < -1e-14):
                active[idx[pa < -1e-14]] = False
                stick[idx[pa < -1e-14]] = False
                changed = True
            p, q = p_new, q_new
            gap = g0 - delta_n - A * (K @ p)
            add = (~active) & (gap < -1e-16)
            if np.any(add):
                active[add] = True
                stick[add] = True
                changed = True
            for i in np.flatnonzero(active & stick):
                lim = mu * max(p[i], 0.0)
                if abs(q[i]) > lim + 1e-9 * max(1.0, lim):
                    stick[i] = False
                    q[i] = np.copysign(lim, q[i] if q[i] != 0 else (delta_t - gamma[i]))
                    changed = True
            if not changed:
                break

        u_el = A * (K @ q)
        for i in range(n):
            if not active[i]:
                gamma[i] = 0.0
            elif not stick[i]:
                gamma[i] = delta_t - u_el[i]

    return p, q, w, active, stick, gamma


def solve_pack_truth(pack_dir: Path) -> dict:
    meta, x, g0, dn, dt = _load(pack_dir)
    E, nu, mu = float(meta["E"]), float(meta["nu"]), float(meta["mu"])
    rref = float(meta.get("rref", 1.0))
    A = 2.0 * (1.0 - nu * nu) / (np.pi * E)
    K, w = assemble_kernel(x, rref)
    p, q, w, active, stick, gamma = _solve_history(x, g0, A, mu, dn, dt, K, w)
    return _metrics(x, p, q, w, stick, gamma, mu)


def _metrics(x, p, q, w, stick, gamma, mu) -> dict:
    pmax = float(np.max(p)) if p.size else 0.0
    contact = p > CONTACT_FRAC * max(pmax, 1e-30)
    slip = contact & (np.abs(q) >= 0.99 * mu * np.maximum(p, 0.0) - 1e-30)
    stick_m = contact & stick & (~slip)
    slip = contact & (~stick_m)
    cl = float(np.sum(w[contact]))
    sl = float(np.sum(w[slip]))
    sf = float(np.sum(w[stick_m]) / cl) if cl > 0 else 0.0
    fn = float(np.sum(p * w))
    ft = float(np.sum(q * w))
    q_peak = float(np.max(np.abs(q))) if q.size else 0.0
    cum = float(np.sum(w * np.abs(gamma)))
    if np.any(contact) and float(np.sum(np.abs(q[contact]) * w[contact])) > 0:
        xq = float(
            np.sum(x[contact] * np.abs(q[contact]) * w[contact])
            / np.sum(np.abs(q[contact]) * w[contact])
        )
    else:
        xq = 0.0
    return {
        "schema_version": 1,
        "p_peak": pmax,
        "fn_total": fn,
        "ft_total": ft,
        "contact_length": cl,
        "q_peak": q_peak,
        "stick_fraction": sf,
        "slip_length": sl,
        "x_q_centroid": xq,
        "cum_slip_l1": cum,
    }


def truth_for_pack(pack_dir: str | Path) -> dict:
    return solve_pack_truth(Path(pack_dir))
