#!/usr/bin/env python3
"""Independent ground truth for frictional IC contact packs.

Same continuum conventions as the instruction; independent numerics:
upsampled interface, load sub-steps, and a separately coded active-set
complementarity loop (not imported from the reference analyzer).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


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


def log_kernel(x: np.ndarray, rref: float) -> tuple[np.ndarray, np.ndarray]:
    n = len(x)
    edges = _edges(x)
    w = edges[1:] - edges[:-1]
    K = np.zeros((n, n), dtype=float)
    lr = float(np.log(rref))
    for i, xi in enumerate(x):
        xi = float(xi)
        for j in range(n):
            L, R = float(edges[j]), float(edges[j + 1])

            def F(z: float) -> float:
                az = abs(z)
                return 0.0 if az < 1e-45 else z * np.log(az) - z

            K[i, j] = (F(xi - L) - F(xi - R)) - lr * (R - L)
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


def _refine(x, g0, factor=2):
    if factor <= 1:
        return x.copy(), g0.copy()
    xn = np.linspace(x[0], x[-1], (len(x) - 1) * factor + 1)
    return xn, np.interp(xn, x, g0)


def _expand_history(fn, ft, substeps: int):
    if substeps <= 1:
        return fn.astype(float), ft.astype(float)
    fn_o = [0.0]
    ft_o = [0.0]
    prev_n, prev_t = 0.0, 0.0
    for a, b in zip(fn, ft):
        for s in range(1, substeps + 1):
            t = s / substeps
            fn_o.append(prev_n + t * (float(a) - prev_n))
            ft_o.append(prev_t + t * (float(b) - prev_t))
        prev_n, prev_t = float(a), float(b)
    return np.asarray(fn_o[1:]), np.asarray(ft_o[1:])


def _solve_step(K, w, A, g0, mu, Fn, Ft, uslip, max_iter=150):
    n = len(g0)
    p = np.zeros(n)
    q = np.zeros(n)
    active = np.zeros(n, dtype=bool)
    stick = np.zeros(n, dtype=bool)
    # seed near minimum gap
    seed_n = max(7, n // 10)
    active[np.argsort(g0)[:seed_n]] = True
    stick[active] = True
    c_n = c_t = 0.0

    for _ in range(max_iter):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            active[np.argsort(g0)[:seed_n]] = True
            stick[:] = False
            stick[active] = True
            idx = np.flatnonzero(active)
        st = np.flatnonzero(active & stick)
        sl = np.flatnonzero(active & ~stick)
        m, ms = idx.size, st.size

        slip_sign = np.zeros(n)
        if sl.size:
            s = np.sign(q[sl])
            default = 1.0 if Ft >= 0 else -1.0
            s[s == 0] = default
            slip_sign[sl] = s

        nv = m + ms + 2
        M = np.zeros((nv, nv))
        rhs = np.zeros(nv)
        st_map = {int(j): k for k, j in enumerate(st)}

        for a, i in enumerate(idx):
            for b, j in enumerate(idx):
                M[a, b] = A * K[i, j]
            M[a, m + ms] = 1.0
            rhs[a] = g0[i]

        for ka, i in enumerate(st):
            row = m + ka
            for b, j in enumerate(idx):
                if stick[j]:
                    M[row, m + st_map[int(j)]] += A * K[i, j]
                else:
                    M[row, b] += A * K[i, j] * slip_sign[j] * mu
            M[row, m + ms + 1] = 1.0
            rhs[row] = -uslip[i]

        rn = m + ms
        rt = m + ms + 1
        for b, j in enumerate(idx):
            M[rn, b] = w[j]
            if stick[j]:
                M[rt, m + st_map[int(j)]] = w[j]
            else:
                M[rt, b] += w[j] * slip_sign[j] * mu
        rhs[rn] = Fn
        rhs[rt] = Ft

        try:
            sol = np.linalg.solve(M, rhs)
        except np.linalg.LinAlgError:
            sol = np.linalg.lstsq(M, rhs, rcond=None)[0]

        pa = sol[:m]
        qs = sol[m : m + ms]
        c_n = float(sol[m + ms])
        c_t = float(sol[m + ms + 1])
        p.fill(0.0)
        q.fill(0.0)
        p[idx] = pa
        q[st] = qs
        q[sl] = slip_sign[sl] * mu * p[sl]

        changed = False
        neg = pa < -1e-14
        if np.any(neg):
            active[idx[neg]] = False
            stick[idx[neg]] = False
            changed = True
        gap = g0 - c_n - A * (K @ p)
        add = (~active) & (gap < -1e-16)
        if np.any(add):
            active[add] = True
            stick[add] = True
            changed = True
        for i in np.flatnonzero(active & stick):
            lim = mu * max(p[i], 0.0)
            if abs(q[i]) > lim + 1e-9 * max(1.0, lim):
                stick[i] = False
                q[i] = np.copysign(lim, q[i] if q[i] != 0 else Ft)
                changed = True
        if not changed:
            break

    uslip_next = -A * (K @ q) - c_t
    return p, q, uslip_next


def solve_pack_truth(pack_dir: Path, upsample_factor: int = 2, substeps: int = 2) -> dict:
    meta, x0, g00, fn0, ft0 = _load(pack_dir)
    x, g0 = _refine(x0, g00, factor=upsample_factor)
    fn, ft = _expand_history(fn0, ft0, substeps=substeps)
    E, nu, mu = float(meta["E"]), float(meta["nu"]), float(meta["mu"])
    rref = float(meta.get("rref", 1.0))
    A = 2.0 * (1.0 - nu * nu) / (np.pi * E)
    K, w = log_kernel(x, rref)
    uslip = np.zeros(len(x))
    p = q = None
    for Fn, Ft in zip(fn, ft):
        p, q, uslip = _solve_step(K, w, A, g0, mu, float(Fn), float(Ft), uslip)
    return _metrics(x, p, q, w, mu)


def _metrics(x, p, q, w, mu) -> dict:
    pmax = float(np.max(p))
    contact = p > 1e-8 * max(pmax, 1e-30)
    cl = float(np.sum(w[contact]))
    fn = float(np.sum(p * w))
    ft = float(np.sum(q * w))
    if np.any(contact):
        p_rms = float(np.sqrt(np.sum((p[contact] ** 2) * w[contact]) / max(cl, 1e-30)))
    else:
        p_rms = 0.0
    xc = float(np.sum(x[contact] * p[contact] * w[contact]) / fn) if cl > 0 and abs(fn) > 0 else 0.0
    mob = float(abs(ft) / (mu * abs(fn))) if abs(fn) > 0 and mu > 0 else 0.0
    return {
        "schema_version": 1,
        "p_peak": pmax,
        "p_rms": p_rms,
        "fn_total": fn,
        "ft_total": ft,
        "contact_length": cl,
        "x_p_centroid": xc,
        "mobilization": mob,
    }


def truth_for_pack(pack_dir: str | Path) -> dict:
    # Upsampled interface + independently coded active-set (not the reference file).
    return solve_pack_truth(Path(pack_dir), upsample_factor=2, substeps=1)
