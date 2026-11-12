#!/usr/bin/env python3
"""Frictional multi-patch contact diagnostics on 2D plane-strain IC packs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

SCHEMA_VERSION = 1
EPS = 1e-15


def _load_pack(pack_dir: Path):
    meta = json.loads((pack_dir / "meta.json").read_text())
    interface = np.loadtxt(pack_dir / "interface.csv", delimiter=",", skiprows=1)
    loads = np.loadtxt(pack_dir / "loads.csv", delimiter=",", skiprows=1)
    if interface.ndim == 1:
        interface = interface.reshape(1, -1)
    if loads.ndim == 1:
        loads = loads.reshape(1, -1)
    x = interface[:, 0].astype(float)
    g0 = interface[:, 1].astype(float)
    fn = loads[:, 1].astype(float)
    ft = loads[:, 2].astype(float)
    return meta, x, g0, fn, ft


def _edges_from_nodes(x: np.ndarray) -> np.ndarray:
    n = len(x)
    edges = np.empty(n + 1, dtype=float)
    if n == 1:
        # degenerate; width from meta not available — unit width
        edges[0], edges[1] = x[0] - 0.5, x[0] + 0.5
        return edges
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - (edges[1] - x[0])
    edges[-1] = x[-1] + (x[-1] - edges[-2])
    return edges


def influence_matrix(x: np.ndarray, rref: float) -> tuple[np.ndarray, np.ndarray]:
    """Integrated log kernel for piecewise-constant tractions (plane-strain)."""
    n = len(x)
    edges = _edges_from_nodes(x)
    w = edges[1:] - edges[:-1]
    K = np.zeros((n, n), dtype=float)
    log_r = float(np.log(rref))
    for i in range(n):
        xi = x[i]
        for j in range(n):
            L, R = edges[j], edges[j + 1]

            def F(z: float) -> float:
                az = abs(z)
                if az < 1e-45:
                    return 0.0
                return z * np.log(az) - z

            K[i, j] = (F(xi - L) - F(xi - R)) - log_r * (R - L)
    return K, w


def _kernel_coeff(E: float, nu: float) -> float:
    return 2.0 * (1.0 - nu * nu) / (np.pi * E)


def solve_history(
    x: np.ndarray,
    g0: np.ndarray,
    A: float,
    mu: float,
    fn_hist: np.ndarray,
    ft_hist: np.ndarray,
    rref: float,
    max_inner: int = 120,
) -> dict:
    K, w = influence_matrix(x, rref=rref)
    n = len(x)
    p = np.zeros(n)
    q = np.zeros(n)
    uslip = np.zeros(n)
    active = np.zeros(n, dtype=bool)
    stick = np.zeros(n, dtype=bool)

    for Fn, Ft in zip(fn_hist, ft_hist):
        Fn = float(Fn)
        Ft = float(Ft)
        if not np.any(active):
            seed = max(5, n // 12)
            active[np.argsort(g0)[:seed]] = True
            stick[active] = True

        c_n = 0.0
        c_t = 0.0
        for _ in range(max_inner):
            idx = np.flatnonzero(active)
            if idx.size == 0:
                # open all; try seeding again under load
                seed = max(5, n // 12)
                active[np.argsort(g0)[:seed]] = True
                stick[active] = True
                idx = np.flatnonzero(active)
            stick_idx = np.flatnonzero(active & stick)
            slip_idx = np.flatnonzero(active & ~stick)
            m = idx.size
            ms = stick_idx.size

            slip_sign = np.zeros(n)
            if slip_idx.size:
                s = np.sign(q[slip_idx])
                s[s == 0] = 1.0 if Ft >= 0 else -1.0
                slip_sign[slip_idx] = s

            nv = m + ms + 2
            M = np.zeros((nv, nv))
            rhs = np.zeros(nv)

            # normal contact equations
            for a, i in enumerate(idx):
                for b, j in enumerate(idx):
                    M[a, b] = A * K[i, j]
                M[a, m + ms] = 1.0
                rhs[a] = g0[i]

            # stick tangential equations
            stick_local = {int(j): s for s, j in enumerate(stick_idx)}
            for s_a, i in enumerate(stick_idx):
                row = m + s_a
                for b, j in enumerate(idx):
                    if stick[j]:
                        M[row, m + stick_local[int(j)]] += A * K[i, j]
                    else:
                        M[row, b] += A * K[i, j] * (slip_sign[j] * mu)
                M[row, m + ms + 1] = 1.0
                rhs[row] = -uslip[i]

            # forces
            row_n = m + ms
            row_t = m + ms + 1
            for b, j in enumerate(idx):
                M[row_n, b] = w[j]
                if stick[j]:
                    M[row_t, m + stick_local[int(j)]] = w[j]
                else:
                    M[row_t, b] += w[j] * (slip_sign[j] * mu)
            rhs[row_n] = Fn
            rhs[row_t] = Ft

            try:
                sol = np.linalg.solve(M, rhs)
            except np.linalg.LinAlgError:
                sol = np.linalg.lstsq(M, rhs, rcond=None)[0]

            pa = sol[:m]
            qs = sol[m : m + ms]
            c_n = float(sol[m + ms])
            c_t = float(sol[m + ms + 1])

            p_new = np.zeros(n)
            q_new = np.zeros(n)
            p_new[idx] = pa
            q_new[stick_idx] = qs
            for i in slip_idx:
                q_new[i] = slip_sign[i] * mu * p_new[i]

            changed = False
            neg = pa < -1e-14
            if np.any(neg):
                active[idx[neg]] = False
                stick[idx[neg]] = False
                changed = True

            p = p_new
            q = q_new
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
                    q[i] = (1.0 if q[i] >= 0 else -1.0) * lim
                    changed = True

            # slip -> stick if traction below limit and incremental consistency
            for i in np.flatnonzero(active & ~stick):
                lim = mu * max(p[i], 0.0)
                # tentative stick residual using current uslip
                # keep slip if |q| at limit is required by force; only stick if clearly inside
                if abs(q[i]) < 0.98 * lim:
                    # check if stick would be feasible next — conservative stay in slip
                    pass

            if not changed:
                break

        uslip = -A * (K @ q) - c_t

    return {
        "p": p,
        "q": q,
        "w": w,
        "active": active,
        "stick": stick,
        "K": K,
        "c_n": c_n,
        "c_t": c_t,
    }


def diagnostics(meta: dict, x: np.ndarray, sol: dict) -> dict:
    p = sol["p"]
    q = sol["q"]
    w = sol["w"]
    mu = float(meta["mu"])
    pmax = float(np.max(p)) if p.size else 0.0
    contact = p > 1e-8 * max(pmax, 1e-30)
    contact_length = float(np.sum(w[contact]))
    fn_total = float(np.sum(p * w))
    ft_total = float(np.sum(q * w))
    if np.any(contact):
        p_rms = float(np.sqrt(np.sum((p[contact] ** 2) * w[contact]) / max(contact_length, EPS)))
    else:
        p_rms = 0.0
    if np.any(contact) and abs(fn_total) > 0.0:
        x_p_centroid = float(np.sum(x[contact] * p[contact] * w[contact]) / fn_total)
    else:
        x_p_centroid = 0.0
    if abs(fn_total) > 0.0 and mu > 0.0:
        mobilization = float(abs(ft_total) / (mu * abs(fn_total)))
    else:
        mobilization = 0.0

    return {
        "schema_version": SCHEMA_VERSION,
        "p_peak": float(pmax),
        "p_rms": p_rms,
        "fn_total": fn_total,
        "ft_total": ft_total,
        "contact_length": contact_length,
        "x_p_centroid": x_p_centroid,
        "mobilization": mobilization,
    }


def analyze_pack(pack_dir: Path) -> dict:
    meta, x, g0, fn, ft = _load_pack(pack_dir)
    E = float(meta["E"])
    nu = float(meta["nu"])
    mu = float(meta["mu"])
    rref = float(meta.get("rref", 1.0))
    A = _kernel_coeff(E, nu)
    sol = solve_history(x, g0, A, mu, fn, ft, rref=rref)
    return diagnostics(meta, x, sol)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Contact pack analyzer")
    parser.add_argument("--pack", required=True, help="Absolute path to pack directory")
    parser.add_argument("--out", required=True, help="Absolute path to output JSON")
    args = parser.parse_args(argv)

    pack_dir = Path(args.pack)
    out_path = Path(args.out)
    if not pack_dir.is_dir():
        print(f"pack not found: {pack_dir}", file=sys.stderr)
        return 2

    result = analyze_pack(pack_dir)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
