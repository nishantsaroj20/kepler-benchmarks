#!/usr/bin/env python3
"""Approach-controlled plane-strain frictional contact diagnostics."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

SCHEMA_VERSION = 1
EPS = 1e-15
CONTACT_FRAC = 1e-8


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
    # columns: step, delta_n, delta_t
    dn = loads[:, 1].astype(float)
    dt = loads[:, 2].astype(float)
    return meta, x, g0, dn, dt


def _edges_from_nodes(x: np.ndarray) -> np.ndarray:
    n = len(x)
    edges = np.empty(n + 1, dtype=float)
    if n == 1:
        edges[0], edges[1] = x[0] - 0.5, x[0] + 0.5
        return edges
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - (edges[1] - x[0])
    edges[-1] = x[-1] + (x[-1] - edges[-2])
    return edges


def influence_matrix(x: np.ndarray, rref: float) -> tuple[np.ndarray, np.ndarray]:
    n = len(x)
    edges = _edges_from_nodes(x)
    w = edges[1:] - edges[:-1]
    K = np.zeros((n, n), dtype=float)
    log_r = float(np.log(rref))
    for i in range(n):
        xi = float(x[i])
        for j in range(n):
            L, R = float(edges[j]), float(edges[j + 1])

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
    dn_hist: np.ndarray,
    dt_hist: np.ndarray,
    rref: float,
    max_inner: int = 160,
) -> dict:
    """Prescribed rigid approach δ_n and shift δ_t; forces are unknowns."""
    K, w = influence_matrix(x, rref=rref)
    n = len(x)
    p = np.zeros(n)
    q = np.zeros(n)
    gamma = np.zeros(n)  # accumulated tangential slip (plastic)
    active = np.zeros(n, dtype=bool)
    stick = np.zeros(n, dtype=bool)

    for delta_n, delta_t in zip(dn_hist, dt_hist):
        delta_n = float(delta_n)
        delta_t = float(delta_t)
        # seed contact where geometric interference is largest
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
                s = np.sign(q[sl])
                # tendency from elastic predictor vs delta_t
                pred = A * (K @ q) + gamma
                s = np.sign(pred[sl] - delta_t)
                s[s == 0] = 1.0 if delta_t >= 0 else -1.0
                slip_sign[sl] = s

            # unknowns: p on active, q on stick
            nv = m + ms
            M = np.zeros((nv, nv))
            rhs = np.zeros(nv)
            st_map = {int(j): k for k, j in enumerate(st)}

            # Normal: A K p = g0 - delta_n  on contact
            for a, i in enumerate(idx):
                for b, j in enumerate(idx):
                    M[a, b] = A * K[i, j]
                rhs[a] = g0[i] - delta_n

            # Stick: A K q + gamma = delta_t
            for ka, i in enumerate(st):
                row = m + ka
                for b, j in enumerate(idx):
                    if stick[j]:
                        M[row, m + st_map[int(j)]] += A * K[i, j]
                    else:
                        M[row, b] += A * K[i, j] * (slip_sign[j] * mu)
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
            neg = pa < -1e-14
            if np.any(neg):
                active[idx[neg]] = False
                stick[idx[neg]] = False
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

            # try return-to-stick when clearly inside cone
            u_t = A * (K @ q) + gamma
            for i in np.flatnonzero(active & ~stick):
                lim = mu * max(p[i], 0.0)
                if lim <= 0:
                    continue
                # if tangential compatibility almost holds at |q|<lim, restore stick
                if abs(q[i]) < 0.95 * lim and abs(u_t[i] - delta_t) < 1e-4 * max(1.0, abs(delta_t)):
                    stick[i] = True
                    changed = True

            if not changed:
                break

        # update plastic slip so stick residual uses accumulated state
        # On slip nodes, gamma tracks so that A K q + gamma ≈ delta_t with q on bound
        u_el = A * (K @ q)
        for i in range(n):
            if not active[i]:
                gamma[i] = 0.0
            elif stick[i]:
                # keep gamma; compatibility already A K q + gamma = delta_t
                pass
            else:
                gamma[i] = delta_t - u_el[i]

    return {
        "p": p,
        "q": q,
        "w": w,
        "active": active,
        "stick": stick,
        "gamma": gamma,
        "K": K,
    }


def diagnostics(meta: dict, x: np.ndarray, sol: dict) -> dict:
    p = sol["p"]
    q = sol["q"]
    w = sol["w"]
    gamma = sol["gamma"]
    mu = float(meta["mu"])
    pmax = float(np.max(p)) if p.size else 0.0
    contact = p > CONTACT_FRAC * max(pmax, 1e-30)
    slip = contact & (np.abs(q) >= 0.99 * mu * np.maximum(p, 0.0) - 1e-30)
    # prefer solver stick flags when available, refined by cone
    stick = contact & (~slip)
    if "stick" in sol:
        stick = contact & sol["stick"] & (~slip)
        slip = contact & (~stick)

    contact_length = float(np.sum(w[contact]))
    slip_length = float(np.sum(w[slip]))
    stick_fraction = float(np.sum(w[stick]) / contact_length) if contact_length > 0 else 0.0
    fn_total = float(np.sum(p * w))
    ft_total = float(np.sum(q * w))
    q_peak = float(np.max(np.abs(q))) if q.size else 0.0
    cum_slip_l1 = float(np.sum(w * np.abs(gamma)))
    if np.any(contact) and float(np.sum(np.abs(q[contact]) * w[contact])) > 0:
        x_q_centroid = float(
            np.sum(x[contact] * np.abs(q[contact]) * w[contact])
            / np.sum(np.abs(q[contact]) * w[contact])
        )
    else:
        x_q_centroid = 0.0

    return {
        "schema_version": SCHEMA_VERSION,
        "p_peak": float(pmax),
        "fn_total": fn_total,
        "ft_total": ft_total,
        "contact_length": contact_length,
        "q_peak": q_peak,
        "stick_fraction": stick_fraction,
        "slip_length": slip_length,
        "x_q_centroid": x_q_centroid,
        "cum_slip_l1": cum_slip_l1,
    }


def analyze_pack(pack_dir: Path) -> dict:
    meta, x, g0, dn, dt = _load_pack(pack_dir)
    E = float(meta["E"])
    nu = float(meta["nu"])
    mu = float(meta["mu"])
    rref = float(meta.get("rref", 1.0))
    A = _kernel_coeff(E, nu)
    sol = solve_history(x, g0, A, mu, dn, dt, rref=rref)
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
    # parent must already exist in verifier sandbox; oracle may create
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
