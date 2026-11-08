"""Closed-loop reconstruction for D_yu = 0.

Two constructions are returned: the block formula and a block-equation
interconnection that never calls the formula. They are compared by the caller.
"""

from __future__ import annotations

import numpy as np


def stack_controller(blocks: list[dict]) -> dict:
    """blocks are in converter order, each with A,B,C,D as ndarrays."""
    Ak = _blkdiag([b["A"] for b in blocks])
    Bk = _blkdiag([b["B"] for b in blocks])
    Ck = _blkdiag([b["C"] for b in blocks])
    Dk = _blkdiag([b["D"] for b in blocks])
    return {"A": Ak, "B": Bk, "C": Ck, "D": Dk}


def _blkdiag(mats: list[np.ndarray]) -> np.ndarray:
    rows = sum(m.shape[0] for m in mats)
    cols = sum(m.shape[1] for m in mats)
    out = np.zeros((rows, cols))
    r = c = 0
    for m in mats:
        out[r : r + m.shape[0], c : c + m.shape[1]] = m
        r += m.shape[0]
        c += m.shape[1]
    return out


def closed_loop_formula(plant: dict, K: dict) -> dict:
    A = plant["A"]
    Bw = plant["Bw"]
    Bu = plant["Bu"]
    Cy = plant["Cy"]
    Dyw = plant["Dyw"]
    Cz = plant["Cz"]
    Dzw = plant["Dzw"]
    Dzu = plant["Dzu"]
    Ak, Bk, Ck, Dk = K["A"], K["B"], K["C"], K["D"]
    if plant["Dyu"].shape != (Cy.shape[0], Bu.shape[1]) or np.any(plant["Dyu"] != 0.0):
        raise ValueError("formula is only for an exactly zero Dyu")
    Acl = np.block(
        [
            [A + Bu @ Dk @ Cy, Bu @ Ck],
            [Bk @ Cy, Ak],
        ]
    )
    Bcl = np.vstack([Bw + Bu @ Dk @ Dyw, Bk @ Dyw])
    Ccl = np.hstack([Cz + Dzu @ Dk @ Cy, Dzu @ Ck])
    Dcl = Dzw + Dzu @ Dk @ Dyw
    return {"A": Acl, "B": Bcl, "C": Ccl, "D": Dcl}


def closed_loop_interconnection(plant: dict, K: dict) -> dict:
    """Assemble the series/parallel blocks without using the compact formula."""
    A = np.array(plant["A"], float)
    Bw = np.array(plant["Bw"], float)
    Bu = np.array(plant["Bu"], float)
    Cy = np.array(plant["Cy"], float)
    Dyw = np.array(plant["Dyw"], float)
    Cz = np.array(plant["Cz"], float)
    Dzw = np.array(plant["Dzw"], float)
    Dzu = np.array(plant["Dzu"], float)
    Ak = np.array(K["A"], float)
    Bk = np.array(K["B"], float)
    Ck = np.array(K["C"], float)
    Dk = np.array(K["D"], float)
    n = A.shape[0]
    nc = Ak.shape[0]
    # State order [x, xc]. Input w. Output z.
    # xdot = A x + Bw w + Bu u
    # xcdot = Ak xc + Bk y
    # y = Cy x + Dyw w
    # u = Ck xc + Dk y
    # z = Cz x + Dzw w + Dzu u
    y_from_x = Cy
    u_from_x = Dk @ y_from_x
    u_from_xc = Ck
    u_from_w = Dk @ Dyw
    Acl = np.zeros((n + nc, n + nc))
    Acl[:n, :n] = A + Bu @ u_from_x
    Acl[:n, n:] = Bu @ u_from_xc
    Acl[n:, :n] = Bk @ y_from_x
    Acl[n:, n:] = Ak
    Bcl = np.zeros((n + nc, Bw.shape[1]))
    Bcl[:n, :] = Bw + Bu @ u_from_w
    Bcl[n:, :] = Bk @ Dyw
    z_from_x = Cz + Dzu @ u_from_x
    z_from_xc = Dzu @ u_from_xc
    z_from_w = Dzw + Dzu @ u_from_w
    Ccl = np.hstack([z_from_x, z_from_xc])
    Dcl = z_from_w
    return {"A": Acl, "B": Bcl, "C": Ccl, "D": Dcl}


def formula_matches_interconnection(plant: dict, K: dict, atol: float = 1e-9) -> float:
    a = closed_loop_formula(plant, K)
    b = closed_loop_interconnection(plant, K)
    err = 0.0
    for key in ("A", "B", "C", "D"):
        err = max(err, float(np.max(np.abs(a[key] - b[key]))))
    if err > atol:
        raise AssertionError(f"closed-loop constructions differ by {err}")
    return err
