"""Admissible local blocks and a few explicit parameterizations.

Input order of each block is [e_p, e_q]. Output order is [i_d_ref_pu, i_q_ref_pu].
Ak is not required to be Hurwitz: a zero eigenvalue is the integrator.
"""

from __future__ import annotations

import numpy as np


def block(A, B, C, D, ident: str) -> dict:
    return {
        "id": ident,
        "A": np.array(A, dtype=float),
        "B": np.array(B, dtype=float),
        "C": np.array(C, dtype=float),
        "D": np.array(D, dtype=float),
    }


def zero_blocks() -> list[dict]:
    z = np.zeros((2, 2))
    return [block(z, z, z, z, "converter_1"), block(z, z, z, z, "converter_2")]


def diagonal_pi(kp_p, ki_p, kp_q, ki_q, kp_p2=None, ki_p2=None, kp_q2=None, ki_q2=None) -> list[dict]:
    if kp_p2 is None:
        kp_p2, ki_p2, kp_q2, ki_q2 = kp_p, ki_p, kp_q, ki_q

    def one(kpp, kip, kpq, kiq, ident):
        return block(
            np.zeros((2, 2)),
            np.eye(2),
            np.diag([kip, kiq]),
            np.diag([kpp, kpq]),
            ident,
        )

    return [
        one(kp_p, ki_p, kp_q, ki_q, "converter_1"),
        one(kp_p2, ki_p2, kp_q2, ki_q2, "converter_2"),
    ]


def static_diagonal(k1: float, k2: float | None = None) -> list[dict]:
    """Static gain with unused stable controller poles at -100 rad/s.

    A literal Ak = 0 realization of zero feedthrough puts eigenvalues on the
    imaginary axis. Those are a different, inadmissible realization of the
    same input-output law and are reported separately.
    """
    if k2 is None:
        k2 = k1
    z = np.zeros((2, 2))
    eye = np.eye(2)
    return [
        block(-100.0 * eye, z, z, k1 * eye, "converter_1"),
        block(-100.0 * eye, z, z, k2 * eye, "converter_2"),
    ]


def pi_full_tied(theta: np.ndarray) -> list[dict]:
    """Ak = 0, Bk = I, and C, D are full 2-by-2. Eight parameters, both blocks tied."""
    C = theta[:4].reshape(2, 2)
    D = theta[4:8].reshape(2, 2)
    z = np.zeros((2, 2))
    b = block(z, np.eye(2), C, D, "converter_1")
    c = block(z, np.eye(2), C.copy(), D.copy(), "converter_2")
    return [b, c]


def general_blocks(theta: np.ndarray, tied: bool = False) -> list[dict]:
    """16 parameters per block: A, B, C, D each 2-by-2. Untied uses 32."""
    need = 16 if tied else 32
    if theta.size != need:
        raise ValueError(f"expected {need} parameters")

    def one(vec, ident):
        return block(vec[0:4].reshape(2, 2), vec[4:8].reshape(2, 2), vec[8:12].reshape(2, 2), vec[12:16].reshape(2, 2), ident)

    if tied:
        return [one(theta, "converter_1"), one(theta, "converter_2")]
    return [one(theta[:16], "converter_1"), one(theta[16:], "converter_2")]
