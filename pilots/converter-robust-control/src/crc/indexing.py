"""State layout shared as a naming convention, not as an equation implementation."""

from __future__ import annotations

N_PER = 18


def base(conv: int) -> int:
    return conv * N_PER


# Offsets inside one converter block.
IC, VC, IO = 0, 2, 4
XI, VT, PLL, DELTA = 6, 8, 10, 11
VF, ICF, IOF = 12, 14, 16
