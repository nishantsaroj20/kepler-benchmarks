#!/usr/bin/env python3
"""Harbor-free local validation of oracle vs Gauss-truth; defective probes fail."""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "solution"))
sys.path.insert(0, str(ROOT / "tests"))

from analyze import analyze_pack, solve_history, _load_pack, _kernel_coeff, diagnostics  # noqa: E402
from hidden_family import materialize_heldout  # noqa: E402
from truth_solver import truth_for_pack  # noqa: E402
import numpy as np  # noqa: E402

REL = {
    "p_peak": 0.05,
    "fn_total": 0.02,
    "ft_total": 0.02,
    "contact_length": 0.10,
    "q_peak": 0.10,
    "slip_length": 0.15,
    "cum_slip_l1": 0.30,
}
ABS = {"stick_fraction": 0.08, "x_q_centroid": 1e-4}


def compare(pred, truth, name):
    for k, tol in REL.items():
        p, t = float(pred[k]), float(truth[k])
        if abs(t) < 1e-30 and abs(p) < 1e-30:
            continue
        err = abs(p - t) / max(abs(t), 1e-30)
        if err > tol:
            raise AssertionError(f"{name} {k} rel={err} > {tol}")
    for k, tol in ABS.items():
        if abs(float(pred[k]) - float(truth[k])) > tol:
            raise AssertionError(f"{name} {k} abs diff too large")


def main() -> int:
    practice = ROOT / "environment/data/practice"
    times = []
    for pack in sorted(practice.glob("practice_*")):
        t0 = time.time()
        pred = analyze_pack(pack)
        times.append(time.time() - t0)
        truth = truth_for_pack(pack)
        compare(pred, truth, pack.name)
        print(f"ORACLE_PRACTICE_PASS {pack.name} wall_s={times[-1]:.4f}")

    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        prog = td_path / "analyze.py"
        shutil.copy2(ROOT / "solution/analyze.py", prog)
        packs = materialize_heldout(prog, td_path / "packs", n_packs=8)
        for pack in packs:
            t0 = time.time()
            pred = analyze_pack(pack)
            dt = time.time() - t0
            truth = truth_for_pack(pack)
            compare(pred, truth, pack.name)
            print(f"ORACLE_HELD_PASS {pack.name} wall_s={dt:.4f}")

        # frictionless probe on a contacting pack
        contact_pack = next(p for p in packs if truth_for_pack(p)["contact_length"] > 0)
        truth = truth_for_pack(contact_pack)
        meta, x, g0, dn, dtv = _load_pack(contact_pack)
        A = _kernel_coeff(meta["E"], meta["nu"])
        sol = solve_history(x, g0, A, meta["mu"], dn, np.zeros_like(dtv), meta["rref"])
        probe = diagnostics(meta, x, sol)
        try:
            compare(probe, truth, "zero_dt_probe")
            print("FRICTIONLESS_PROBE_UNEXPECTED_PASS")
            return 2
        except AssertionError:
            print("FRICTIONLESS_PROBE_FAIL_AS_EXPECTED")

        zeros = {k: (1 if k == "schema_version" else 0.0) for k in truth}
        try:
            compare(zeros, truth, "zeros")
            print("BASELINE_ZERO_FAIL_UNEXPECTED")
            return 3
        except AssertionError:
            print("BASELINE_ZERO_FAIL_AS_EXPECTED")

    print(f"ORACLE_MAX_PRACTICE_S {max(times):.4f}")
    print("LOCAL_VALIDATE_REWARD 1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
