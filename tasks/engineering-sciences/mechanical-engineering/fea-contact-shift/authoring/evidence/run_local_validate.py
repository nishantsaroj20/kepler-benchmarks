#!/usr/bin/env python3
"""Harbor-free local validation: oracle-like analyze vs truth; nop/baselines fail."""

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

from analyze import analyze_pack, main as analyze_main  # noqa: E402
from hidden_family import materialize_heldout  # noqa: E402
from truth_solver import truth_for_pack  # noqa: E402

REL = {
    "p_peak": 0.12,
    "p_rms": 0.12,
    "fn_total": 0.02,
    "ft_total": 0.08,
    "contact_length": 0.12,
}
ABS = {"mobilization": 0.05, "x_p_centroid": 5e-4}


def compare(pred, truth, name):
    for k, tol in REL.items():
        err = abs(pred[k] - truth[k]) / max(abs(truth[k]), 1e-30)
        if err > tol:
            raise AssertionError(f"{name} {k} rel={err} > {tol}")
    for k, tol in ABS.items():
        if abs(pred[k] - truth[k]) > tol:
            raise AssertionError(f"{name} {k} abs={abs(pred[k]-truth[k])} > {tol}")


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
        # Simulate oracle artifact
        prog = td_path / "analyze.py"
        shutil.copy2(ROOT / "solution/analyze.py", prog)
        packs = materialize_heldout(prog, td_path / "packs", n_packs=8)
        for pack in packs:
            out = td_path / f"{pack.name}.json"
            t0 = time.time()
            rc = analyze_main(["--pack", str(pack), "--out", str(out)])
            dt = time.time() - t0
            assert rc == 0
            pred = json.loads(out.read_text())
            truth = truth_for_pack(pack)
            compare(pred, truth, pack.name)
            print(f"ORACLE_HELD_PASS {pack.name} wall_s={dt:.4f}")

        # nop: missing program
        nop_ok = not (td_path / "missing.py").exists()
        print(f"NOP_MISSING_PROGRAM {nop_ok}")

        # baseline zero fails
        truth = truth_for_pack(packs[0])
        zeros = {
            "schema_version": 1,
            "p_peak": 0.0,
            "p_rms": 0.0,
            "fn_total": 0.0,
            "ft_total": 0.0,
            "contact_length": 0.0,
            "x_p_centroid": 0.0,
            "mobilization": 0.0,
        }
        try:
            compare(zeros, truth, "zeros")
            print("BASELINE_ZERO_FAIL_UNEXPECTED")
            return 2
        except AssertionError:
            print("BASELINE_ZERO_FAIL_AS_EXPECTED")

        # copy-input style nonsense
        bogus = dict(zeros)
        bogus.update({"fn_total": 1.0, "p_peak": 1.0, "contact_length": 1.0})
        try:
            compare(bogus, truth, "bogus")
            print("BASELINE_CONSTANT_FAIL_UNEXPECTED")
            return 3
        except AssertionError:
            print("BASELINE_CONSTANT_FAIL_AS_EXPECTED")

    print(f"ORACLE_MAX_PRACTICE_S {max(times):.4f}")
    print("LOCAL_VALIDATE_REWARD 1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
