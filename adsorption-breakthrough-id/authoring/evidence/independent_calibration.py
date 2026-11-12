#!/usr/bin/env python3
"""Independent checks for hidden-shift oracle margins and timeout headroom."""

from __future__ import annotations

import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "tests" / "hidden" / "manifest.json"
SOLVER = ROOT / "solution" / "column_solver.py"


def main() -> None:
    import subprocess
    import sys

    manifest = json.loads(MANIFEST.read_text())
    timeout_s = float(manifest["per_run_timeout_s"])
    rows = []
    for entry in manifest["hidden_cases"]:
        case_id = entry["case_id"]
        case_path = ROOT / "tests" / "hidden" / case_id / "case.json"
        out_path = Path("/tmp") / f"calib_{case_id}.json"
        t0 = time.monotonic()
        proc = subprocess.run(
            [sys.executable, str(SOLVER), "--case", str(case_path), "--out", str(out_path)],
            capture_output=True,
            text=True,
            env={
                **dict(**{k: "1" for k in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )})
            },
        )
        elapsed = time.monotonic() - t0
        assert proc.returncode == 0, proc.stderr
        payload = json.loads(out_path.read_text())
        tol = entry["tolerances"]
        ok = True
        for key, rel in (
            ("qm", tol["qm_rel"]),
            ("K", tol["K_rel"]),
            ("k_LDF", tol["k_LDF_rel"]),
            ("D_ax", tol["D_ax_rel"]),
        ):
            err = abs(payload[key] - entry["parameters"][key]) / abs(entry["parameters"][key])
            ok = ok and err <= rel
        Lerr = abs(payload["L_min_m"] - entry["L_min_m"]) / abs(entry["L_min_m"])
        ok = ok and Lerr <= tol["L_min_rel"] and payload["L_min_m"] >= tol["L_min_undersize_floor"] * entry["L_min_m"]
        rows.append(
            {
                "case_id": case_id,
                "elapsed_s": round(elapsed, 2),
                "timeout_s": timeout_s,
                "headroom_s": round(timeout_s - elapsed, 2),
                "pass": ok,
                "L_min_rel_err": round(Lerr, 4),
                "shift": entry.get("shift"),
            }
        )
    report = {
        "oracle_all_pass": all(r["pass"] for r in rows),
        "min_headroom_s": min(r["headroom_s"] for r in rows),
        "rows": rows,
        "hidden_shift_builder": "authoring/provenance/generate_all_cases.py",
    }
    out = Path(__file__).with_name("independent_calibration_report.json")
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
