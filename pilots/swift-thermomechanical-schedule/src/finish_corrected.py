#!/usr/bin/env python3
"""Assemble saved corrected histories, then finish the joint fit. Does not rerun the sweep."""
from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path

import numpy as np

from problem_spec import F_CAP_N, GAUGE_WINDOW_MM, R_FLOOR_MM, T_NECK_CAP_K, TAU_TARGET_S, Schedule
import run_corrected_pilot as rc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "corrected"
HIST = OUT / "histories"


def schedule_of(name):
    if name.startswith("P") and name[1:].isdigit():
        tau = float(name[1:])
        return Schedule(0.35 * tau, 0.375 * tau, 0.275 * tau, 0.0)
    if name == "M2_manual":
        return Schedule(1.5, 1.6, 1.3, 1.0)
    raise KeyError(name)


def row_from_npz(path: Path):
    m = re.match(
        r"(?P<tag>.+?)_(?P<case>P\d+|M2_manual)_(?P<sc>S\d+_\w+)_m(?P<nr>\d+)x(?P<nz>\d+)_n(?P<nsub>\d+)_k(?P<k>[0-9.]+)(?:_tol(?P<tol>[0-9.]+))?$",
        path.stem,
    )
    if not m:
        return None
    z = np.load(path)
    t, tn, force, radius, gauge = z["t_s"], z["T_neck_K"], z["force_N"], z["r_min_mm"], z["gauge_mm"]
    sch = schedule_of(m.group("case"))
    tn_max, g_final, rmin, fmax = float(np.max(tn)), float(gauge[-1]), float(np.min(radius)), float(np.max(force))
    fl = {
        "T": tn_max <= T_NECK_CAP_K, "tau": sch.tau <= TAU_TARGET_S, "r": rmin >= R_FLOOR_MM,
        "F": fmax <= 14000.0, "g": GAUGE_WINDOW_MM[0] <= g_final <= GAUGE_WINDOW_MM[1],
    }
    return {
        "tag": m.group("tag"), "name": m.group("case"), "scenario": m.group("sc"),
        "mesh_r": int(m.group("nr")), "mesh_z": int(m.group("nz")), "nsub": int(m.group("nsub")),
        "k": float(m.group("k")), "stagger_tol_K": float(m.group("tol") or 0.5),
        "tau_s": sch.tau, "T_neck_max_K": tn_max,
        "T_neck_max_after_t0_K": float(np.max(tn[t > 0])) if np.any(t > 0) else tn_max,
        "force_max_N": fmax, "r_min_mm": rmin, "gauge_mm": g_final,
        "T_margin_K": T_NECK_CAP_K - tn_max, "r_margin_mm": rmin - R_FLOOR_MM,
        "g_margin_low_mm": g_final - GAUGE_WINDOW_MM[0], "F_margin_N": F_CAP_N - fmax,
        "tau_margin_s": TAU_TARGET_S - sch.tau, "n_steps": int(t.size),
        "pass_all": all(fl.values()), "history_file": path.name,
    }


def main():
    rows = [r for p in sorted(HIST.glob("*.npz")) if (r := row_from_npz(p))]
    failures = [ln.split(" ", 1)[1] for ln in Path("/tmp/corrected_stdout.txt").read_text().splitlines() if " FAIL " in ln]
    (OUT / "summary_partial.json").write_text(json.dumps({"n": len(rows), "failures": failures}, indent=2))
    print("assembled", len(rows), "logged failures", len(failures), flush=True)
    os.environ.setdefault("FIT_STARTS", "0,2,3")
    cal = rc.calibration()
    payload = {"rows": rows, "logged_failures": failures, "fits": cal["fits"],
               "discretization_vs_noise": cal["discretization_vs_noise"],
               "calibration_sha256": cal["calibration_sha256"], "replays": cal["replays"]}
    (OUT / "assembled_rows.json").write_text(json.dumps(payload, indent=2, default=float))
    keys = []
    for r in rows + cal["replays"]:
        for k in r:
            if k != "flags" and k not in keys:
                keys.append(k)
    with (OUT / "summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows + cal["replays"]:
            w.writerow(r)
    print("fits", [(f["start"], round(f["rms_standardized"], 3), f["card"]) for f in cal["fits"]], flush=True)


if __name__ == "__main__":
    main()
