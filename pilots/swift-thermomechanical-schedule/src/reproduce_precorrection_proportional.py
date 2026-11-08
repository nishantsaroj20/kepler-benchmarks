#!/usr/bin/env python3
"""Reproduce the one-variable no-hold proportional family on the PRE-CORRECTION solver.

This file is the provenance command for results/precorrection_immutable/proportional_family.json.
It must be run before the acceptance-rule fix. Durations are (0.35, 0.375, 0.275)*tau with h=0.
Label: one-variable baseline, not a 4-D optimized design.
"""
from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

from coupled_current import quick_card, make_solver
from problem_spec import GAUGE_WINDOW_MM, F_CAP_N, R_FLOOR_MM, SCENARIOS, T_NECK_CAP_K, TAU_TARGET_S, Schedule

OUT = Path(__file__).resolve().parents[1] / "results" / "precorrection_immutable" / "proportional_family.json"

TAUS = (4.0, 5.0, 6.0, 7.0)


def sched(tau):
    return Schedule(0.35 * tau, 0.375 * tau, 0.275 * tau, 0.0)


def margins(row):
    g_lo, g_hi = GAUGE_WINDOW_MM
    return {
        "T_margin_K": T_NECK_CAP_K - row["T_neck_max_K"],
        "tau_margin_s": TAU_TARGET_S - row["tau_s"],
        "r_margin_mm": row["r_min_mm"] - R_FLOOR_MM,
        "F_margin_N": F_CAP_N - row["force_max_N"],
        "gauge_margin_low_mm": row["gauge_mm"] - g_lo,
        "gauge_margin_high_mm": g_hi - row["gauge_mm"],
    }


def main():
    card = quick_card()
    rows = []
    t_all = time.perf_counter()
    for tau in TAUS:
        sch = sched(tau)
        for sc in SCENARIOS:
            t0 = time.perf_counter()
            sol = make_solver(card, sc["dR"], sc["T_init"], 4, 8, k_cond=45.0)
            try:
                out = sol.run(sch, nsub=8, mode="coupled", max_stagger=2)
                err = None
            except Exception as exc:
                out = None
                err = str(exc)
            wall = time.perf_counter() - t0
            if err:
                rec = {"tau_s": tau, "scenario": sc["id"], "error": err, "wall_s": wall, "pass_all": False,
                       "t1": sch.t1, "t2": sch.t2, "t3": sch.t3, "h": 0.0}
                print(f"FAIL tau={tau} {sc['id']} {err}", flush=True)
            else:
                e, f = out["extrema"], out["final"]
                hist = out["history"]
                flags = {
                    "T": e["T_neck_K"] <= T_NECK_CAP_K,
                    "tau": sch.tau <= TAU_TARGET_S,
                    "r": e["r_min_mm"] >= R_FLOOR_MM,
                    "F": e["force_N"] <= F_CAP_N,
                    "g": GAUGE_WINDOW_MM[0] <= f["gauge_mm"] <= GAUGE_WINDOW_MM[1],
                }
                rec = {
                    "label": "one_variable_no_hold_proportional_NOT_4D_optimized",
                    "solver": "precorrection_coupled_current",
                    "mesh": [4, 8], "nsub": 8, "max_stagger": 2, "k_cond": 45.0, "mode": "coupled",
                    "tau_s": sch.tau, "t1": sch.t1, "t2": sch.t2, "t3": sch.t3, "h": sch.h,
                    "scenario": sc["id"], "T_init": sc["T_init"], "dR": sc["dR"],
                    "force_max_N": e["force_N"], "force_end_N": f["force_N"],
                    "T_neck_max_K": e["T_neck_K"], "T_neck_end_K": f["T_neck_K"],
                    "T_shoulder_max_K": e["T_shoulder_K"], "T_shoulder_end_K": f["T_shoulder_K"],
                    "r_min_mm": e["r_min_mm"], "gauge_mm": f["gauge_mm"],
                    "ep_max": e["ep_max"], "ep_neck": e.get("ep_neck"),
                    "n_steps": len(hist), "n_newton_retries": len(out["failures"]),
                    "failures": out["failures"],
                    "heat": out["heat"],
                    "flags": flags, "pass_all": all(flags.values()),
                    "wall_s": wall,
                    "history": hist,
                }
                rec.update(margins(rec))
                print(
                    f"tau={tau:.0f} {sc['id']:16} T={e['T_neck_K']:.2f} r={e['r_min_mm']:.4f} "
                    f"g={f['gauge_mm']:.4f} F={e['force_N']:.1f} retries={len(out['failures'])} "
                    f"pass={rec['pass_all']} {wall:.1f}s",
                    flush=True,
                )
            rows.append(rec)
    payload = {
        "command": "cd src && OMP_NUM_THREADS=1 python3 reproduce_precorrection_proportional.py",
        "python": sys.version,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "note": "Run on the solver BEFORE the acceptance-rule correction. M2/M3 in the old sweep were manually selected; this family is the one-variable proportional baseline with h=0.",
        "wall_s": time.perf_counter() - t_all,
        "n_pass": sum(1 for r in rows if r.get("pass_all")),
        "n_rows": len(rows),
        "rows": rows,
    }
    OUT.write_text(json.dumps(payload))
    print(f"wrote {OUT} pass {payload['n_pass']}/{payload['n_rows']} in {payload['wall_s']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
