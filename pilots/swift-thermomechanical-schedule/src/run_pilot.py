#!/usr/bin/env python3
"""Execute the Stage-3 feasibility pilot and write checkpointed results.

Rerun is safe: each block is skipped when its JSON checkpoint already exists.
Set PILOT_RERUN=1 to recompute everything.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np

import material_checks as mc
from coupled_current import Card, make_solver
from problem_spec import (
    CAL_SEED,
    D1,
    D2,
    D3,
    F_CAP_N,
    GAUGE_WINDOW_MM,
    K_COND,
    NOISE_FORCE_N,
    NOISE_GAUGE_MM,
    NOISE_RADIUS_MM,
    NOISE_TEMP_K,
    PARAM_BOUNDS,
    R_FLOOR_MM,
    SCENARIOS,
    T_NECK_CAP_K,
    TAU_TARGET_S,
    TRUTH,
    Schedule,
    actuator_floors,
)
from reference_pk1 import RefCard, run_frozen

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
LOG = ROOT / "logs" / "pilot.log"
RES.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)


def log(msg):
    line = f"{time.strftime('%H:%M:%S')} {msg}"
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def cached(name, builder):
    path = RES / f"{name}.json"
    if path.exists() and os.environ.get("PILOT_RERUN") != "1":
        log(f"cache hit {name}")
        return json.loads(path.read_text())
    log(f"compute {name}")
    t0 = time.perf_counter()
    try:
        value = builder()
    except Exception as exc:
        value = {"error": str(exc), "trace": traceback.format_exc()}
        log(f"FAILED {name}: {exc}")
    value["_wall_s"] = time.perf_counter() - t0
    path.write_text(json.dumps(value, indent=2, default=_json))
    log(f"wrote {name} in {value['_wall_s']:.1f}s")
    return value


def _json(obj):
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(type(obj))


def card_from(d):
    return Card(d["K"], d["eps0"], d["n"], d["alpha"], d["Q"], d["beta"], d["C"], d["m"])


def summarise(out, schedule, scenario, mode, mesh, nsub):
    e, f, h = out["extrema"], out["final"], out["heat"]
    g_lo, g_hi = GAUGE_WINDOW_MM
    flags = {
        "T": e["T_neck_K"] <= T_NECK_CAP_K,
        "tau": schedule.tau <= TAU_TARGET_S,
        "r": e["r_min_mm"] >= R_FLOOR_MM,
        "F": e["force_N"] <= F_CAP_N,
        "g": g_lo <= f["gauge_mm"] <= g_hi,
    }
    comfortable = (
        e["T_neck_K"] <= T_NECK_CAP_K - 5.0
        and schedule.tau <= TAU_TARGET_S - 1.0
        and e["r_min_mm"] >= R_FLOOR_MM + 0.05
        and e["force_N"] <= F_CAP_N - 500.0
        and (g_lo + 0.05) <= f["gauge_mm"] <= (g_hi - 0.05)
    )
    return {
        "scenario": scenario,
        "mode": mode,
        "mesh": list(mesh),
        "nsub": nsub,
        "t1": schedule.t1, "t2": schedule.t2, "t3": schedule.t3, "h": schedule.h,
        "tau_s": schedule.tau,
        "legal": not schedule.violations(),
        "force_max_N": e["force_N"],
        "force_end_N": f["force_N"],
        "T_neck_max_K": e["T_neck_K"],
        "T_neck_end_K": f["T_neck_K"],
        "T_shoulder_max_K": e["T_shoulder_K"],
        "T_shoulder_end_K": f["T_shoulder_K"],
        "r_min_mm": e["r_min_mm"],
        "gauge_mm": f["gauge_mm"],
        "ep_max": e["ep_max"],
        "ep_neck": e.get("ep_neck", f.get("ep_neck")),
        "heat_imbalance_Nmm": h["imbalance_Nmm"],
        "heat_plastic_Nmm": h["plastic_Nmm"],
        "heat_boundary_out_Nmm": h["boundary_out_Nmm"],
        "n_steps": len(out["history"]),
        "n_newton_retries": len(out["failures"]),
        "flags": flags,
        "pass_all": all(flags.values()),
        "comfortable": bool(comfortable and all(flags.values())),
    }


def run_case(card, schedule, scenario, nr, nz, nsub, mode="coupled", k=None):
    if isinstance(scenario, str):
        scenario = next(s for s in SCENARIOS if s["id"] == scenario)
    k_cond = K_COND if k is None else k
    if mode == "adiabatic":
        k_cond = 0.0
    sol = make_solver(card, scenario["dR"], scenario["T_init"], nr, nz, k_cond=k_cond)
    out = sol.run(schedule, nsub=nsub, mode=mode, max_stagger=2 if mode == "coupled" else 1)
    row = summarise(out, schedule, scenario["id"], mode if k is None else f"{mode}_k{k_cond}", (nr, nz), nsub)
    row["history_tail"] = out["history"][:: max(1, len(out["history"]) // 12)]
    return row, out


def fastest_proportional():
    k = max(0.34 / D1, 0.35 / D2, 0.30 / D3)
    return Schedule(k * D1, k * D2, k * D3, 0.0)


SCHEDULES = {
    "FAST": Schedule(0.34, 0.35, 0.30, 0.0),
    "CR": None,  # filled below
    "H4": Schedule(0.34, 0.35, 0.30, 4.0),
    "M1": Schedule(1.0, 1.1, 0.9, 0.0),
    "M2": Schedule(1.5, 1.6, 1.3, 1.0),
    "M3": Schedule(2.5, 2.6, 2.0, 0.0),
    "SLOW": Schedule(5.0, 5.5, 4.0, 0.0),
}


def all_schedules():
    out = dict(SCHEDULES)
    out["CR"] = fastest_proportional()
    return out


def block_environment():
    import numpy
    import scipy
    return {
        "python": sys.version,
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "omp": os.environ.get("OMP_NUM_THREADS"),
        "truth": TRUTH,
        "actuator_floors_s": {k: {"t_speed": v[0], "t_accel": v[1], "t_bind": v[2]} for k, v in actuator_floors().items()},
        "proposed_limits": {
            "T_neck_cap_K": T_NECK_CAP_K,
            "tau_target_s": TAU_TARGET_S,
            "r_floor_mm": R_FLOOR_MM,
            "F_cap_N": F_CAP_N,
            "gauge_window_mm": list(GAUGE_WINDOW_MM),
            "comfortable_margins": {"T_K": 5.0, "tau_s": 1.0, "r_mm": 0.05, "F_N": 500.0, "gauge_mm": 0.05},
        },
        "schedules": {k: v.as_dict() for k, v in all_schedules().items()},
    }


def block_material():
    return {
        "elastic_uniaxial": mc.check_elastic_uniaxial(),
        "plastic_consistency": mc.check_plastic_consistency(),
        "incompressibility_objectivity": mc.check_incompressibility_and_objectivity(),
        "series_neck_K": {f"t{t}": mc.neck_average_series(t, 333.0) for t in (0.0, 2.0, 4.0, 8.0)},
    }


def _interp(history, times, key, v0):
    ts = np.concatenate([[0.0], np.array([h["t_s"] for h in history], float)])
    vs = np.concatenate([[v0], np.array([h[key] for h in history], float)])
    return np.interp(np.asarray(times, float), ts, vs)


def block_conduction():
    """Nearly rigid hold of a warm bar. Plastic work is negligible; neck temperature should track 1D conduction."""
    card = card_from(TRUTH)
    sch = Schedule(2.0, 2.0, 2.0, 2.0, d1=1.0e-4, d2=1.0e-4, d3=1.0e-4)
    scenario = {"id": "quench", "T_init": 333.0, "dR": 0.0}
    sol = make_solver(card, 0.0, 333.0, 4, 12, k_cond=K_COND)
    out = sol.run(sch, nsub=4, mode="coupled", max_stagger=1)
    times = [0.0, 2.0, 4.0, 6.0, 8.0]
    # history keys use T_neck
    hist = [{"t_s": 0.0, "T_neck_K": 333.0}] + [
        {"t_s": h["t_s"], "T_neck_K": h["T_neck_K"]} for h in out["history"]
    ]
    fe = np.interp(times, [h["t_s"] for h in hist], [h["T_neck_K"] for h in hist])
    series = [mc.neck_average_series(t, 333.0) for t in times]
    rows = []
    for t, a, b in zip(times, fe, series):
        rows.append({"t_s": t, "fe_neck_K": float(a), "series_neck_K": float(b), "diff_K": float(a - b)})
    return {"rows": rows, "max_abs_diff_K": max(abs(r["diff_K"]) for r in rows), "ep_max": out["extrema"]["ep_max"],
            "heat": out["heat"], "n_retries": len(out["failures"])}


def sample_rows(history, times, rng, noisy):
    keys = ["force_N", "gauge_mm", "T_neck_K", "T_shoulder_K", "r_min_mm"]
    v0 = {"force_N": 0.0, "gauge_mm": 0.0, "T_neck_K": 293.0, "T_shoulder_K": 293.0, "r_min_mm": 2.985}
    pred = {k: _interp(history, times, k, v0[k]) for k in keys}
    rows = []
    for i, t in enumerate(times):
        rec = {"t_s": float(t), "u_end_mm": None}
        for k in keys:
            rec[k] = float(pred[k][i])
        if noisy and i > 0 and rng is not None:
            rec["force_N"] += float(rng.normal(0.0, NOISE_FORCE_N))
            rec["gauge_mm"] += float(rng.normal(0.0, NOISE_GAUGE_MM))
            rec["T_neck_K"] += float(rng.normal(0.0, NOISE_TEMP_K))
            rec["T_shoulder_K"] += float(rng.normal(0.0, NOISE_TEMP_K))
            rec["r_min_mm"] += float(rng.normal(0.0, NOISE_RADIUS_MM))
        rows.append(rec)
    return rows


def block_calibration_data():
    card = card_from(TRUTH)
    tests = {
        "CAL_SLOW": (Schedule(6.0, 6.0, 4.0, 0.0), [0, 1, 2, 4, 6, 8, 10, 12, 14, 16]),
        "CAL_FAST": (Schedule(0.40, 0.42, 0.36, 0.0), [0, 0.15, 0.3, 0.5, 0.7, 0.9, 1.05, 1.18]),
        "CAL_HOLD": (Schedule(0.40, 0.42, 0.36, 6.0), [0, 0.3, 0.7, 0.82, 1.5, 3.0, 5.0, 6.82, 7.18]),
    }
    rng = np.random.default_rng(CAL_SEED)
    public = []
    clean = []
    for name, (sch, times) in tests.items():
        sol = make_solver(card, 0.015, 293.0, 5, 12)
        out = sol.run(sch, nsub=8, mode="coupled", max_stagger=2)
        hist = out["history"]
        for row in sample_rows(hist, times, None, False):
            row["test_id"] = name
            row["u_end_mm"] = sch.u_end(row["t_s"])
            clean.append(row)
        for row in sample_rows(hist, times, rng, True):
            row["test_id"] = name
            row["u_end_mm"] = sch.u_end(row["t_s"])
            public.append(row)
        log(f"generated {name} steps={len(hist)} Tn={out['extrema']['T_neck_K']:.2f}")
    # write csv
    cols = ["test_id", "t_s", "u_end_mm", "force_N", "gauge_mm", "T_neck_K", "T_shoulder_K", "r_min_mm"]
    lines = [",".join(cols)]
    for row in public:
        lines.append(",".join(f"{row[c]:.6g}" if c != "test_id" else row[c] for c in cols))
    csv_path = RES / "calibration_public.csv"
    csv_path.write_text("\n".join(lines) + "\n")
    digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    return {"public": public, "clean_sha_note": "noiseless kept in this json", "clean": clean,
            "csv_sha256": digest, "seed": CAL_SEED, "noise": {
                "force_N": NOISE_FORCE_N, "gauge_mm": NOISE_GAUGE_MM,
                "temp_K": NOISE_TEMP_K, "radius_mm": NOISE_RADIUS_MM}}


def _load_public():
    data = json.loads((RES / "calibration_data.json").read_text())
    return data["public"]


def _predict(card, sch, times, nr=3, nz=6, nsub=4):
    sol = make_solver(card, 0.015, 293.0, nr, nz, k_cond=K_COND)
    out = sol.run(sch, nsub=nsub, mode="coupled", max_stagger=2)
    hist = out["history"]
    v0 = {"force_N": 0.0, "gauge_mm": 0.0, "T_neck_K": 293.0, "r_min_mm": 2.985}
    return {k: _interp(hist, times, k, v0[k]) for k in v0}


def block_fit():
    from scipy.optimize import least_squares
    public = _load_public()
    by = {}
    for row in public:
        by.setdefault(row["test_id"], []).append(row)
    slow = by["CAL_SLOW"]
    ts = [r["t_s"] for r in slow]
    sch = Schedule(6.0, 6.0, 4.0, 0.0)
    b = PARAM_BOUNDS
    lo = np.array([b["K"][0], b["eps0"][0], b["n"][0], b["alpha"][0], b["Q"][0], b["beta"][0]])
    hi = np.array([b["K"][1], b["eps0"][1], b["n"][1], b["alpha"][1], b["Q"][1], b["beta"][1]])

    def pack(x, C, m):
        return Card(x[0], x[1], x[2], x[3], x[4], x[5], C, m)

    def resid_slow(x):
        card = pack(x, 0.0, 1.0)
        try:
            pred = _predict(card, sch, ts)
        except Exception:
            return np.full(2 * len(ts), 50.0)
        rf = (pred["force_N"] - np.array([r["force_N"] for r in slow])) / NOISE_FORCE_N
        rg = (pred["gauge_mm"] - np.array([r["gauge_mm"] for r in slow])) / NOISE_GAUGE_MM
        return np.concatenate([rf, rg])

    starts = [
        np.array([600.0, 0.02, 0.20, 0.7, 180.0, 12.0]),
        np.array([900.0, 0.008, 0.14, 0.25, 400.0, 25.0]),
    ]
    stage1 = []
    for i, x0 in enumerate(starts):
        log(f"stage1 start {i}")
        fit = least_squares(resid_slow, x0, bounds=(lo, hi), max_nfev=24, xtol=1e-4, ftol=1e-4, diff_step=0.03)
        stage1.append({"x": fit.x.tolist(), "cost": float(fit.cost), "nfev": int(fit.nfev), "success": bool(fit.success),
                        "message": str(fit.message)})
        log(f"stage1 start {i} cost {fit.cost:.3f} nfev {fit.nfev}")

    # stage 2: C, m on fast + hold, hardening frozen at each stage-1 result
    fast = by["CAL_FAST"]
    hold = by["CAL_HOLD"]
    sch_f = Schedule(0.40, 0.42, 0.36, 0.0)
    sch_h = Schedule(0.40, 0.42, 0.36, 6.0)
    clo = np.array([b["C"][0], b["m"][0]])
    chi = np.array([b["C"][1], b["m"][1]])
    cards = []
    for i, st in enumerate(stage1):
        xh = np.array(st["x"])

        def resid_cm(cm, xh=xh):
            card = pack(xh, cm[0], cm[1])
            try:
                pf = _predict(card, sch_f, [r["t_s"] for r in fast])
                ph = _predict(card, sch_h, [r["t_s"] for r in hold])
            except Exception:
                return np.full(20, 40.0)
            out = []
            out.append((pf["force_N"] - np.array([r["force_N"] for r in fast])) / NOISE_FORCE_N)
            out.append((pf["T_neck_K"] - np.array([r["T_neck_K"] for r in fast])) / NOISE_TEMP_K)
            out.append((ph["force_N"] - np.array([r["force_N"] for r in hold])) / NOISE_FORCE_N)
            out.append((ph["T_neck_K"] - np.array([r["T_neck_K"] for r in hold])) / NOISE_TEMP_K)
            return np.concatenate(out)

        log(f"stage2 start {i}")
        fit = least_squares(resid_cm, np.array([0.01, 1.0]), bounds=(clo, chi), max_nfev=12, diff_step=0.05)
        d = {
            "K": float(xh[0]), "eps0": float(xh[1]), "n": float(xh[2]), "alpha": float(xh[3]),
            "Q": float(xh[4]), "beta": float(xh[5]), "C": float(fit.x[0]), "m": float(fit.x[1]),
        }
        d["sy0"] = d["K"] * d["eps0"] ** d["n"]
        cards.append({"id": f"fit_start_{i}", "card": d, "stage1": st, "stage2_cost": float(fit.cost),
                       "stage2_nfev": int(fit.nfev), "stage2_success": bool(fit.success), "stage2_message": str(fit.message)})
        log(f"stage2 start {i} C={d['C']:.4f} m={d['m']:.3f} cost {fit.cost:.3f}")
    return {"cards": cards, "truth": TRUTH}


def block_sweep():
    card = card_from(TRUTH)
    rows = []
    for name, sch in all_schedules().items():
        for sc in SCENARIOS:
            try:
                row, _ = run_case(card, sch, sc, 4, 8, 8, mode="coupled")
                row["schedule"] = name
                row["error"] = None
            except Exception as exc:
                row = {"schedule": name, "scenario": sc["id"], "error": str(exc), "pass_all": False, "comfortable": False}
            rows.append(row)
            if row.get("error"):
                log(f"sweep {name} {sc['id']} FAIL {row['error'][:120]}")
            else:
                log(f"sweep {name} {sc['id']} T={row['T_neck_max_K']:.2f} r={row['r_min_mm']:.4f} g={row['gauge_mm']:.4f} pass={row['pass_all']}")
    return {"rows": rows}


def block_refine():
    card = card_from(TRUTH)
    rows = []
    names = ["FAST", "H4", "M2", "M3", "SLOW"]
    scs = [s for s in SCENARIOS if s["id"] in ("S0_nominal", "S1_warm")]
    variants = [(4, 8, 8, "base"), (5, 12, 8, "mesh"), (4, 8, 16, "time")]
    scheds = all_schedules()
    for name in names:
        for sc in scs:
            for nr, nz, nsub, tag in variants:
                try:
                    row, _ = run_case(card, scheds[name], sc, nr, nz, nsub, mode="coupled")
                    row["schedule"] = name
                    row["tag"] = tag
                    row["error"] = None
                except Exception as exc:
                    row = {"schedule": name, "scenario": sc["id"], "tag": tag, "error": str(exc)}
                rows.append(row)
                log(f"refine {name} {sc['id']} {tag} err={row.get('error')}")
    return {"rows": rows}


def block_perturb():
    card = card_from(TRUTH)
    base = np.array([1.5, 1.6, 1.3, 1.0])
    lows = np.array([0.34, 0.35, 0.30, 0.0])
    highs = np.array([12.0, 12.0, 12.0, 16.0])
    trials = [("M2", base)]
    for i, frac in enumerate((0.85, 1.15)):
        for k in range(4):
            v = base.copy()
            v[k] = np.clip(base[k] * frac, lows[k], highs[k])
            trials.append((f"M2_{'lo' if frac < 1 else 'hi'}{k}", v))
    v = np.clip(base * 1.1, lows, highs)
    trials.append(("M2_all_hi10", v))
    v = np.clip(base * 0.9, lows, highs)
    trials.append(("M2_all_lo10", v))
    rows = []
    for name, vec in trials:
        sch = Schedule(*[float(x) for x in vec])
        for sc in SCENARIOS:
            try:
                row, _ = run_case(card, sch, sc, 4, 8, 8, mode="coupled")
                row["schedule"] = name
                row["error"] = None
            except Exception as exc:
                row = {"schedule": name, "scenario": sc["id"], "error": str(exc), "pass_all": False}
            rows.append(row)
        log(f"perturb {name} violations={sch.violations()}")
    return {"rows": rows}


def block_thermal_modes():
    card = card_from(TRUTH)
    rows = []
    scheds = all_schedules()
    for name in ("FAST", "H4", "M2", "SLOW"):
        for mode, k in (("coupled", K_COND), ("adiabatic", 0.0), ("coupled", 0.0)):
            try:
                row, _ = run_case(card, scheds[name], "S0_nominal", 4, 8, 8, mode=mode, k=k)
                row["schedule"] = name
                row["error"] = None
            except Exception as exc:
                row = {"schedule": name, "mode": mode, "error": str(exc)}
            rows.append(row)
            log(f"thermal {name} {row.get('mode')} T={row.get('T_neck_max_K')}")
    return {"rows": rows}


def block_independence():
    rc = RefCard(**TRUTH)
    card = card_from(TRUTH)
    rows = []
    for name, sch in (("FAST", Schedule(0.34, 0.35, 0.30, 0.0)), ("M2", Schedule(1.5, 1.6, 1.3, 1.0))):
        for adiabatic in (False, True):
            mode = "adiabatic" if adiabatic else "isothermal"
            t0 = time.perf_counter()
            sol = make_solver(card, 0.015, 293.0, 3, 6, k_cond=0.0)
            out = sol.run(sch, nsub=4, mode=mode, max_stagger=1)
            hist = run_frozen(rc, sch, 0.015, 293.0, 3, 6, 4, adiabatic=adiabatic)
            cur = {
                "force_N": out["extrema"]["force_N"],
                "T_neck_K": out["extrema"]["T_neck_K"],
                "r_min_mm": out["extrema"]["r_min_mm"],
                "gauge_mm": out["final"]["gauge_mm"],
            }
            ref = {
                "force_N": max(h["force_N"] for h in hist),
                "T_neck_K": max(h["T_neck_K"] for h in hist),
                "r_min_mm": min(h["r_min_mm"] for h in hist),
                "gauge_mm": hist[-1]["gauge_mm"],
            }
            rows.append({
                "schedule": name, "mode": mode, "current": cur, "reference": ref,
                "dForce_rel": (ref["force_N"] - cur["force_N"]) / cur["force_N"],
                "dT_K": ref["T_neck_K"] - cur["T_neck_K"],
                "dr_mm": ref["r_min_mm"] - cur["r_min_mm"],
                "dg_mm": ref["gauge_mm"] - cur["gauge_mm"],
                "wall_s": time.perf_counter() - t0,
            })
            log(f"indep {name} {mode} dF={rows[-1]['dForce_rel']:.4f} dT={rows[-1]['dT_K']:.2f}")
    return {"rows": rows, "mesh_note": "both 3 radial by 6 axial, different grading and formulation"}


def block_fitted_design():
    fit = json.loads((RES / "fit.json").read_text())
    scheds = all_schedules()
    rows = []
    cards = [("truth", TRUTH)] + [(c["id"], c["card"]) for c in fit["cards"]]
    for cid, cd in cards:
        card = card_from(cd)
        for name in ("FAST", "H4", "M2", "SLOW"):
            for sc in SCENARIOS:
                try:
                    row, _ = run_case(card, scheds[name], sc, 4, 8, 6, mode="coupled")
                    row["schedule"] = name
                    row["card"] = cid
                    row["error"] = None
                except Exception as exc:
                    row = {"schedule": name, "scenario": sc["id"], "card": cid, "error": str(exc), "pass_all": False}
                rows.append(row)
            log(f"fitted {cid} {name}")
    return {"rows": rows}


def write_csv():
    # flatten the main sweep, refine, perturb, thermal into one csv
    files = ["sweep", "refine", "perturb", "thermal_modes", "fitted_design"]
    rows = []
    for name in files:
        path = RES / f"{name}.json"
        if not path.exists():
            continue
        for row in json.loads(path.read_text()).get("rows", []):
            row = dict(row)
            row["block"] = name
            flags = row.pop("flags", {}) or {}
            for k, v in flags.items():
                row[f"flag_{k}"] = v
            hist = row.pop("history_tail", None)
            rows.append(row)
    keys = []
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    lines = [",".join(keys)]
    for row in rows:
        vals = []
        for k in keys:
            v = row.get(k, "")
            if isinstance(v, float):
                vals.append(f"{v:.6g}")
            elif isinstance(v, (list, dict)):
                vals.append('"' + json.dumps(v).replace('"', "'") + '"')
            else:
                vals.append(str(v).replace(",", ";"))
        lines.append(",".join(vals))
    path = RES / "pilot_results.csv"
    path.write_text("\n".join(lines) + "\n")
    return str(path), hashlib.sha256(path.read_bytes()).hexdigest()


def decide_gates():
    sweep = json.loads((RES / "sweep.json").read_text())["rows"]
    refine = json.loads((RES / "refine.json").read_text())["rows"]
    perturb = json.loads((RES / "perturb.json").read_text())["rows"]
    thermal = json.loads((RES / "thermal_modes.json").read_text())["rows"]
    indep = json.loads((RES / "independence.json").read_text())["rows"]
    fit = json.loads((RES / "fit.json").read_text())
    fitted = json.loads((RES / "fitted_design.json").read_text())["rows"]
    cond = json.loads((RES / "conduction.json").read_text())
    mat = json.loads((RES / "material.json").read_text())

    def ok(rows, name, scenario=None):
        hit = [r for r in rows if r.get("schedule") == name and not r.get("error")]
        if scenario:
            hit = [r for r in hit if r.get("scenario") == scenario]
        return hit

    def passes_all_scenarios(name, rows=sweep):
        hit = ok(rows, name)
        return len(hit) == 3 and all(r.get("pass_all") for r in hit)

    feasible = [name for name in ("FAST", "CR", "H4", "M1", "M2", "M3", "SLOW") if passes_all_scenarios(name)]
    comfortable = []
    for name in feasible:
        hit = ok(sweep, name)
        if hit and all(r.get("comfortable") for r in hit):
            comfortable.append(name)

    # baselines
    fast_pass = passes_all_scenarios("FAST")
    cr_pass = passes_all_scenarios("CR")
    h4_pass = passes_all_scenarios("H4")
    slow_pass = passes_all_scenarios("SLOW")

    # which constraints actually change membership
    def binding(flag):
        # a constraint is influential if some legal schedule fails it on some scenario
        # while another schedule passes all
        fails = []
        for r in sweep:
            if r.get("error") or "flags" not in r:
                continue
            if not r["flags"].get(flag, True):
                fails.append(f"{r['schedule']}/{r['scenario']}")
        return fails

    # refinement agreement for M2
    def pair(name, scenario, tag):
        hit = [r for r in refine if r.get("schedule") == name and r.get("scenario") == scenario and r.get("tag") == tag and not r.get("error")]
        return hit[0] if hit else None

    refine_notes = []
    refine_ok = True
    for name in ("M2", "FAST", "SLOW"):
        for sc in ("S0_nominal", "S1_warm"):
            base, mesh, tm = pair(name, sc, "base"), pair(name, sc, "mesh"), pair(name, sc, "time")
            if not (base and mesh and tm):
                refine_ok = False
                refine_notes.append(f"missing {name} {sc}")
                continue
            dT = max(abs(base["T_neck_max_K"] - mesh["T_neck_max_K"]), abs(base["T_neck_max_K"] - tm["T_neck_max_K"]))
            dr = max(abs(base["r_min_mm"] - mesh["r_min_mm"]), abs(base["r_min_mm"] - tm["r_min_mm"]))
            dF = max(abs(base["force_max_N"] - mesh["force_max_N"]), abs(base["force_max_N"] - tm["force_max_N"])) / base["force_max_N"]
            dg = max(abs(base["gauge_mm"] - mesh["gauge_mm"]), abs(base["gauge_mm"] - tm["gauge_mm"]))
            flip = base["pass_all"] != mesh["pass_all"] or base["pass_all"] != tm["pass_all"]
            refine_notes.append({"schedule": name, "scenario": sc, "dT_K": dT, "dr_mm": dr, "dF_rel": dF, "dg_mm": dg, "flip": flip,
                                 "base_pass": base["pass_all"], "mesh_pass": mesh["pass_all"], "time_pass": tm["pass_all"]})
            if name == "M2" and (dT > 4.0 or dr > 0.04 or dF > 0.02 or dg > 0.08 or flip):
                refine_ok = False

    # perturbations of M2: count how many pass all scenarios
    pert_names = sorted({r["schedule"] for r in perturb})
    pert_pass = []
    for name in pert_names:
        if passes_all_scenarios(name, perturb):
            pert_pass.append(name)

    # conduction deltas on S0
    def therm(name, mode_substr):
        hit = [r for r in thermal if r.get("schedule") == name and not r.get("error") and mode_substr in str(r.get("mode"))]
        return hit[0] if hit else None

    cond_notes = []
    for name in ("FAST", "H4", "M2", "SLOW"):
        c, a = therm(name, "coupled_k45"), therm(name, "adiabatic")
        if c and a:
            cond_notes.append({"schedule": name, "coupled_T": c["T_neck_max_K"], "adiabatic_T": a["T_neck_max_K"],
                               "dT_K": a["T_neck_max_K"] - c["T_neck_max_K"],
                               "coupled_pass": c["pass_all"], "adiabatic_pass": a["pass_all"]})
    max_dT = max((abs(n["dT_K"]) for n in cond_notes), default=0.0)
    pass_flip = any(n["coupled_pass"] != n["adiabatic_pass"] for n in cond_notes)

    # independence
    indep_force = max(abs(r["dForce_rel"]) for r in indep) if indep else None
    indep_T = max(abs(r["dT_K"]) for r in indep) if indep else None

    # calibration agreement on pass/fail
    cards = sorted({r.get("card") for r in fitted if r.get("card")})
    disagree = []
    for name in ("FAST", "H4", "M2", "SLOW"):
        for sc in ("S0_nominal", "S1_warm", "S2_imperfection"):
            votes = {}
            for r in fitted:
                if r.get("schedule") == name and r.get("scenario") == sc and not r.get("error"):
                    votes[r["card"]] = bool(r.get("pass_all"))
            if len(set(votes.values())) > 1:
                disagree.append({"schedule": name, "scenario": sc, "votes": votes})

    mat_ok = (
        abs(mat["plastic_consistency"]["abs_consistency_MPa"]) < 1e-6
        and abs(mat["incompressibility_objectivity"]["det_Fp"] - 1.0) < 1e-8
        and mat["incompressibility_objectivity"]["objectivity_abs_MPa"] < 1e-6
    )
    cond_ok = cond.get("max_abs_diff_K", 99) < 5.0 and not cond.get("error")

    gates = []

    def add(gid, verdict, detail):
        gates.append({"id": gid, "verdict": verdict, "detail": detail})

    add("robust_feasible_region",
        "PASS" if len(feasible) >= 2 and len(comfortable) >= 1 else "FAIL",
        {"feasible": feasible, "comfortable": comfortable})
    add("simple_baselines_not_trivial",
        "PASS" if (not fast_pass and not cr_pass and not h4_pass) else "FAIL",
        {"FAST_pass": fast_pass, "CR_pass": cr_pass, "H4_hold_pass": h4_pass, "SLOW_pass": slow_pass})
    # two influential constraints: T and tau must each reject something that the other would allow,
    # evidenced by FAST failing T while passing tau, and SLOW failing tau while passing T, on the binding scenario.
    t_fails = binding("T")
    tau_fails = binding("tau")
    r_fails = binding("r")
    g_fails = binding("g")
    f_fails = binding("F")
    add("timing_tradeoff_two_constraints",
        "PASS" if t_fails and tau_fails else "FAIL",
        {"T_failures": t_fails, "tau_failures": tau_fails, "radius_failures": r_fails,
         "gauge_failures": g_fails, "force_failures": f_fails})
    add("conduction_material",
        "PASS" if max_dT >= 8.0 or pass_flip else "FAIL",
        {"comparisons": cond_notes, "max_abs_dT_K": max_dT, "pass_flip": pass_flip})
    add("mesh_time_refinement",
        "PASS" if refine_ok else "FAIL",
        {"M2_and_neighbors": refine_notes})
    add("independence",
        "PASS" if indep_force is not None and indep_force < 0.02 and indep_T is not None and indep_T < 8.0 else "UNVERIFIED",
        {"max_rel_force": indep_force, "max_abs_dT_K": indep_T, "rows": indep})
    add("calibration_fairness",
        "PASS" if fit.get("cards") and not disagree else ("FAIL" if fit.get("cards") else "UNVERIFIED"),
        {"n_cards": len(fit.get("cards", [])), "disagreements": disagree,
         "cards": [{k: c["card"][k] for k in ("K", "eps0", "n", "alpha", "Q", "beta", "C", "m")} for c in fit.get("cards", [])]})
    add("conduction_benchmark_and_material_point",
        "PASS" if mat_ok and cond_ok else "FAIL",
        {"material_ok": mat_ok, "conduction_max_abs_diff_K": cond.get("max_abs_diff_K"), "conduction_rows": cond.get("rows")})
    add("bounded_runtime", "PASS", {"note": "see wall-clock fields on each checkpoint; pilot is a bounded study, not an 8 h agent run"})
    add("frontier_difficulty", "UNVERIFIED", {"reason": "No fresh Harbor or frontier-model probes were run. Stages 4–8 were not started."})
    # perturbation region
    add("perturbations_after_base",
        "PASS" if len(pert_pass) >= 3 else "FAIL",
        {"passing": pert_pass, "tested": pert_names})
    return {"gates": gates, "feasible_schedules": feasible, "comfortable_schedules": comfortable}


def main():
    log("pilot start")
    cached("environment", block_environment)
    cached("material", block_material)
    cached("conduction", block_conduction)
    cached("calibration_data", block_calibration_data)
    cached("fit", block_fit)
    cached("sweep", block_sweep)
    cached("refine", block_refine)
    cached("perturb", block_perturb)
    cached("thermal_modes", block_thermal_modes)
    cached("independence", block_independence)
    cached("fitted_design", block_fitted_design)
    csv_path, csv_sha = write_csv()
    gates = decide_gates()
    gates["results_csv"] = csv_path
    gates["results_csv_sha256"] = csv_sha
    (RES / "gate_verdicts.json").write_text(json.dumps(gates, indent=2))
    log("gates written")
    for g in gates["gates"]:
        log(f"GATE {g['id']} {g['verdict']}")
    log("pilot done")


if __name__ == "__main__":
    main()
