#!/usr/bin/env python3
"""One fresh command for the corrected pilot. Does not read the old results cache.

Outputs go to results/corrected/. The pre-correction proportional run stays in
results/precorrection_immutable/ and is not overwritten.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from coupled_current import Card, make_solver, mesh_bar, quick_card
from problem_spec import (
    CAL_SEED,
    F_CAP_N,
    GAUGE_WINDOW_MM,
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
)
from reference_pk1 import RefCard, mesh_uniform, run_conducting

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "corrected"
HIST = OUT / "histories"
LOG = ROOT / "logs" / "corrected_pilot.log"
OUT.mkdir(parents=True, exist_ok=True)
HIST.mkdir(parents=True, exist_ok=True)


def log(msg):
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def props(tau):
    return Schedule(0.35 * tau, 0.375 * tau, 0.275 * tau, 0.0)


def m2():
    return Schedule(1.5, 1.6, 1.3, 1.0)


def flags_of(force, tn, radius, gauge, tau):
    return {
        "T": tn <= T_NECK_CAP_K,
        "tau": tau <= TAU_TARGET_S,
        "r": radius >= R_FLOOR_MM,
        "F": force <= F_CAP_N,
        "g": GAUGE_WINDOW_MM[0] <= gauge <= GAUGE_WINDOW_MM[1],
    }


def summarize(name, sch, sc, out, mesh, nsub, k, tol, wall, tag):
    e, f = out["extrema"], out["final"]
    after = [h for h in out["history"] if h["t_s"] > 0.0]
    tn_after = max(h["T_neck_K"] for h in after)
    fl = flags_of(e["force_N"], e["T_neck_K"], e["r_min_mm"], f["gauge_mm"], sch.tau)
    return {
        "tag": tag, "name": name, "scenario": sc["id"], "T_init": sc["T_init"], "dR": sc["dR"],
        "t1": sch.t1, "t2": sch.t2, "t3": sch.t3, "h": sch.h, "tau_s": sch.tau,
        "mesh_r": mesh[0], "mesh_z": mesh[1], "nsub": nsub, "k": k, "stagger_tol_K": tol,
        "T_neck_max_K": e["T_neck_K"], "T_neck_t_s": e["T_neck_t_s"],
        "T_neck_max_after_t0_K": tn_after,
        "T_shoulder_max_K": e["T_shoulder_K"],
        "force_max_N": e["force_N"], "force_t_s": e["force_t_s"],
        "r_min_mm": e["r_min_mm"], "r_min_t_s": e["r_min_t_s"],
        "gauge_mm": f["gauge_mm"], "ep_max": e["ep_max"],
        "T_margin_K": T_NECK_CAP_K - e["T_neck_K"],
        "tau_margin_s": TAU_TARGET_S - sch.tau,
        "r_margin_mm": e["r_min_mm"] - R_FLOOR_MM,
        "F_margin_N": F_CAP_N - e["force_N"],
        "g_margin_low_mm": f["gauge_mm"] - GAUGE_WINDOW_MM[0],
        "g_margin_high_mm": GAUGE_WINDOW_MM[1] - f["gauge_mm"],
        "n_steps": len(out["history"]), "n_rejected": len(out["failures"]),
        "max_stagger_delta_K": max(h["stagger_delta_K"] for h in out["history"]),
        "heat_imbalance_Nmm": out["heat"]["imbalance_Nmm"],
        "pass_all": all(fl.values()), "flags": fl, "wall_s": wall, "error": "",
    }


def save_history(key, out):
    rows = out["history"]
    cols = ["t_s", "u_end_mm", "force_N", "gauge_mm", "r_min_mm", "T_neck_K", "T_shoulder_K",
            "ep_max", "ep_neck", "residual_N", "stagger_delta_K", "dt_s", "newton_iter", "stagger_passes"]
    arr = {c: np.array([h.get(c, np.nan) for h in rows], float) for c in cols}
    np.savez_compressed(HIST / f"{key}.npz", **arr)


def run_primary(name, sch, sc, nr, nz, nsub, k, tol, tag, card=None):
    t0 = time.perf_counter()
    sol = make_solver(quick_card() if card is None else card, sc["dR"], sc["T_init"], nr, nz, k_cond=k)
    try:
        out = sol.run(sch, nsub=nsub, mode="coupled", max_stagger=4, stagger_tol=tol)
        err = None
    except Exception as exc:
        err = str(exc)
        out = None
    wall = time.perf_counter() - t0
    if err:
        log(f"FAIL {tag} {name} {sc['id']} {err}")
        return {"tag": tag, "name": name, "scenario": sc["id"], "mesh_r": nr, "mesh_z": nz, "nsub": nsub,
                "k": k, "stagger_tol_K": tol, "pass_all": False, "error": err, "wall_s": wall, "tau_s": sch.tau}
    key = f"{tag}_{name}_{sc['id']}_m{nr}x{nz}_n{nsub}_k{k:g}_tol{tol:g}"
    save_history(key, out)
    row = summarize(name, sch, sc, out, (nr, nz), nsub, k, tol, wall, tag)
    row["history_file"] = key + ".npz"
    log(f"{tag} {name} {sc['id']} m{nr}x{nz} n{nsub} k={k:g} T={row['T_neck_max_K']:.2f} "
        f"g={row['gauge_mm']:.4f} r={row['r_min_mm']:.4f} pass={row['pass_all']} rej={row['n_rejected']} {wall:.1f}s")
    return row


def primary_studies():
    rows = []
    family = [(f"P{int(tau)}", props(tau)) for tau in (4, 5, 6, 7)]
    family.append(("M2_manual", m2()))
    # Base decision mesh. P* is the one-variable family. M2 was selected by hand in the old sweep.
    for name, sch in family:
        for sc in SCENARIOS:
            rows.append(run_primary(name, sch, sc, 4, 8, 8, 45.0, 0.25, "base"))
    # Mesh only, time only, for P4, P7, M2, all scenarios.
    focus = [("P4", props(4)), ("P7", props(7)), ("M2_manual", m2())]
    for name, sch in focus:
        for sc in SCENARIOS:
            for nr, nz in ((3, 6), (5, 12)):
                rows.append(run_primary(name, sch, sc, nr, nz, 8, 45.0, 0.25, "mesh"))
            for nsub in (4, 16):
                rows.append(run_primary(name, sch, sc, 4, 8, nsub, 45.0, 0.25, "time"))
    # Combined finer grid on the temperature-critical warm member and M2.
    for name, sch in (("P4", props(4)), ("P7", props(7)), ("M2_manual", m2())):
        for sc in SCENARIOS:
            if not (name == "P4" or sc["id"] == "S1_warm"):
                continue
            rows.append(run_primary(name, sch, sc, 5, 12, 16, 45.0, 0.25, "combined"))
    # Coupling tolerance, same grid, warm critical cases.
    warm = SCENARIOS[1]
    for name, sch in (("P4", props(4)), ("M2_manual", m2())):
        for tol in (0.5, 0.05):
            rows.append(run_primary(name, sch, warm, 4, 8, 8, 45.0, tol, "coupling_tol"))
    # Same integrator, vary only k. End Dirichlet stays T_grip.
    for name, sch in (("P4", props(4)), ("M2_manual", m2()), ("P7", props(7))):
        for sc in SCENARIOS:
            rows.append(run_primary(name, sch, sc, 4, 8, 8, 0.0, 0.25, "k0"))
    return rows


def independent_studies():
    rows = []
    card = RefCard(**TRUTH)
    jobs = []
    for name, sch in (("P4", props(4)), ("M2_manual", m2())):
        for sc in SCENARIOS:
            for k in (45.0, 0.0):
                jobs.append((name, sch, sc, 3, 6, 4, k, "indep_coarse"))
    # One refined pair on the warm P4 case, both conductivities, own mesh.
    jobs.append(("P4", props(4), SCENARIOS[1], 4, 8, 4, 45.0, "indep_finer"))
    jobs.append(("P4", props(4), SCENARIOS[1], 4, 8, 4, 0.0, "indep_finer"))
    for name, sch, sc, nr, nz, nsub, k, tag in jobs:
        t0 = time.perf_counter()
        try:
            hist = run_conducting(card, sch, sc["dR"], sc["T_init"], nr, nz, nsub, k=k, stagger_tol=0.5, max_stagger=3)
            err = ""
        except Exception as exc:
            log(f"FAIL {tag} {name} {sc['id']} k={k} {exc}")
            rows.append({"tag": tag, "name": name, "scenario": sc["id"], "k": k, "mesh_r": nr, "mesh_z": nz,
                         "nsub": nsub, "pass_all": False, "error": str(exc), "wall_s": time.perf_counter() - t0})
            continue
        after = [h for h in hist if h["t_s"] > 0]
        tn = max(h["T_neck_K"] for h in hist)
        tn_after = max(h["T_neck_K"] for h in after)
        force = max(h["force_N"] for h in hist)
        radius = min(h["r_min_mm"] for h in hist)
        gauge = hist[-1]["gauge_mm"]
        fl = flags_of(force, tn, radius, gauge, sch.tau)
        key = f"{tag}_{name}_{sc['id']}_m{nr}x{nz}_n{nsub}_k{k:g}"
        np.savez_compressed(
            HIST / f"{key}.npz",
            t_s=np.array([h["t_s"] for h in hist]),
            force_N=np.array([h["force_N"] for h in hist]),
            gauge_mm=np.array([h["gauge_mm"] for h in hist]),
            r_min_mm=np.array([h["r_min_mm"] for h in hist]),
            T_neck_K=np.array([h["T_neck_K"] for h in hist]),
        )
        row = {
            "tag": tag, "name": name, "scenario": sc["id"], "tau_s": sch.tau, "k": k,
            "mesh_r": nr, "mesh_z": nz, "nsub": nsub, "stagger_tol_K": 0.5,
            "T_neck_max_K": tn, "T_neck_max_after_t0_K": tn_after,
            "force_max_N": force, "r_min_mm": radius, "gauge_mm": gauge,
            "T_margin_K": T_NECK_CAP_K - tn, "r_margin_mm": radius - R_FLOOR_MM,
            "g_margin_low_mm": gauge - GAUGE_WINDOW_MM[0],
            "n_steps": len(hist), "pass_all": all(fl.values()), "error": "",
            "wall_s": time.perf_counter() - t0, "history_file": key + ".npz",
            "solver": "reference_pk1_conducting",
        }
        rows.append(row)
        log(f"{tag} {name} {sc['id']} k={k:g} T={tn:.2f} g={gauge:.4f} r={radius:.4f} F={force:.0f} pass={row['pass_all']} {row['wall_s']:.1f}s")
    # Same-node diagnosis: primary mesh, both conductivities, warm P4, reference solver only.
    xy, conn, _ = mesh_bar(3.0, 0.015, 10.0, 4, 8)
    for k in (45.0, 0.0):
        t0 = time.perf_counter()
        try:
            hist = run_conducting(card, props(4), 0.015, 333.0, 4, 8, 4, k=k, stagger_tol=0.5, max_stagger=3, xy=xy, conn=conn)
            tn = max(h["T_neck_K"] for h in hist)
            row = {"tag": "same_nodes_as_primary", "name": "P4", "scenario": "S1_warm", "k": k,
                   "mesh_r": 4, "mesh_z": 8, "nsub": 4, "T_neck_max_K": tn,
                   "gauge_mm": hist[-1]["gauge_mm"], "r_min_mm": min(h["r_min_mm"] for h in hist),
                   "force_max_N": max(h["force_N"] for h in hist), "n_steps": len(hist),
                   "wall_s": time.perf_counter() - t0, "error": ""}
        except Exception as exc:
            row = {"tag": "same_nodes_as_primary", "name": "P4", "scenario": "S1_warm", "k": k, "error": str(exc),
                   "wall_s": time.perf_counter() - t0}
        rows.append(row)
        log(f"same-nodes k={k:g} {row}")
    return rows


def _interp_hist(hist, times, key, v0):
    ts = np.concatenate([[0.0], np.array([h["t_s"] for h in hist], float)])
    vs = np.concatenate([[v0], np.array([h[key] for h in hist], float)])
    return np.interp(times, ts, vs)


def calibration():
    """Regenerate observations with the corrected solver and fit all channels."""
    from scipy.optimize import least_squares
    card = quick_card()
    tests = {
        "CAL_SLOW": (Schedule(6.0, 6.0, 4.0, 0.0), [0, 1, 2, 4, 6, 8, 10, 12, 14, 16]),
        "CAL_FAST": (Schedule(0.40, 0.42, 0.36, 0.0), [0, 0.15, 0.3, 0.5, 0.7, 0.9, 1.05, 1.18]),
        "CAL_HOLD": (Schedule(0.40, 0.42, 0.36, 6.0), [0, 0.3, 0.7, 0.82, 1.5, 3.0, 5.0, 6.82, 7.18]),
    }
    levels = {"gen_5x12": (5, 12, 8), "fitmesh_4x8": (4, 8, 6), "coarse_3x6": (3, 6, 4)}
    clean = {}
    for level, (nr, nz, nsub) in levels.items():
        clean[level] = {}
        for tid, (sch, times) in tests.items():
            sol = make_solver(card, 0.015, 293.0, nr, nz)
            out = sol.run(sch, nsub=nsub, mode="coupled", max_stagger=4, stagger_tol=0.25)
            v0 = {"force_N": 0.0, "gauge_mm": 0.0, "T_neck_K": 293.0, "T_shoulder_K": 293.0, "r_min_mm": 2.985}
            clean[level][tid] = {k: _interp_hist(out["history"], times, k, v0[k]) for k in v0}
            log(f"cal {level} {tid} Tn={out['extrema']['T_neck_K']:.2f}")
    # discretization discrepancy of the fit mesh against the generating mesh
    disc = {}
    for tid in tests:
        for key, noise in (("force_N", NOISE_FORCE_N), ("gauge_mm", NOISE_GAUGE_MM),
                           ("T_neck_K", NOISE_TEMP_K), ("T_shoulder_K", NOISE_TEMP_K),
                           ("r_min_mm", NOISE_RADIUS_MM)):
            d = clean["fitmesh_4x8"][tid][key] - clean["gen_5x12"][tid][key]
            d = d[1:]  # drop the exact t=0 row
            disc[f"{tid}.{key}"] = {"rms": float(np.sqrt(np.mean(d * d))), "max_abs": float(np.max(np.abs(d))), "noise": noise}
    rng = np.random.default_rng(CAL_SEED)
    public_lines = ["test_id,t_s,u_end_mm,force_N,gauge_mm,T_neck_K,T_shoulder_K,r_min_mm"]
    obs = {}
    for tid, (sch, times) in tests.items():
        obs[tid] = {"t": np.array(times, float)}
        base = clean["gen_5x12"][tid]
        noisy = {}
        for key, sig in (("force_N", NOISE_FORCE_N), ("gauge_mm", NOISE_GAUGE_MM),
                         ("T_neck_K", NOISE_TEMP_K), ("T_shoulder_K", NOISE_TEMP_K),
                         ("r_min_mm", NOISE_RADIUS_MM)):
            v = base[key].copy()
            if len(v) > 1:
                v[1:] = v[1:] + rng.normal(0.0, sig, len(v) - 1)
            noisy[key] = v
        obs[tid].update(noisy)
        for i, t in enumerate(times):
            public_lines.append(",".join([
                tid, f"{t:.6g}", f"{sch.u_end(t):.6g}",
                f"{noisy['force_N'][i]:.6g}", f"{noisy['gauge_mm'][i]:.6g}",
                f"{noisy['T_neck_K'][i]:.6g}", f"{noisy['T_shoulder_K'][i]:.6g}",
                f"{noisy['r_min_mm'][i]:.6g}",
            ]))
    csv_path = OUT / "calibration_public.csv"
    csv_path.write_text("\n".join(public_lines) + "\n")

    def sigma(tid, key, noise):
        rms = disc[f"{tid}.{key}"]["rms"]
        return float(np.hypot(noise, rms))

    def predict(card_obj, tid, nr=4, nz=8, nsub=6):
        sch, times = tests[tid]
        sol = make_solver(card_obj, 0.015, 293.0, nr, nz)
        out = sol.run(sch, nsub=nsub, mode="coupled", max_stagger=3, stagger_tol=0.5)
        v0 = {"force_N": 0.0, "gauge_mm": 0.0, "T_neck_K": 293.0, "T_shoulder_K": 293.0, "r_min_mm": 2.985}
        return {k: _interp_hist(out["history"], times, k, v0[k]) for k in v0}

    b = PARAM_BOUNDS
    lo = np.array([b["K"][0], b["eps0"][0], b["n"][0], b["alpha"][0], b["Q"][0], b["beta"][0], b["C"][0], b["m"][0]])
    hi = np.array([b["K"][1], b["eps0"][1], b["n"][1], b["alpha"][1], b["Q"][1], b["beta"][1], b["C"][1], b["m"][1]])
    starts = [
        np.array([680, 0.015, 0.18, 0.5, 220, 14, 0.012, 0.8]),
        np.array([600, 0.02, 0.20, 0.7, 180, 12, 0.01, 1.0]),
        np.array([900, 0.008, 0.15, 0.25, 400, 25, 0.02, 0.6]),
        np.array([500, 0.04, 0.25, 0.9, 80, 8, 0.005, 1.4]),
    ]

    def resid(x):
        card_obj = Card(*[float(v) for v in x])
        pieces = []
        try:
            for tid in tests:
                pred = predict(card_obj, tid)
                for key, noise in (("force_N", NOISE_FORCE_N), ("gauge_mm", NOISE_GAUGE_MM),
                                   ("T_neck_K", NOISE_TEMP_K), ("T_shoulder_K", NOISE_TEMP_K),
                                   ("r_min_mm", NOISE_RADIUS_MM)):
                    scale = sigma(tid, key, noise)
                    pieces.append((pred[key] - obs[tid][key]) / scale)
        except Exception:
            return np.full(135, 20.0)
        out = np.concatenate(pieces)
        if out.size != 135:
            return np.full(135, 20.0)
        return out

    fits = []
    for i, x0 in enumerate(starts):
        wanted = os.environ.get("FIT_STARTS")
        if wanted and str(i) not in {p.strip() for p in wanted.split(",")}:
            continue
        log(f"joint fit start {i}")
        t0 = time.perf_counter()
        fit = least_squares(resid, x0, bounds=(lo, hi), max_nfev=18, diff_step=0.08, xtol=1e-3, ftol=1e-3)
        n = len(fit.fun)
        rms = float(np.sqrt(2.0 * fit.cost / max(n, 1)))
        d = {"K": float(fit.x[0]), "eps0": float(fit.x[1]), "n": float(fit.x[2]), "alpha": float(fit.x[3]),
             "Q": float(fit.x[4]), "beta": float(fit.x[5]), "C": float(fit.x[6]), "m": float(fit.x[7])}
        fits.append({"start": i, "x0": x0.tolist(), "card": d, "cost": float(fit.cost), "rms_standardized": rms, "n_residual": n,
                     "nfev": int(fit.nfev), "success": bool(fit.success), "message": str(fit.message),
                     "wall_s": time.perf_counter() - t0})
        (OUT / "fits_partial.json").write_text(json.dumps(fits, indent=2))
        log(f"fit start {i} rms={rms:.3f} nfev={fit.nfev} K={d['K']:.1f} alpha={d['alpha']:.3f} C={d['C']:.4f} m={d['m']:.3f}")
    # Replay the simple family member closest to the cap, on each card, decision mesh.
    replays = []
    cards = [("truth", TRUTH)] + [(f"fit{f['start']}", f["card"]) for f in fits]
    for cid, cd in cards:
        cobj = Card(cd["K"], cd["eps0"], cd["n"], cd["alpha"], cd["Q"], cd["beta"], cd["C"], cd["m"])
        for sc in SCENARIOS:
            replays.append(run_primary(f"card_{cid}", props(4), sc, 4, 8, 8, 45.0, 0.25, "card_replay", card=cobj))
            replays[-1]["card"] = cid
    return {"discretization_vs_noise": disc, "fits": fits, "replays": replays,
            "calibration_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest()}


def decide(rows, indep, cal):
    def base_family(name):
        hit = [r for r in rows if r.get("tag") == "base" and r.get("name") == name and not r.get("error")]
        return hit

    p_pass = {}
    for tau in (4, 5, 6, 7):
        hit = base_family(f"P{tau}")
        p_pass[str(tau)] = len(hit) == 3 and all(r.get("pass_all") for r in hit)
    family_broad = all(p_pass.values())
    # refinement of P4 warm
    def find(tag, name, scenario, **kw):
        for r in rows:
            if r.get("tag") == tag and r.get("name") == name and r.get("scenario") == scenario and not r.get("error"):
                if all(r.get(k) == v for k, v in kw.items()):
                    return r
        return None

    p4_warm_levels = []
    for tag, extra in (
        ("base", {"mesh_r": 4, "nsub": 8, "k": 45.0}),
        ("mesh", {"mesh_r": 3, "nsub": 8}),
        ("mesh", {"mesh_r": 5, "nsub": 8}),
        ("time", {"nsub": 4}),
        ("time", {"nsub": 16}),
        ("combined", {"mesh_r": 5, "nsub": 16}),
        ("k0", {"k": 0.0}),
    ):
        rec = find(tag, "P4", "S1_warm", **{k: extra[k] for k in extra})
        if rec:
            p4_warm_levels.append({k: rec.get(k) for k in ("tag", "mesh_r", "mesh_z", "nsub", "k", "T_neck_max_K", "T_neck_max_after_t0_K", "gauge_mm", "r_min_mm", "pass_all", "error")})
    if family_broad:
        decision = "FOUR_VARIABLE_SCHEDULE_NOT_DEMONSTRATED"
        text = (
            "The one-variable no-hold family t=(0.35,0.375,0.275)*tau at tau=4,5,6,7 s "
            "passes every public scenario on the corrected base mesh. A user can pick a moderate "
            "duration with h=0. No independent ramp-duration optimization was required, and none was run. "
            "M2 and M3 in the previous report were manual selections. Do not tighten the 8 s cap or the "
            "other limits to push this family out. The four-variable hardening is not demonstrated."
        )
    else:
        decision = "SIMPLE_FAMILY_DOES_NOT_PASS_BROADLY"
        text = (
            "The corrected base mesh does not pass every member of the proportional no-hold family. "
            "See p_pass. This is not a claim that a four-variable search is difficult, and it is not a "
            "fresh model-probe result."
        )
    return {
        "decision": decision,
        "text": text,
        "proportional_base_pass": p_pass,
        "p4_warm_levels": p4_warm_levels,
        "independent_rows": len(indep),
        "fit_rms": [f.get("rms_standardized") for f in cal.get("fits", [])],
    }


def main():
    LOG.write_text("")
    t0 = time.perf_counter()
    log("corrected pilot start")
    log(f"python {sys.version.split()[0]} numpy {np.__version__}")
    import scipy
    log(f"scipy {scipy.__version__} platform {platform.platform()}")
    # Acceptance tests are part of the command.
    import test_acceptance
    test_acceptance.test_restore_bitwise()
    test_acceptance.test_rejected_newton_rolls_back()
    test_acceptance.test_rejected_stagger_does_not_keep_temperature()
    test_acceptance.test_clean_step_records_convergence()
    log("acceptance tests passed")
    rows = primary_studies()
    indep = independent_studies()
    cal = calibration()
    decision = decide(rows, indep, cal)
    # compact summary csv
    keys = []
    for r in rows + indep + cal["replays"]:
        for k in r:
            if k != "flags" and k not in keys:
                keys.append(k)
    with (OUT / "summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows + indep + cal["replays"]:
            w.writerow({k: r.get(k, "") for k in keys})
    payload = {
        "command": "cd pilots/swift-thermomechanical-schedule/src && OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 run_corrected_pilot.py",
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "decision": decision,
        "calibration_discretization": cal["discretization_vs_noise"],
        "fits": cal["fits"],
        "calibration_sha256": cal["calibration_sha256"],
        "n_primary_rows": len(rows),
        "n_independent_rows": len(indep),
        "wall_s": time.perf_counter() - t0,
    }
    (OUT / "gates_and_decision.json").write_text(json.dumps(payload, indent=2, default=float))
    log(f"DECISION {decision['decision']}")
    log(decision["text"])
    log(f"done in {payload['wall_s']:.1f}s")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        log("FATAL\n" + traceback.format_exc())
        sys.exit(1)
