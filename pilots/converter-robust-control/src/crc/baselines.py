"""Bounded baseline campaign. Thresholds are read from params and are not edited here."""

from __future__ import annotations

import time

import numpy as np
from scipy.optimize import differential_evolution, minimize

from .closed_loop import closed_loop_formula, stack_controller
from .controllers import diagonal_pi, general_blocks, pi_full_tied, static_diagonal, zero_blocks, block
from .evaluate import evaluate_blocks, physical_closed_loop
from .norms import classify_stability, slicot_linf, stability_report
from .params import GAMMA_PASS, STAB_PASS
from .synthesis import full_order_hinf, state_feedback_feasible


def _jsonify(value):
    if isinstance(value, dict):
        return {k: _jsonify(v) for k, v in value.items() if k != "blocks"}
    if isinstance(value, list):
        return [_jsonify(v) for v in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def fast_case(plant: dict, blocks: list[dict]) -> dict:
    K = stack_controller(blocks)
    A_phys = physical_closed_loop(plant, K)
    phys = stability_report(A_phys)
    phys_class = classify_stability(phys, STAB_PASS)
    if phys_class != "pass":
        return {
            "case_id": plant["case_id"],
            "pass_stability": False,
            "stability_class": phys_class,
            "abscissa": phys["abscissa"],
            "uncertainty": phys["uncertainty"],
            "slicot_linf": None,
            "pass_hinf": False,
        }
    cl = closed_loop_formula(plant, K)
    # Hurwitz gate for the weighted loop uses a loose threshold; the margin is physical.
    aug = stability_report(cl["A"])
    aug_class = classify_stability(aug, -1e-3)
    if aug_class != "pass":
        return {
            "case_id": plant["case_id"],
            "pass_stability": False,
            "stability_class": phys_class,
            "augmented_hurwitz": aug_class,
            "abscissa": phys["abscissa"],
            "uncertainty": phys["uncertainty"],
            "slicot_linf": None,
            "pass_hinf": False,
        }
    value = float(slicot_linf(cl["A"], cl["B"], cl["C"], cl["D"])["value"])
    return {
        "case_id": plant["case_id"],
        "pass_stability": True,
        "stability_class": phys_class,
        "abscissa": phys["abscissa"],
        "uncertainty": phys["uncertainty"],
        "slicot_linf": value,
        "pass_hinf": value <= GAMMA_PASS,
    }


def fast_eval(blocks, plants) -> dict:
    rows = [fast_case(p, blocks) for p in plants]
    passed = all(r["pass_stability"] and r["pass_hinf"] for r in rows)
    worst = max((r["slicot_linf"] if r["slicot_linf"] is not None else 1e6) for r in rows)
    worst_a = max(r["abscissa"] for r in rows)
    return {"pass_all": passed, "worst_slicot": worst, "worst_abscissa": worst_a, "per_case": rows}


def _record(family, start_id, seed, theta, result, wall, n_eval, note, search_cost):
    return {
        "family": family,
        "start_id": start_id,
        "seed": seed,
        "theta": None if theta is None else np.asarray(theta, float).tolist(),
        "search_cost": search_cost,
        "wall_s": wall,
        "n_evaluations": n_eval,
        "note": note,
        "fast": _jsonify(result),
    }


def pi_grid(plants) -> list[dict]:
    rows = []
    kps = [0.1, 0.2, 0.4, 0.8, 1.5, 3.0]
    kis = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 40.0]
    n = 0
    t0 = time.perf_counter()
    best = None
    for kp in kps:
        for ki in kis:
            n += 1
            ev = fast_eval(diagonal_pi(kp, ki, kp, ki), plants)
            rec = _record(
                "tied_diagonal_pi_grid",
                n,
                None,
                [kp, ki, kp, ki],
                ev,
                time.perf_counter() - t0,
                n,
                "tied diagonal PI, identical blocks",
                ev["worst_slicot"],
            )
            rows.append(rec)
            if best is None or ev["worst_slicot"] < best[0]:
                best = (ev["worst_slicot"], rec, diagonal_pi(kp, ki, kp, ki))
    rows.append(
        _record(
            "tied_diagonal_pi_grid_best_pointer",
            0,
            None,
            best[1]["theta"],
            best[1]["fast"] if False else fast_eval(best[2], plants),
            time.perf_counter() - t0,
            n,
            "pointer restates the best grid entry",
            best[0],
        )
    )
    return rows, best[2], best[0]


def static_grid(plants) -> list[dict]:
    rows = []
    t0 = time.perf_counter()
    n = 0
    for k in list(np.geomspace(0.05, 8.0, 12)) + [-0.4, -1.0, 0.0]:
        n += 1
        ev = fast_eval(static_diagonal(float(k)), plants)
        rows.append(
            _record(
                "static_diagonal",
                n,
                None,
                [float(k)],
                ev,
                time.perf_counter() - t0,
                n,
                "static tied gain, stable unused poles at -100",
                ev["worst_slicot"],
            )
        )
    return rows


def optimize_family(plants, family, theta0_fn, bounds, seeds, maxiter, tied_general=False, pi_cross=False):
    rows = []
    counter = {"n": 0}

    def cost(theta):
        counter["n"] += 1
        if pi_cross:
            blocks = pi_full_tied(np.asarray(theta, float))
        elif tied_general or len(bounds) == 16:
            blocks = general_blocks(np.asarray(theta, float), tied=True)
        elif len(bounds) == 32:
            blocks = general_blocks(np.asarray(theta, float), tied=False)
        else:
            kp, ki, kq, kq_i = theta
            blocks = diagonal_pi(kp, ki, kq, kq_i)
        ev = fast_eval(blocks, plants)
        # Stable failures stay large so the search prefers Hurwitz controllers.
        if not all(r["pass_stability"] for r in ev["per_case"]):
            return 100.0 + max(0.0, ev["worst_abscissa"])
        return float(ev["worst_slicot"])

    best_blocks = None
    best_cost = 1e9
    for i, seed in enumerate(seeds):
        t0 = time.perf_counter()
        counter["n"] = 0
        if len(bounds) <= 8:
            result = differential_evolution(
                cost,
                bounds,
                seed=seed,
                maxiter=maxiter,
                popsize=6,
                mutation=0.5,
                recombination=0.7,
                workers=1,
                updating="immediate",
                atol=0.01,
                tol=0.01,
                polish=True,
            )
            theta = result.x
            message = result.message
        else:
            rng = np.random.default_rng(seed)
            x0 = np.array([rng.uniform(lo, hi) for lo, hi in bounds])
            result = minimize(
                cost,
                x0,
                method="Powell",
                bounds=bounds,
                options={"maxfev": maxiter, "xtol": 1e-3, "ftol": 1e-3},
            )
            theta = result.x
            message = result.message
        if pi_cross:
            blocks = pi_full_tied(theta)
        elif len(bounds) >= 16:
            blocks = general_blocks(theta, tied=(len(bounds) == 16))
        else:
            blocks = diagonal_pi(*theta)
        ev = fast_eval(blocks, plants)
        rows.append(
            _record(family, i, seed, theta, ev, time.perf_counter() - t0, counter["n"], str(message), ev["worst_slicot"])
        )
        if ev["worst_slicot"] < best_cost:
            best_cost = ev["worst_slicot"]
            best_blocks = blocks
    return rows, best_blocks


def nominal_only(plants) -> tuple[list[dict], list[dict]]:
    """Minimize the nominal SLICOT value over the tied PI grid, then replay on every case."""
    kps = [0.1, 0.2, 0.4, 0.8, 1.5, 3.0]
    kis = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 40.0]
    best = None
    rows = []
    n = 0
    t0 = time.perf_counter()
    for kp in kps:
        for ki in kis:
            n += 1
            blocks = diagonal_pi(kp, ki, kp, ki)
            ev = fast_eval(blocks, [plants[0]])
            if best is None or ev["worst_slicot"] < best[0]:
                best = (ev["worst_slicot"], blocks, [kp, ki])
    full = fast_eval(best[1], plants)
    rows.append(
        _record(
            "nominal_only_tied_pi",
            1,
            None,
            best[2] + best[2],
            full,
            time.perf_counter() - t0,
            n,
            f"selected on nominal slicot {best[0]}",
            full["worst_slicot"],
        )
    )
    return rows, best[1]


def per_case(plants) -> list[dict]:
    kps = [0.1, 0.2, 0.4, 0.8, 1.5, 3.0]
    kis = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 40.0]
    rows = []
    for plant in plants:
        best = None
        n = 0
        t0 = time.perf_counter()
        for kp in kps:
            for ki in kis:
                n += 1
                blocks = diagonal_pi(kp, ki, kp, ki)
                ev = fast_eval(blocks, [plant])
                if best is None or ev["worst_slicot"] < best[0]:
                    best = (ev["worst_slicot"], [kp, ki], ev)
        rows.append(
            _record(
                "per_case_tied_pi",
                plant["case_id"],
                None,
                best[1],
                best[2],
                time.perf_counter() - t0,
                n,
                "best tied PI for this case alone; not required to be shared",
                best[0],
            )
        )
    return rows


def _decouple(plant: dict) -> dict:
    out = dict(plant)
    phys = dict(plant["physical"])
    A = phys["A_scaled_physical"].copy()
    Bu = phys["Bu_pu"].copy()
    A[:18, 18:36] = 0.0
    A[18:36, :18] = 0.0
    Bu[:18, 2:] = 0.0
    Bu[18:36, :2] = 0.0
    phys["A_scaled_physical"] = A
    phys["Bu_pu"] = Bu
    Ag = plant["A"].copy()
    Bug = plant["Bu"].copy()
    Ag[:18, 18:36] = 0.0
    Ag[18:36, :18] = 0.0
    Bug[:18, 2:] = 0.0
    Bug[18:36, :2] = 0.0
    out["physical"] = phys
    out["A"] = Ag
    out["Bu"] = Bug
    out["case_id"] = plant["case_id"] + "_decoupled_design"
    return out


def decoupled_design(plants) -> tuple[list[dict], list[dict]]:
    design = [_decouple(p) for p in plants]
    rows, blocks, cost = pi_grid(design)
    for row in rows:
        row["family"] = "decoupled_design_" + row["family"]
    full = fast_eval(blocks, plants)
    rows.append(
        _record(
            "decoupled_design_replay_on_coupled",
            0,
            None,
            None,
            full,
            0.0,
            0,
            "best decoupled-grid controller replayed on the coupled plants",
            full["worst_slicot"],
        )
    )
    return rows, blocks


def centralized(plants) -> list[dict]:
    rows = []
    for plant in plants:
        t0 = time.perf_counter()
        synth = full_order_hinf(plant, gamma0=2.0)
        entry = {
            "family": "centralized_full_order_sb10ad",
            "case_id": plant["case_id"],
            "wall_s": time.perf_counter() - t0,
            "status": synth["status"],
            "gamma_reported": synth.get("gamma_reported"),
            "error": synth.get("error"),
            "n_controller": synth.get("n_controller"),
            "admissible": False,
        }
        if synth["status"] == "returned":
            K = {"A": synth["Ak"], "B": synth["Bk"], "C": synth["Ck"], "D": synth["Dk"]}
            cl = closed_loop_formula(plant, K)
            A_phys = physical_closed_loop(plant, K)
            phys = stability_report(A_phys)
            try:
                value = float(slicot_linf(cl["A"], cl["B"], cl["C"], cl["D"])["value"])
            except Exception as exc:
                value = None
                entry["slicot_error"] = str(exc)
            entry["independent_slicot"] = value
            entry["physical_abscissa"] = phys["abscissa"]
            entry["physical_class"] = classify_stability(phys, STAB_PASS)
        t1 = time.perf_counter()
        sf = state_feedback_feasible(plant, 1.0)
        entry["state_feedback_gamma_1_feasible"] = sf["feasible"]
        entry["state_feedback_status"] = sf["status"]
        entry["state_feedback_wall_s"] = time.perf_counter() - t1
        rows.append(entry)
    return rows


def run_campaign(plants) -> dict:
    t_all = time.perf_counter()
    records = []
    # Zero control law, stable unused poles.
    z = np.zeros((2, 2))
    zero = [
        block(-100.0 * np.eye(2), z, z, z, "converter_1"),
        block(-100.0 * np.eye(2), z, z, z, "converter_2"),
    ]
    records.append(_record("zero_stable_realization", 0, None, None, fast_eval(zero, plants), 0.0, 1, "Dk=0, Ak=-100 I", None))
    # Pathological realization: Ak=0, disconnected states on the axis.
    records.append(
        _record(
            "zero_ak_on_imaginary_axis",
            0,
            None,
            None,
            fast_eval(zero_blocks(), plants),
            0.0,
            1,
            "same zero law realized with Ak=0; internal stability must fail",
            None,
        )
    )
    records.extend(static_grid(plants))
    pi_rows, pi_best, pi_cost = pi_grid(plants)
    records.extend(pi_rows)
    # Negative sign, one competent-magnitude attempt.
    records.append(
        _record(
            "tied_diagonal_pi_wrong_sign",
            0,
            None,
            [-0.4, -4, -0.4, -4],
            fast_eval(diagonal_pi(-0.4, -4, -0.4, -4), plants),
            0.0,
            1,
            "sign flipped relative to the passing positive PI",
            None,
        )
    )
    de_rows, _ = optimize_family(
        plants,
        "tied_diagonal_pi_de",
        None,
        [(-1.0, 4.0), (-5.0, 60.0), (-1.0, 4.0), (-5.0, 60.0)],
        seeds=[1101, 1102],
        maxiter=8,
    )
    records.extend(de_rows)
    cross_rows, _ = optimize_family(
        plants,
        "tied_pi_with_cross_terms",
        None,
        [(-20.0, 20.0)] * 8,
        seeds=[1201],
        maxiter=6,
        pi_cross=True,
    )
    records.extend(cross_rows)
    gen_rows, _ = optimize_family(
        plants,
        "general_2state_multistart",
        None,
        [(-5.0, 5.0)] * 32,
        seeds=[1301, 1302, 1303, 1304],
        maxiter=25,
    )
    records.extend(gen_rows)
    nom_rows, _ = nominal_only(plants)
    records.extend(nom_rows)
    records.extend(per_case(plants))
    dec_rows, _ = decoupled_design(plants)
    records.extend(dec_rows)
    central = centralized(plants)
    return {
        "records": records,
        "centralized": central,
        "pi_grid_best_cost": pi_cost,
        "pi_grid_best_blocks": pi_best,
        "zero_blocks": zero,
        "wall_s": time.perf_counter() - t_all,
    }
