#!/usr/bin/env python3
"""Fit ablations and naive baselines against the sealed latent truth.

Reads the shipped campaign files and tests/truth/expected.json. Does not
rewrite either.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "environment" / "data"
TRUTH_PATH = ROOT / "tests" / "truth" / "expected.json"
OUT_DIR = Path(__file__).resolve().parent

R_GAS = 8.314462618
T_REF = 348.15
REL_TOL = 0.040
ABS_TOL = 0.0006
FRAC_TOL = 0.030
KEYS = (
    "resistance_fraction_gas_liquid",
    "resistance_fraction_liquid_solid",
    "resistance_fraction_kinetic",
)


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def rate_closed(Cstar, kLa, ks, w, a, k, K) -> float:
    Rmt = 1.0 / kLa + 1.0 / ks
    B = w * a * k
    Aq = K * Rmt
    Bq = -(1.0 + K * Cstar + B * K * Rmt)
    Cq = B * K * Cstar
    disc = max(0.0, Bq * Bq - 4.0 * Aq * Cq)
    return (2.0 * Cq) / (-Bq + math.sqrt(disc))


def load(use_tc=True, use_gc=True, use_pvap=True, use_activity=True, use_mt=True):
    vessel = json.loads((DATA / "vessel.json").read_text())
    solvent = json.loads((DATA / "solvent.json").read_text())
    catalyst = json.loads((DATA / "catalyst.json").read_text())
    volume = float(vessel["slurry_volume_m3"])
    p_atm = float(vessel["site_atmospheric_pressure_bar"])
    rho = float(catalyst["particle_density_kg_m3"])
    fresh = float(catalyst["fresh_co_uptake_umol_per_g"])
    A, B, C = (float(solvent["A"]), float(solvent["B"]), float(solvent["C"]))
    cal = read_csv("thermocouple_calibration.csv")
    tic = np.array([float(row["tic201_C"]) for row in cal])
    ref = np.array([float(row["reference_C"]) for row in cal])
    order = np.argsort(tic)
    tic, ref = tic[order], ref[order]
    sol = read_csv("h2_solubility.csv")
    sol_T = np.array([float(row["reference_C"]) + 273.15 for row in sol])
    sol_C = np.array([float(row["cstar_mol_m3_at_1bar"]) for row in sol])
    inv = 1.0 / sol_T
    order = np.argsort(inv)
    inv, lnC = inv[order], np.log(sol_C[order])
    standards = [
        (row["start_date"], row["end_date"], float(row["response_factor_mol_per_count"]))
        for row in read_csv("gc_standards.csv")
    ]
    kla, kc, uptake = {}, {}, {}
    for row in read_csv("kla_measurements.csv"):
        if row["N_rps"].strip():
            kla[round(float(row["N_rps"]), 5)] = float(row["kla_per_s"])
    for row in read_csv("kc_measurements.csv"):
        if row["N_rps"].strip() and row["note"].strip() != "baffles removed":
            kc[(round(float(row["N_rps"]), 5), round(float(row["dp_um"]), 5))] = float(row["kc_m_per_s"])
    for row in read_csv("catalyst_sites.csv"):
        uptake[round(float(row["N_rps"]), 5)] = float(row["co_uptake_umol_per_g"])

    def temperature(tic201: float) -> float:
        if use_tc:
            return float(np.interp(tic201, tic, ref))
        return tic201

    def interface(T_C: float, p_gauge: float) -> float:
        T_K = T_C + 273.15
        vapor = 10 ** (A - B / (T_K + C)) if use_pvap else 0.0
        p_h2 = p_gauge + p_atm - vapor
        return float(math.exp(float(np.interp(1.0 / T_K, inv, lnC))) * p_h2)

    def pack(T_C, p_gauge, speed, diameter, mass, state, rate):
        w = (mass / 1000.0) / volume
        area = 6.0 * w / (rho * diameter * 1e-6)
        if use_activity and state == "conditioned":
            activity = uptake[round(speed, 5)] / fresh
        else:
            activity = 1.0
        kLa = kla[round(speed, 5)] if use_mt else 1e6
        solid = kc[(round(speed, 5), round(diameter, 5))] * area if use_mt else 1e6
        return {
            "T_K": T_C + 273.15,
            "Cstar": interface(T_C, p_gauge),
            "kLa": kLa,
            "ks": solid,
            "w": w,
            "a": activity,
            "rate": rate,
        }

    runs = []
    for row in read_csv("rate_runs.csv"):
        T_C = temperature(float(row["tic201_C"]))
        if use_gc:
            factor = [value for start, end, value in standards if start <= row["run_date"] <= end][0]
            rate = float(row["h2_peak_area"]) * factor / float(row["window_s"]) / volume
        else:
            rate = float(row["sheet_rate_mol_m3_s"])
        runs.append(
            pack(T_C, float(row["P_total_barg"]), float(row["N_rps"]), float(row["dp_um"]), float(row["m_cat_g"]), "conditioned", rate)
        )
    cases = []
    for row in read_csv("predict_cases.csv"):
        packed = pack(
            float(row["bath_C"]),
            float(row["P_total_barg"]),
            float(row["N_rps"]),
            float(row["dp_um"]),
            float(row["m_cat_g"]),
            row["catalyst_state"].strip(),
            0.0,
        )
        packed["case_id"] = row["case_id"]
        cases.append(packed)
    return runs, cases


def fit_lh(runs, fix_enthalpy=None):
    def residuals(theta):
        if fix_enthalpy is None:
            ln_k, e_kJ, ln_K, dh_kJ = theta
        else:
            ln_k, e_kJ, ln_K = theta
            dh_kJ = fix_enthalpy
        k_ref, K_ref = math.exp(ln_k), math.exp(ln_K)
        activation, enthalpy = e_kJ * 1000.0, dh_kJ * 1000.0
        out = []
        for run in runs:
            k = k_ref * math.exp(-activation / R_GAS * (1.0 / run["T_K"] - 1.0 / T_REF))
            K = K_ref * math.exp(-enthalpy / R_GAS * (1.0 / run["T_K"] - 1.0 / T_REF))
            pred = rate_closed(run["Cstar"], run["kLa"], run["ks"], run["w"], run["a"], k, K)
            out.append((pred - run["rate"]) / run["rate"])
        return np.asarray(out)

    start = [math.log(0.01), 55.0, math.log(0.3)]
    if fix_enthalpy is None:
        start.append(-35.0)
    solved = least_squares(residuals, np.array(start), method="trf", xtol=1e-12, ftol=1e-12, gtol=1e-12)
    theta = [float(v) for v in solved.x]
    if fix_enthalpy is None:
        ln_k, e_kJ, ln_K, dh_kJ = theta
    else:
        ln_k, e_kJ, ln_K = theta
        dh_kJ = fix_enthalpy
    return math.exp(ln_k), e_kJ * 1000.0, math.exp(ln_K), dh_kJ * 1000.0


def predict_lh(cases, params):
    k_ref, activation, K_ref, enthalpy = params
    rows = []
    for case in cases:
        k = k_ref * math.exp(-activation / R_GAS * (1.0 / case["T_K"] - 1.0 / T_REF))
        K = K_ref * math.exp(-enthalpy / R_GAS * (1.0 / case["T_K"] - 1.0 / T_REF))
        rate = rate_closed(case["Cstar"], case["kLa"], case["ks"], case["w"], case["a"], k, K)
        rows.append(fractions(case["case_id"], case["Cstar"], case["kLa"], case["ks"], rate))
    return rows


def fit_first(runs):
    def residuals(theta):
        k_ref = math.exp(theta[0])
        activation = theta[1] * 1000.0
        out = []
        for run in runs:
            k = k_ref * math.exp(-activation / R_GAS * (1.0 / run["T_K"] - 1.0 / T_REF))
            Rmt = 1.0 / run["kLa"] + 1.0 / run["ks"]
            B = run["w"] * run["a"] * k
            pred = B * run["Cstar"] / (1.0 + B * Rmt)
            out.append((pred - run["rate"]) / run["rate"])
        return np.asarray(out)

    solved = least_squares(residuals, np.array([-6.5, 60.0]), method="trf")
    return float(solved.x[0]), float(solved.x[1])


def predict_first(cases, theta):
    ln_k, e_kJ = theta
    rows = []
    for case in cases:
        k = math.exp(ln_k) * math.exp(-e_kJ * 1000.0 / R_GAS * (1.0 / case["T_K"] - 1.0 / T_REF))
        Rmt = 1.0 / case["kLa"] + 1.0 / case["ks"]
        B = case["w"] * case["a"] * k
        rate = B * case["Cstar"] / (1.0 + B * Rmt)
        rows.append(fractions(case["case_id"], case["Cstar"], case["kLa"], case["ks"], rate))
    return rows


def fit_power(runs):
    def residuals(theta):
        k_ref = math.exp(theta[0])
        activation = theta[1] * 1000.0
        order = theta[2]
        out = []
        for run in runs:
            k = k_ref * math.exp(-activation / R_GAS * (1.0 / run["T_K"] - 1.0 / T_REF))
            Rmt = 1.0 / run["kLa"] + 1.0 / run["ks"]
            coef = run["w"] * run["a"] * k
            rate = run["rate"]
            ceiling = 0.999 * run["Cstar"] / Rmt
            for _ in range(40):
                surface = run["Cstar"] - rate * Rmt
                if surface <= 1e-14:
                    rate *= 0.5
                    continue
                catalytic = coef * surface**order
                slope = coef * order * surface ** (order - 1.0) * (-Rmt)
                updated = rate - (catalytic - rate) / (slope - 1.0)
                if updated <= 0.0 or updated >= ceiling:
                    updated = 0.5 * (rate + min(catalytic, ceiling))
                if abs(updated - rate) < 1e-14:
                    rate = updated
                    break
                rate = updated
            out.append((rate - run["rate"]) / run["rate"])
        return np.asarray(out)

    solved = least_squares(
        residuals,
        np.array([-5.0, 50.0, 0.5]),
        bounds=([-20.0, 5.0, 0.05], [5.0, 150.0, 1.5]),
        method="trf",
    )
    return [float(v) for v in solved.x]


def predict_power(cases, theta):
    ln_k, e_kJ, order = theta
    rows = []
    for case in cases:
        k = math.exp(ln_k) * math.exp(-e_kJ * 1000.0 / R_GAS * (1.0 / case["T_K"] - 1.0 / T_REF))
        Rmt = 1.0 / case["kLa"] + 1.0 / case["ks"]
        coef = case["w"] * case["a"] * k
        rate = 0.02
        ceiling = 0.999 * case["Cstar"] / Rmt
        for _ in range(40):
            surface = case["Cstar"] - rate * Rmt
            if surface <= 1e-14:
                rate *= 0.5
                continue
            catalytic = coef * surface**order
            slope = coef * order * surface ** (order - 1.0) * (-Rmt)
            updated = rate - (catalytic - rate) / (slope - 1.0)
            if updated <= 0.0 or updated >= ceiling:
                updated = 0.5 * (rate + min(catalytic, ceiling))
            if abs(updated - rate) < 1e-14:
                rate = updated
                break
            rate = updated
        rows.append(fractions(case["case_id"], case["Cstar"], case["kLa"], case["ks"], rate))
    return rows


def fractions(case_id, Cstar, kLa, ks, rate):
    drop_gl = rate / kLa
    drop_ls = rate / ks
    surface = Cstar - drop_gl - drop_ls
    return {
        "case_id": case_id,
        "rate_mol_m3_s": rate,
        "resistance_fraction_gas_liquid": drop_gl / Cstar,
        "resistance_fraction_liquid_solid": drop_ls / Cstar,
        "resistance_fraction_kinetic": surface / Cstar,
    }


def judge(prediction, truth_by_id):
    per_case = []
    for row in prediction:
        truth = truth_by_id[row["case_id"]]
        rate_err = abs(row["rate_mol_m3_s"] - truth["rate_mol_m3_s"])
        rate_tol = max(REL_TOL * abs(truth["rate_mol_m3_s"]), ABS_TOL)
        frac_err = {key: abs(row[key] - truth[key]) for key in KEYS}
        rate_ok = math.isfinite(row["rate_mol_m3_s"]) and rate_err <= rate_tol
        frac_ok = all(math.isfinite(row[key]) and frac_err[key] <= FRAC_TOL for key in KEYS)
        per_case.append(
            {
                "case_id": row["case_id"],
                "abs_rate_error": rate_err,
                "rel_rate_error": rate_err / abs(truth["rate_mol_m3_s"]),
                "rate_tolerance": rate_tol,
                "fraction_abs_error": frac_err,
                "case_passes": bool(rate_ok and frac_ok),
            }
        )
    return {
        "max_rel_rate_error": max(item["rel_rate_error"] for item in per_case),
        "max_fraction_abs_error": max(max(item["fraction_abs_error"].values()) for item in per_case),
        "cases_passed": [item["case_id"] for item in per_case if item["case_passes"]],
        "cases_failed": [item["case_id"] for item in per_case if not item["case_passes"]],
        "every_case_fails": all(not item["case_passes"] for item in per_case),
        "cases": per_case,
    }


def constant_prediction(case_ids, rate_value, fractions):
    return [
        {
            "case_id": case_id,
            "rate_mol_m3_s": rate_value,
            "resistance_fraction_gas_liquid": fractions[0],
            "resistance_fraction_liquid_solid": fractions[1],
            "resistance_fraction_kinetic": fractions[2],
        }
        for case_id in case_ids
    ]


def main() -> None:
    truth = json.loads(TRUTH_PATH.read_text())["cases"]
    truth_by_id = {row["case_id"]: row for row in truth}
    case_ids = [row["case_id"] for row in truth]
    rates = np.array([row["rate_mol_m3_s"] for row in truth])
    frac_means = [float(np.mean([row[key] for row in truth])) for key in KEYS]
    frac_meds = [float(np.median([row[key] for row in truth])) for key in KEYS]

    oracle = json.loads(Path("/tmp/oracle.json").read_text())["cases"]
    indep = json.loads(Path("/tmp/indep.json").read_text())["cases"]

    results = {
        "tolerances": {
            "rate_relative": REL_TOL,
            "rate_absolute_mol_m3_s": ABS_TOL,
            "fraction_absolute": FRAC_TOL,
            "rule": "abs(rate error) <= max(rate_relative * abs(truth), rate_absolute); abs(fraction error) <= fraction_absolute",
        },
        "correct": {
            "output_error_nls": judge(oracle, truth_by_id),
            "hanes_newton": judge(indep, truth_by_id),
        },
    }

    ablations = {}
    specs = [
        ("no_thermocouple_calibration", dict(use_tc=False)),
        ("spreadsheet_rates", dict(use_gc=False)),
        ("no_solvent_vapor_pressure", dict(use_pvap=False)),
        ("no_uptake_scaling", dict(use_activity=False)),
        ("no_film_resistances", dict(use_mt=False)),
    ]
    for name, flags in specs:
        runs, cases = load(**flags)
        ablations[name] = judge(predict_lh(cases, fit_lh(runs)), truth_by_id)
    runs, cases = load()
    ablations["first_order"] = judge(predict_first(cases, fit_first(runs)), truth_by_id)
    ablations["power_law"] = judge(predict_power(cases, fit_power(runs)), truth_by_id)
    ablations["adsorption_coefficient_independent_of_temperature"] = judge(
        predict_lh(cases, fit_lh(runs, fix_enthalpy=0.0)), truth_by_id
    )
    results["ablations"] = ablations

    baselines = {}
    mean_rate = float(rates.mean())
    median_rate = float(np.median(rates))
    named = {
        "constant_0": constant_prediction(case_ids, 0.0, (0.0, 0.0, 0.0)),
        "constant_1": constant_prediction(case_ids, 1.0, (1.0, 1.0, 1.0)),
        "mean_of_true_values": constant_prediction(case_ids, mean_rate, tuple(frac_means)),
        "median_of_true_values": constant_prediction(case_ids, median_rate, tuple(frac_meds)),
    }
    for index, rate in enumerate(rates):
        named[f"constant_equal_to_truth_rate_{case_ids[index]}"] = constant_prediction(
            case_ids, float(rate), (float(rate), float(rate), float(rate))
        )
    predict_rows = read_csv("predict_cases.csv")
    for column in ("bath_C", "P_total_barg", "N_rps", "m_cat_g", "dp_um"):
        copied = []
        for row in predict_rows:
            value = float(row[column])
            copied.append(
                {
                    "case_id": row["case_id"],
                    "rate_mol_m3_s": value,
                    "resistance_fraction_gas_liquid": value,
                    "resistance_fraction_liquid_solid": value,
                    "resistance_fraction_kinetic": value,
                }
            )
        named[f"copy_{column}"] = copied
    raw_baselines = {}
    for name, prediction in named.items():
        baselines[name] = judge(prediction, truth_by_id)
        raw_baselines[name] = {"cases": prediction}
        if not baselines[name]["every_case_fails"]:
            print(f"BASELINE LEAK {name} passed {baselines[name]['cases_passed']}", file=sys.stderr)
    results["baselines"] = baselines
    (OUT_DIR / "calibration_results.json").write_text(json.dumps(results["correct"], indent=2) + "\n")
    (OUT_DIR / "ablation_outputs.json").write_text(json.dumps(ablations, indent=2) + "\n")
    (OUT_DIR / "baseline_failures.json").write_text(json.dumps(baselines, indent=2) + "\n")
    (OUT_DIR / "baseline_predictions.json").write_text(json.dumps(raw_baselines, indent=2) + "\n")
    print("correct NLS max rel", results["correct"]["output_error_nls"]["max_rel_rate_error"])
    print("correct Hanes max rel", results["correct"]["hanes_newton"]["max_rel_rate_error"])
    print("correct NLS max frac", results["correct"]["output_error_nls"]["max_fraction_abs_error"])
    print("correct Hanes max frac", results["correct"]["hanes_newton"]["max_fraction_abs_error"])
    for name, summary in ablations.items():
        print(
            f"ablation {name}: max rel {summary['max_rel_rate_error']:.4f} "
            f"fail {summary['cases_failed']} pass {summary['cases_passed']}"
        )
    leaks = [name for name, summary in baselines.items() if not summary["every_case_fails"]]
    print("baseline leaks", leaks)


if __name__ == "__main__":
    main()
