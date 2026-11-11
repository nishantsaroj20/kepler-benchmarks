#!/usr/bin/env python3
"""Oracle: output-error fit of a two-film Langmuir-Hinshelwood rate, closed form."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

R_GAS = 8.314462618
T_REF = 348.15
DATA = Path("/app/data")
OUT = Path("/app/output/predictions.json")


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def rate_closed(Cstar: float, kLa: float, ks: float, w: float, a: float, k: float, K: float) -> float:
    Rmt = 1.0 / kLa + 1.0 / ks
    B = w * a * k
    Aq = K * Rmt
    Bq = -(1.0 + K * Cstar + B * K * Rmt)
    Cq = B * K * Cstar
    disc = Bq * Bq - 4.0 * Aq * Cq
    if disc < 0.0 and disc > -1e-18:
        disc = 0.0
    return (2.0 * Cq) / (-Bq + math.sqrt(disc))


def build():
    vessel = json.loads((DATA / "vessel.json").read_text())
    solvent = json.loads((DATA / "solvent.json").read_text())
    catalyst = json.loads((DATA / "catalyst.json").read_text())
    volume = float(vessel["slurry_volume_m3"])
    p_atm = float(vessel["site_atmospheric_pressure_bar"])
    rho = float(catalyst["particle_density_kg_m3"])
    fresh = float(catalyst["fresh_co_uptake_umol_per_g"])
    antoine = (float(solvent["A"]), float(solvent["B"]), float(solvent["C"]))

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

    def reference_temperature(tic201: float) -> float:
        return float(np.interp(tic201, tic, ref))

    def cstar(T_C: float, p_h2: float) -> float:
        T_K = T_C + 273.15
        return float(math.exp(np.interp(1.0 / T_K, inv, lnC)) * p_h2)

    def vapor(T_C: float) -> float:
        A, B, C = antoine
        T_K = T_C + 273.15
        return 10 ** (A - B / (T_K + C))

    def partial(T_C: float, p_gauge: float) -> float:
        return p_gauge + p_atm - vapor(T_C)

    standards = []
    for row in read_csv("gc_standards.csv"):
        standards.append(
            (row["start_date"], row["end_date"], float(row["response_factor_mol_per_count"]))
        )

    def response_factor(day: str) -> float:
        hits = [factor for start, end, factor in standards if start <= day <= end]
        if len(hits) != 1:
            raise RuntimeError(f"response factor is not unique on {day}")
        return hits[0]

    kla = {}
    for row in read_csv("kla_measurements.csv"):
        if row["N_rps"].strip() == "":
            continue
        kla[round(float(row["N_rps"]), 5)] = float(row["kla_per_s"])

    kc = {}
    for row in read_csv("kc_measurements.csv"):
        if row["N_rps"].strip() == "" or row["note"].strip() == "baffles removed":
            continue
        key = (round(float(row["N_rps"]), 5), round(float(row["dp_um"]), 5))
        if key in kc:
            raise RuntimeError(f"duplicate liquid-solid row {key}")
        kc[key] = float(row["kc_m_per_s"])

    uptake = {}
    for row in read_csv("catalyst_sites.csv"):
        uptake[round(float(row["N_rps"]), 5)] = float(row["co_uptake_umol_per_g"])

    runs = []
    for row in read_csv("rate_runs.csv"):
        T_C = reference_temperature(float(row["tic201_C"]))
        p_h2 = partial(T_C, float(row["P_total_barg"]))
        C = cstar(T_C, p_h2)
        speed = round(float(row["N_rps"]), 5)
        diameter = round(float(row["dp_um"]), 5)
        w = (float(row["m_cat_g"]) / 1000.0) / volume
        area = 6.0 * w / (rho * (float(row["dp_um"]) * 1e-6))
        moles = float(row["h2_peak_area"]) * response_factor(row["run_date"])
        rate = moles / float(row["window_s"]) / volume
        runs.append(
            {
                "T_K": T_C + 273.15,
                "Cstar": C,
                "kLa": kla[speed],
                "ks": kc[(speed, diameter)] * area,
                "w": w,
                "a": uptake[speed] / fresh,
                "rate": rate,
            }
        )

    cases = []
    for row in read_csv("predict_cases.csv"):
        T_C = float(row["bath_C"])
        p_h2 = partial(T_C, float(row["P_total_barg"]))
        speed = round(float(row["N_rps"]), 5)
        diameter = round(float(row["dp_um"]), 5)
        w = (float(row["m_cat_g"]) / 1000.0) / volume
        area = 6.0 * w / (rho * (float(row["dp_um"]) * 1e-6))
        state = row["catalyst_state"].strip()
        activity = 1.0 if state == "fresh" else uptake[speed] / fresh
        cases.append(
            {
                "case_id": row["case_id"],
                "T_K": T_C + 273.15,
                "Cstar": cstar(T_C, p_h2),
                "kLa": kla[speed],
                "ks": kc[(speed, diameter)] * area,
                "w": w,
                "a": activity,
            }
        )
    return runs, cases


def fit(runs: list[dict]) -> tuple[float, float, float, float]:
    def residuals(theta: np.ndarray) -> np.ndarray:
        ln_k, e_kJ, ln_K, dh_kJ = theta
        k_ref = math.exp(ln_k)
        K_ref = math.exp(ln_K)
        activation = e_kJ * 1000.0
        enthalpy = dh_kJ * 1000.0
        out = []
        for run in runs:
            k = k_ref * math.exp(-activation / R_GAS * (1.0 / run["T_K"] - 1.0 / T_REF))
            K = K_ref * math.exp(-enthalpy / R_GAS * (1.0 / run["T_K"] - 1.0 / T_REF))
            pred = rate_closed(run["Cstar"], run["kLa"], run["ks"], run["w"], run["a"], k, K)
            out.append(pred - run["rate"])
        return np.asarray(out)

    start = np.array([math.log(0.01), 55.0, math.log(0.3), -35.0])
    solved = least_squares(
        residuals,
        start,
        method="trf",
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=200,
    )
    if not solved.success:
        raise RuntimeError(solved.message)
    ln_k, e_kJ, ln_K, dh_kJ = (float(v) for v in solved.x)
    return math.exp(ln_k), e_kJ * 1000.0, math.exp(ln_K), dh_kJ * 1000.0


def predict(cases: list[dict], k_ref: float, activation: float, K_ref: float, enthalpy: float) -> dict:
    rows = []
    for case in cases:
        k = k_ref * math.exp(-activation / R_GAS * (1.0 / case["T_K"] - 1.0 / T_REF))
        K = K_ref * math.exp(-enthalpy / R_GAS * (1.0 / case["T_K"] - 1.0 / T_REF))
        rate = rate_closed(case["Cstar"], case["kLa"], case["ks"], case["w"], case["a"], k, K)
        drop_gl = rate / case["kLa"]
        drop_ls = rate / case["ks"]
        surface = case["Cstar"] - drop_gl - drop_ls
        rows.append(
            {
                "case_id": case["case_id"],
                "rate_mol_m3_s": rate,
                "resistance_fraction_gas_liquid": drop_gl / case["Cstar"],
                "resistance_fraction_liquid_solid": drop_ls / case["Cstar"],
                "resistance_fraction_kinetic": surface / case["Cstar"],
            }
        )
    return {"cases": rows}


def main() -> None:
    runs, cases = build()
    params = fit(runs)
    payload = predict(cases, *params)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
