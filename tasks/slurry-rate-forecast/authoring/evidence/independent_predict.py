#!/usr/bin/env python3
"""Second route: Hanes lines at each temperature, then a Newton forward solve.

This file is not the oracle. It does not call a nonlinear least-squares solver
and it does not use the closed-form quadratic.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

R_GAS = 8.314462618
T_REF = 348.15
DATA = Path("/app/data")
OUT = Path("/app/output/predictions.json")


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def rate_newton(Cstar: float, kLa: float, ks: float, w: float, a: float, k: float, K: float) -> float:
    Rmt = 1.0 / kLa + 1.0 / ks
    ceiling = 0.999999 * Cstar / Rmt
    B = w * a * k
    rate = min(B * K * Cstar / (1.0 + K * Cstar), ceiling)
    for _ in range(60):
        surface = Cstar - rate * Rmt
        if surface <= 1e-14:
            rate *= 0.5
            continue
        catalytic = B * K * surface / (1.0 + K * surface)
        slope = B * K * (-Rmt) / (1.0 + K * surface) ** 2
        updated = rate - (catalytic - rate) / (slope - 1.0)
        if not math.isfinite(updated) or updated <= 0.0 or updated >= ceiling:
            updated = 0.5 * (rate + min(max(catalytic, 1e-16), ceiling))
        if abs(updated - rate) <= 1e-14 * max(1.0, abs(rate)):
            return float(updated)
        rate = updated
    return float(rate)


def build():
    vessel = json.loads((DATA / "vessel.json").read_text())
    solvent = json.loads((DATA / "solvent.json").read_text())
    catalyst = json.loads((DATA / "catalyst.json").read_text())
    volume = float(vessel["slurry_volume_m3"])
    p_atm = float(vessel["site_atmospheric_pressure_bar"])
    rho = float(catalyst["particle_density_kg_m3"])
    fresh = float(catalyst["fresh_co_uptake_umol_per_g"])
    A, B, C = (float(solvent["A"]), float(solvent["B"]), float(solvent["C"]))

    cal = read_csv("thermocouple_calibration.csv")
    tic = np.array([float(row["tic201_C"]) for row in cal], dtype=float)
    ref = np.array([float(row["reference_C"]) for row in cal], dtype=float)
    order = np.argsort(tic)
    tic, ref = tic[order], ref[order]

    sol = read_csv("h2_solubility.csv")
    sol_T = np.array([float(row["reference_C"]) + 273.15 for row in sol], dtype=float)
    sol_C = np.array([float(row["cstar_mol_m3_at_1bar"]) for row in sol], dtype=float)
    inv = 1.0 / sol_T
    order = np.argsort(inv)
    inv, lnC = inv[order], np.log(sol_C[order])

    standards = [
        (row["start_date"], row["end_date"], float(row["response_factor_mol_per_count"]))
        for row in read_csv("gc_standards.csv")
    ]
    kla = {}
    for row in read_csv("kla_measurements.csv"):
        if not row["N_rps"].strip():
            continue
        kla[round(float(row["N_rps"]), 5)] = float(row["kla_per_s"])
    kc = {}
    for row in read_csv("kc_measurements.csv"):
        if not row["N_rps"].strip() or row["note"].strip() == "baffles removed":
            continue
        kc[(round(float(row["N_rps"]), 5), round(float(row["dp_um"]), 5))] = float(row["kc_m_per_s"])
    uptake = {
        round(float(row["N_rps"]), 5): float(row["co_uptake_umol_per_g"])
        for row in read_csv("catalyst_sites.csv")
    }

    def temperature(tic201: float) -> float:
        return float(np.interp(tic201, tic, ref))

    def partial(T_C: float, p_gauge: float) -> float:
        T_K = T_C + 273.15
        return p_gauge + p_atm - 10 ** (A - B / (T_K + C))

    def interface(T_C: float, p_h2: float) -> float:
        return float(math.exp(float(np.interp(1.0 / (T_C + 273.15), inv, lnC))) * p_h2)

    runs = []
    for row in read_csv("rate_runs.csv"):
        T_C = temperature(float(row["tic201_C"]))
        p_h2 = partial(T_C, float(row["P_total_barg"]))
        speed = round(float(row["N_rps"]), 5)
        diameter = round(float(row["dp_um"]), 5)
        w = (float(row["m_cat_g"]) / 1000.0) / volume
        area = 6.0 * w / (rho * diameter * 1e-6)
        factor = [value for start, end, value in standards if start <= row["run_date"] <= end]
        rate = float(row["h2_peak_area"]) * factor[0] / float(row["window_s"]) / volume
        runs.append(
            (
                T_C + 273.15,
                interface(T_C, p_h2),
                kla[speed],
                kc[(speed, diameter)] * area,
                w,
                uptake[speed] / fresh,
                rate,
            )
        )

    cases = []
    for row in read_csv("predict_cases.csv"):
        T_C = float(row["bath_C"])
        speed = round(float(row["N_rps"]), 5)
        diameter = round(float(row["dp_um"]), 5)
        w = (float(row["m_cat_g"]) / 1000.0) / volume
        area = 6.0 * w / (rho * diameter * 1e-6)
        activity = 1.0 if row["catalyst_state"].strip() == "fresh" else uptake[speed] / fresh
        p_h2 = partial(T_C, float(row["P_total_barg"]))
        cases.append(
            {
                "case_id": row["case_id"],
                "T_K": T_C + 273.15,
                "Cstar": interface(T_C, p_h2),
                "kLa": kla[speed],
                "ks": kc[(speed, diameter)] * area,
                "w": w,
                "a": activity,
            }
        )
    return runs, cases


def hanes(runs):
    groups: dict[float, list[tuple[float, float]]] = defaultdict(list)
    for T_K, Cstar, kLa, ks, w, activity, rate in runs:
        surface = Cstar - rate * (1.0 / kLa + 1.0 / ks)
        specific = rate / (w * activity)
        groups[round(T_K, 3)].append((surface, surface / specific))
    temperatures = []
    surface_coefficients = []
    adsorption = []
    for T_K in sorted(groups):
        xs = np.array([pair[0] for pair in groups[T_K]])
        ys = np.array([pair[1] for pair in groups[T_K]])
        slope, intercept = np.polyfit(xs, ys, 1)
        temperatures.append(T_K)
        surface_coefficients.append(1.0 / slope)
        adsorption.append(slope / intercept)
    temperatures = np.array(temperatures)
    reduced = 1.0 / temperatures - 1.0 / T_REF
    slope_k, intercept_k = np.polyfit(reduced, np.log(np.array(surface_coefficients)), 1)
    slope_K, intercept_K = np.polyfit(reduced, np.log(np.array(adsorption)), 1)
    return math.exp(intercept_k), -slope_k * R_GAS, math.exp(intercept_K), -slope_K * R_GAS


def main() -> None:
    runs, cases = build()
    k_ref, activation, K_ref, enthalpy = hanes(runs)
    rows = []
    for case in cases:
        k = k_ref * math.exp(-activation / R_GAS * (1.0 / case["T_K"] - 1.0 / T_REF))
        K = K_ref * math.exp(-enthalpy / R_GAS * (1.0 / case["T_K"] - 1.0 / T_REF))
        rate = rate_newton(case["Cstar"], case["kLa"], case["ks"], case["w"], case["a"], k, K)
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
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"cases": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
