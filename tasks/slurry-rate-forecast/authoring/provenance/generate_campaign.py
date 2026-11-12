#!/usr/bin/env python3
"""Write the visible campaign and the shifted hidden campaign.

The hidden forecast truth is this latent forward model on a finer grid than
the reference program. It is not the reference fit. The hidden grid is a
different regime, solvent, instrument model, and operating map. A new noise
seed is not the shift.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("forecast", ROOT / "solution" / "forecast.py")
forecast = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(forecast)

R_GAS = 8.314462618
T_REF = 348.15
TRUTH_NODES = 80
TRUTH_STEPS = 200

VISIBLE_SEED = 20261111
HIDDEN_SEED = 20261112


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def dump_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def visible_campaign() -> None:
    data = ROOT / "environment" / "data"
    k_ref, kads_ref, energy, enthalpy = 0.0125, 0.28, 62000.0, -42000.0
    cstar_ref, dh_sol = 2.85, 4800.0
    rho, fresh = 1900.0, 186.0
    slurry, p_atm = 0.0015, 1.01325
    antoine = (4.07827, 1343.943, -53.773)
    diffusivity = 2.5e-9
    shell_um = 3.0

    def psat(temp_k: float) -> float:
        a_coef, b_coef, c_coef = antoine
        return 10 ** (a_coef - b_coef / (temp_k + c_coef))

    def c_one(temp_c: float) -> float:
        temp_k = temp_c + 273.15
        return cstar_ref * math.exp(-dh_sol / R_GAS * (1.0 / temp_k - 1.0 / T_REF))

    def tic_of(temp_c: float) -> float:
        return temp_c + 4.6 + 0.018 * (temp_c - 75.0)

    def kla(speed: float) -> float:
        return 0.048 * (speed / 12.0) ** 0.82

    def kc(speed: float, diameter: float) -> float:
        return 6.8e-5 * math.sqrt(speed / 12.0) * math.sqrt(80.0 / diameter)

    def conditioned(speed: float) -> float:
        return fresh / (1.0 + (speed / 21.0) ** 2)

    dump_json(data / "vessel.json", {
        "slurry_volume_m3": slurry,
        "site_atmospheric_pressure_bar": p_atm,
        "baffled": True,
        "working_volume_description": "1.5 L slurry, gas headspace above a single pitched-blade impeller",
    })
    dump_json(data / "solvent.json", {
        "name": "toluene",
        "antoine_form": "log10(P_sat_bar) = A - B / (T_K + C)",
        "A": antoine[0], "B": antoine[1], "C": antoine[2],
        "T_unit": "kelvin", "P_unit": "bar",
    })
    dump_json(data / "catalyst.json", {
        "metal": "palladium",
        "support": "activated carbon",
        "metal_location": "eggshell",
        "eggshell_depth_um": shell_um,
        "effective_diffusivity_m2_per_s": diffusivity,
        "particle_density_kg_m3": rho,
        "fresh_co_uptake_umol_per_g": fresh,
    })
    speeds = (6.0, 9.0, 12.0, 16.0)
    diameters = (40.0, 80.0, 150.0)
    write_csv(data / "kla_measurements.csv", ["N_rps", "kla_per_s", "note"], [
        {"N_rps": f"{speed:g}", "kla_per_s": f"{kla(speed):.8e}", "note": ""} for speed in speeds
    ] + [{"N_rps": "", "kla_per_s": "4.10000000e-02", "note": "tachometer blank"}])
    kc_rows = []
    for speed in speeds:
        for diameter in diameters:
            kc_rows.append({"N_rps": f"{speed:g}", "dp_um": f"{diameter:g}", "kc_m_per_s": f"{kc(speed, diameter):.8e}", "note": ""})
    kc_rows.append({"N_rps": "12", "dp_um": "80", "kc_m_per_s": f"{kc(12, 80) * 1.65:.8e}", "note": "baffles removed"})
    write_csv(data / "kc_measurements.csv", ["N_rps", "dp_um", "kc_m_per_s", "note"], kc_rows)
    write_csv(data / "catalyst_sites.csv", ["catalyst_state", "N_rps", "co_uptake_umol_per_g"], [
        {"catalyst_state": "conditioned", "N_rps": f"{speed:g}", "co_uptake_umol_per_g": f"{conditioned(speed):.6f}"}
        for speed in speeds
    ] + [{"catalyst_state": "fresh", "N_rps": "", "co_uptake_umol_per_g": f"{fresh:.6f}"}])
    write_csv(data / "thermocouple_calibration.csv", ["tic201_C", "reference_C"], [
        {"tic201_C": f"{tic_of(temp):.6f}", "reference_C": f"{temp:.4f}"} for temp in range(40, 135, 5)
    ])
    write_csv(data / "h2_solubility.csv", ["reference_C", "cstar_mol_m3_at_1bar"], [
        {"reference_C": f"{temp:.4f}", "cstar_mol_m3_at_1bar": f"{c_one(temp):.8e}"} for temp in range(40, 132, 2)
    ])
    write_csv(data / "gc_standards.csv", ["start_date", "end_date", "response_factor_mol_per_count"], [
        {"start_date": "2024-04-01", "end_date": "2024-04-15", "response_factor_mol_per_count": "2.62000000e-06"},
        {"start_date": "2024-04-16", "end_date": "2024-04-30", "response_factor_mol_per_count": "2.02000000e-06"},
    ])
    runs = [
        ("R01", "2024-04-18", 60, 0.8, 12, 22, 80),
        ("R02", "2024-04-09", 60, 1.6, 12, 22, 80),
        ("R03", "2024-04-19", 60, 3.0, 12, 22, 80),
        ("R04", "2024-04-08", 60, 5.0, 12, 22, 80),
        ("R05", "2024-04-11", 60, 7.0, 12, 22, 80),
        ("R06", "2024-04-17", 75, 0.8, 12, 22, 80),
        ("R07", "2024-04-10", 75, 1.6, 12, 22, 80),
        ("R08", "2024-04-18", 75, 3.0, 12, 22, 80),
        ("R09", "2024-04-07", 75, 5.0, 12, 22, 80),
        ("R10", "2024-04-12", 75, 7.0, 12, 22, 80),
        ("R11", "2024-04-20", 90, 0.8, 12, 22, 80),
        ("R12", "2024-04-08", 90, 1.6, 12, 22, 80),
        ("R13", "2024-04-22", 90, 3.0, 12, 22, 80),
        ("R14", "2024-04-09", 90, 5.0, 12, 22, 80),
        ("R15", "2024-04-14", 90, 7.0, 12, 22, 80),
        ("R16", "2024-04-11", 75, 3.0, 12, 8, 80),
        ("R17", "2024-04-23", 75, 3.0, 12, 14, 80),
        ("R18", "2024-04-06", 75, 3.0, 12, 35, 80),
        ("R19", "2024-04-24", 75, 3.0, 12, 50, 80),
        ("R20", "2024-04-13", 75, 3.0, 6, 22, 80),
        ("R21", "2024-04-16", 75, 3.0, 9, 22, 80),
        ("R22", "2024-04-05", 75, 3.0, 16, 22, 80),
        ("R23", "2024-04-15", 75, 3.0, 12, 22, 40),
        ("R24", "2024-04-25", 75, 3.0, 12, 22, 150),
        ("R25", "2024-04-04", 75, 3.0, 6, 8, 40),
        ("R26", "2024-04-26", 75, 3.0, 16, 8, 40),
        ("R27", "2024-04-03", 75, 2.0, 6, 50, 150),
        ("R28", "2024-04-27", 60, 3.0, 6, 40, 150),
        ("R29", "2024-04-02", 90, 2.0, 16, 10, 40),
        ("R30", "2024-04-28", 75, 3.0, 12, 22, 80),
    ]
    run_rows = []
    for run_id, date, temp, gauge, speed, mass, diameter in runs:
        run_rows.append({
            "run_id": run_id, "run_date": date, "status": "ok",
            "tic201_C": f"{tic_of(temp):.6f}", "P_total_barg": f"{gauge:.4f}",
            "N_rps": f"{speed:g}", "m_cat_g": f"{mass:.4f}", "dp_um": f"{diameter:g}",
            "h2_peak_area": "1.0", "window_s": "600.0", "sheet_rate_mol_m3_s": "0.0",
        })
    run_rows.append({
        "run_id": "R31", "run_date": "2024-04-12", "status": "aborted",
        "tic201_C": f"{tic_of(75):.6f}", "P_total_barg": "3.0000",
        "N_rps": "12", "m_cat_g": "22.0000", "dp_um": "80",
        "h2_peak_area": "9.50000000e+06", "window_s": "600.0", "sheet_rate_mol_m3_s": "1.20000000e+01",
    })
    write_csv(data / "rate_runs.csv", list(run_rows[0].keys()), run_rows)
    cases = [
        ("H1", 75.0, 11.5, 12.0, 40.0, 40.0, "conditioned"),
        ("L1", 75.0, 0.25, 12.0, 22.0, 80.0, "conditioned"),
        ("T1", 114.0, 1.55, 12.0, 10.0, 80.0, "conditioned"),
        ("W1", 75.0, 2.0, 6.0, 85.0, 80.0, "conditioned"),
        ("D1", 60.0, 5.5, 9.0, 36.0, 150.0, "conditioned"),
        ("F1", 75.0, 4.5, 16.0, 10.0, 40.0, "fresh"),
        ("X1", 108.0, 6.5, 6.0, 78.0, 150.0, "fresh"),
        ("K1", 60.0, 0.85, 16.0, 7.0, 40.0, "conditioned"),
    ]
    write_csv(data / "predict_cases.csv", ["case_id", "bath_C", "P_total_barg", "N_rps", "m_cat_g", "dp_um", "catalyst_state", "hold_s", "headspace_m3"], [
        {
            "case_id": case_id, "bath_C": f"{temp:.4f}", "P_total_barg": f"{gauge:.4f}",
            "N_rps": f"{speed:g}", "m_cat_g": f"{mass:.4f}", "dp_um": f"{diameter:g}",
            "catalyst_state": state, "hold_s": "600.0", "headspace_m3": "0.008000",
        }
        for case_id, temp, gauge, speed, mass, diameter, state in cases
    ])
    (data / "operator_note.txt").write_text(
        "Weekly note, 28 Apr 2024. Across the 12 rps and 9 rps sheet rates the uptake looked flat, "
        "so I treated the hydrogenation as kinetic and first order in hydrogen, same as the 2022 anisole charges. "
        "The sheet rates are what went into the report.\n"
    )
    fill_measurements(
        data, (k_ref, energy, kads_ref, enthalpy), VISIBLE_SEED, 0.011, "2.62000000e-06", slurry,
    )
    write_effectiveness(data, slurry, fresh, conditioned(16.0), 75.0, tic_of(75.0), 16.0, 8.0, 40.0, 3.0, "2024-04-27")
    for temp, gauge, speed, mass, diameter, _state in [(row[1], row[2], row[3], row[4], row[5], row[6]) for row in cases]:
        partial = gauge + p_atm - psat(temp + 273.15)
        if partial <= 0.05:
            raise SystemExit(f"visible hydrogen partial pressure {partial}")


def hidden_reference(tic: float) -> float:
    return tic - 1.8 - 0.00042 * (tic - 80.0) ** 2


def hidden_campaign() -> list[dict]:
    data = ROOT / "tests" / "hidden_campaign"
    fresh = 174.0
    rho = 1760.0
    slurry = 0.00125
    p_atm = 0.991
    antoine = (5.887036663360744, 2567.350921121411, 9.576306829370747)
    diffusivity = 7.5e-10
    k_ref, kads_ref, energy, enthalpy = 0.0134, 0.31, 64000.0, -39000.0
    cstar_ref, dh_sol = 1.62, 3100.0

    def c_one(temp_c: float) -> float:
        temp_k = temp_c + 273.15
        return cstar_ref * math.exp(-dh_sol / R_GAS * (1.0 / temp_k - 1.0 / T_REF))

    def kla(speed: float) -> float:
        return 0.031 * (speed / 10.0) ** 1.15

    def kc(speed: float, diameter: float) -> float:
        return 4.1e-5 * math.sqrt(speed / 10.0) * math.sqrt(300.0 / diameter)

    def conditioned(speed: float) -> float:
        return fresh / (1.0 + (speed / 14.0) ** 2)

    dump_json(data / "vessel.json", {
        "slurry_volume_m3": slurry,
        "site_atmospheric_pressure_bar": p_atm,
        "baffled": True,
        "working_volume_description": "1.25 L anisole slurry in the same autoclave, smaller fill than the toluene week",
    })
    dump_json(data / "solvent.json", {
        "name": "anisole",
        "antoine_form": "log10(P_sat_bar) = A - B / (T_K + C)",
        "A": antoine[0], "B": antoine[1], "C": antoine[2],
        "T_unit": "kelvin", "P_unit": "bar",
    })
    dump_json(data / "catalyst.json", {
        "metal": "palladium",
        "support": "activated carbon",
        "metal_location": "uniform",
        "eggshell_depth_um": None,
        "effective_diffusivity_m2_per_s": diffusivity,
        "particle_density_kg_m3": rho,
        "fresh_co_uptake_umol_per_g": fresh,
    })
    speeds = (3.0, 8.0, 14.0, 20.0)
    diameters = (110.0, 230.0, 250.0, 310.0, 390.0, 420.0, 470.0)
    write_csv(data / "kla_measurements.csv", ["N_rps", "kla_per_s", "note"], [
        {"N_rps": f"{speed:g}", "kla_per_s": f"{kla(speed):.8e}", "note": ""} for speed in speeds
    ] + [{"N_rps": "", "kla_per_s": "1.80000000e-02", "note": "tachometer blank"}])
    kc_rows = []
    for speed in speeds:
        for diameter in diameters:
            kc_rows.append({"N_rps": f"{speed:g}", "dp_um": f"{diameter:g}", "kc_m_per_s": f"{kc(speed, diameter):.8e}", "note": ""})
    kc_rows.append({"N_rps": "8", "dp_um": "250", "kc_m_per_s": f"{kc(8, 250) * 1.7:.8e}", "note": "baffles removed"})
    write_csv(data / "kc_measurements.csv", ["N_rps", "dp_um", "kc_m_per_s", "note"], kc_rows)
    write_csv(data / "catalyst_sites.csv", ["catalyst_state", "N_rps", "co_uptake_umol_per_g"], [
        {"catalyst_state": "conditioned", "N_rps": f"{speed:g}", "co_uptake_umol_per_g": f"{conditioned(speed):.6f}"}
        for speed in speeds
    ] + [
        {"catalyst_state": "fresh", "N_rps": "", "co_uptake_umol_per_g": f"{fresh:.6f}"},
        {"catalyst_state": "regenerated", "N_rps": "", "co_uptake_umol_per_g": "104.000000"},
    ])
    tics = list(range(25, 170, 5))
    write_csv(data / "thermocouple_calibration.csv", ["tic201_C", "reference_C"], [
        {"tic201_C": f"{tic:.4f}", "reference_C": f"{hidden_reference(tic):.6f}"} for tic in tics
    ])
    write_csv(data / "h2_solubility.csv", ["reference_C", "cstar_mol_m3_at_1bar"], [
        {"reference_C": f"{temp:.4f}", "cstar_mol_m3_at_1bar": f"{c_one(float(temp)):.8e}"} for temp in range(30, 161)
    ])
    factors = [
        ("2025-06-01", "2025-06-10", "3.15000000e-06"),
        ("2025-06-11", "2025-06-20", "2.40000000e-06"),
        ("2025-06-21", "2025-06-30", "1.88000000e-06"),
    ]
    write_csv(data / "gc_standards.csv", ["start_date", "end_date", "response_factor_mol_per_count"], [
        {"start_date": start, "end_date": end, "response_factor_mol_per_count": factor} for start, end, factor in factors
    ])
    train = [
        (45, 0.9, 8, 18, 250, "2025-06-04"),
        (45, 2.2, 8, 18, 250, "2025-06-12"),
        (45, 5.0, 8, 18, 250, "2025-06-22"),
        (45, 9.5, 8, 18, 250, "2025-06-08"),
        (75, 0.9, 8, 18, 250, "2025-06-14"),
        (75, 2.2, 8, 18, 250, "2025-06-24"),
        (75, 5.0, 8, 18, 250, "2025-06-06"),
        (75, 9.5, 8, 18, 250, "2025-06-16"),
        (105, 0.9, 8, 18, 250, "2025-06-26"),
        (105, 2.2, 8, 18, 250, "2025-06-03"),
        (105, 5.0, 8, 18, 250, "2025-06-13"),
        (105, 9.5, 8, 18, 250, "2025-06-23"),
        (130, 2.2, 8, 18, 250, "2025-06-09"),
        (130, 5.0, 8, 18, 250, "2025-06-18"),
        (130, 9.5, 8, 18, 250, "2025-06-27"),
        (155, 2.4, 8, 18, 250, "2025-06-05"),
        (155, 6.0, 8, 18, 250, "2025-06-15"),
        (155, 11.0, 8, 18, 250, "2025-06-25"),
        (105, 4.0, 8, 5, 250, "2025-06-02"),
        (105, 4.0, 8, 16, 250, "2025-06-11"),
        (105, 4.0, 8, 40, 250, "2025-06-21"),
        (105, 4.0, 8, 85, 250, "2025-06-07"),
        (105, 4.0, 3, 18, 250, "2025-06-17"),
        (105, 4.0, 14, 18, 250, "2025-06-28"),
        (105, 4.0, 20, 18, 250, "2025-06-04"),
        (105, 4.0, 14, 18, 110, "2025-06-19"),
        (105, 4.0, 14, 18, 420, "2025-06-29"),
        (45, 1.4, 20, 6, 110, "2025-06-10"),
        (155, 8.0, 3, 70, 420, "2025-06-20"),
        (130, 1.6, 20, 9, 110, "2025-06-30"),
        (75, 7.5, 3, 55, 420, "2025-06-01"),
        (155, 1.8, 20, 11, 110, "2025-06-12"),
    ]
    run_rows = []
    for index, (tic, gauge, speed, mass, diameter, date) in enumerate(train, start=1):
        run_rows.append({
            "run_id": f"H{index:02d}", "run_date": date, "status": "ok",
            "tic201_C": f"{tic:.4f}", "P_total_barg": f"{gauge:.4f}",
            "N_rps": f"{speed:g}", "m_cat_g": f"{mass:.4f}", "dp_um": f"{diameter:g}",
            "h2_peak_area": "1.0", "window_s": "600.0", "sheet_rate_mol_m3_s": "0.0",
        })
    aborted = [
        ("X1", "2025-06-08", 75, 4.0, 8, 18, 250, "8.40000000e+06"),
        ("X2", "2025-06-18", 105, 4.0, 14, 18, 250, "-1.20000000e+03"),
        ("X3", "2025-06-25", 130, 5.0, 8, 40, 250, "2.75000000e+07"),
        ("X4", "2025-06-14", 45, 2.2, 3, 18, 420, "4.40000000e+05"),
    ]
    for run_id, date, tic, gauge, speed, mass, diameter, area in aborted:
        run_rows.append({
            "run_id": run_id, "run_date": date, "status": "aborted",
            "tic201_C": f"{tic:.4f}", "P_total_barg": f"{gauge:.4f}",
            "N_rps": f"{speed:g}", "m_cat_g": f"{mass:.4f}", "dp_um": f"{diameter:g}",
            "h2_peak_area": area, "window_s": "600.0", "sheet_rate_mol_m3_s": "4.00000000e+00",
        })
    write_csv(data / "rate_runs.csv", list(run_rows[0].keys()), run_rows)

    cases = []
    serial = 1
    for temp in (41.0, 67.0, 93.0, 119.0, 137.0, 148.0):
        for gauge in (1.1, 3.4, 7.2, 12.5):
            for speed in speeds:
                for diameter in (230.0, 310.0, 390.0, 470.0):
                    for mass, state in ((2.4, "fresh"), (26.0, "conditioned"), (78.0, "regenerated"), (128.0, "conditioned")):
                        # This one cold, slow, lightly loaded charge sits on a flat surface-concentration branch.
                        if temp == 41.0 and gauge == 12.5 and speed == 3.0 and mass == 2.4 and diameter == 310.0 and state == "fresh":
                            continue
                        cases.append({
                            "case_id": f"C{serial:04d}",
                            "bath_C": f"{temp:.4f}",
                            "P_total_barg": f"{gauge:.4f}",
                            "N_rps": f"{speed:g}",
                            "m_cat_g": f"{mass:.4f}",
                            "dp_um": f"{diameter:g}",
                            "catalyst_state": state,
                            "hold_s": "4500.0",
                            "headspace_m3": "0.000330",
                        })
                        serial += 1
    edges = [
        ("E01", 153.0, 0.20, 8, 14.0, 310.0, "fresh", 4200.0, 0.000180),
        ("E02", 34.0, 13.5, 3, 1.6, 470.0, "conditioned", 6000.0, 0.000400),
        ("E03", 142.0, 11.0, 20, 150.0, 230.0, "regenerated", 3600.0, 0.000120),
        ("E04", 101.0, 0.55, 3, 120.0, 470.0, "fresh", 5400.0, 0.002000),
    ]
    for case_id, temp, gauge, speed, mass, diameter, state, hold, head in edges:
        cases.append({
            "case_id": case_id, "bath_C": f"{temp:.4f}", "P_total_barg": f"{gauge:.4f}",
            "N_rps": f"{speed:g}", "m_cat_g": f"{mass:.4f}", "dp_um": f"{diameter:g}",
            "catalyst_state": state, "hold_s": f"{hold:.4f}", "headspace_m3": f"{head:.6f}",
        })
    write_csv(
        data / "predict_cases.csv",
        ["case_id", "bath_C", "P_total_barg", "N_rps", "m_cat_g", "dp_um", "catalyst_state", "hold_s", "headspace_m3"],
        cases,
    )
    (data / "operator_note.txt").write_text(
        "Week of 30 Jun 2025. The anisole charges were logged from the spreadsheet rate. "
        "Internal gradients were ignored because the toluene week had matched a crushed aliquot, "
        "and the hold was treated as the initial rate times the clock time.\n"
    )
    fill_measurements(data, (k_ref, energy, kads_ref, enthalpy), HIDDEN_SEED, 0.012, "3.15000000e-06", slurry)
    write_effectiveness(
        data, slurry, fresh, conditioned(14.0), hidden_reference(105.0), 105.0, 14.0, 18.0, 420.0, 4.0, "2025-06-18",
    )
    return cases


def write_effectiveness(data, slurry, fresh, uptake, temperature_c, tic, speed, mass, diameter, gauge, date):
    campaign = forecast.load_campaign(data)
    rows = []
    for label, size, note in (
        ("whole", diameter, "side aliquot; crushed sample was mortared"),
        ("crushed", 12.0, "side aliquot; crushed sample was mortared"),
    ):
        loading = (mass / 1000.0) / slurry
        activity = uptake / fresh
        temperature = temperature_c + 273.15
        partial = gauge + campaign["p_atm"] - float(forecast.vapor_pressure(campaign["antoine"], temperature))
        c_one = float(forecast.solubility_1bar(campaign, np.array([temperature_c]))[0])
        kla = campaign["kla"][round(speed, 6)]
        # Crushed diameter is not a campaign sieve cut; scale the measured coefficient by the usual size dependence.
        base_diameter = diameter
        kc_value = campaign["kc"][(round(speed, 6), round(base_diameter, 6))] * math.sqrt(base_diameter / size)
        ks = kc_value * 6.0 * loading / (campaign["rho"] * size * 1e-6)
        k_value = forecast.k_of(np.array([temperature]), 0.0134 if campaign["distribution"] == "uniform" else 0.0125, 64000.0 if campaign["distribution"] == "uniform" else 62000.0)
        k_ads = forecast.kads_of(np.array([temperature]), 0.31 if campaign["distribution"] == "uniform" else 0.28, -39000.0 if campaign["distribution"] == "uniform" else -42000.0)
        shell = None if campaign["distribution"] == "uniform" else campaign["shell_um"] * 1e-6
        rate, *_ = forecast.initial_rates(
            np.array([partial]), np.array([temperature]), np.array([kla]), np.array([ks]),
            np.array([loading]), np.array([activity]), k_value, k_ads, np.array([c_one]),
            size * 0.5e-6, campaign["diffusivity"], campaign["rho"], shell, 36,
        )
        factor = forecast.response_factor(campaign, date)
        area = float(rate[0]) * 600.0 * slurry / factor
        rows.append({
            "aliquot": label, "run_date": date, "tic201_C": f"{tic:.6f}",
            "P_total_barg": f"{gauge:.4f}", "N_rps": f"{speed:g}", "m_cat_g": f"{mass:.4f}",
            "dp_um": f"{size:g}", "co_uptake_umol_per_g": f"{uptake:.6f}",
            "h2_peak_area": f"{area:.8e}", "window_s": "600.0", "note": note,
        })
    write_csv(data / "effectiveness_check.csv", list(rows[0].keys()), rows)


def fill_measurements(data: Path, params, seed: int, sigma: float, sheet_factor: str, slurry: float) -> None:
    campaign = forecast.load_campaign(data)
    bundle = forecast.observed_runs(campaign)
    latent = forecast.rates_for_params(campaign, bundle, *params, TRUTH_NODES)
    rng = np.random.default_rng(seed)
    noisy = latent * (1.0 + rng.normal(0.0, sigma, size=latent.size))
    noisy = np.maximum(noisy, 1e-8)
    path = data / "rate_runs.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
        fieldnames = list(rows[0].keys())
    ok_index = 0
    sheet = float(sheet_factor)
    for row in rows:
        if row["status"].strip() != "ok":
            continue
        factor = forecast.response_factor(campaign, row["run_date"])
        window = float(row["window_s"])
        area = noisy[ok_index] * window * slurry / factor
        row["h2_peak_area"] = f"{area:.8e}"
        row["sheet_rate_mol_m3_s"] = f"{(area * sheet / (window * slurry)):.8e}"
        ok_index += 1
    if ok_index != latent.size:
        raise SystemExit("ok-row count did not match the latent rates")
    write_csv(path, fieldnames, rows)


def assert_shift() -> None:
    visible = list(csv.DictReader((ROOT / "environment" / "data" / "predict_cases.csv").open()))
    hidden = list(csv.DictReader((ROOT / "tests" / "hidden_campaign" / "predict_cases.csv").open()))
    visible_keys = {
        (round(float(row["bath_C"]), 4), round(float(row["P_total_barg"]), 4), round(float(row["N_rps"]), 4),
         round(float(row["m_cat_g"]), 4), round(float(row["dp_um"]), 4), row["catalyst_state"])
        for row in visible
    }
    hidden_keys = set()
    for row in hidden:
        key = (round(float(row["bath_C"]), 4), round(float(row["P_total_barg"]), 4), round(float(row["N_rps"]), 4),
               round(float(row["m_cat_g"]), 4), round(float(row["dp_um"]), 4), row["catalyst_state"])
        if key in hidden_keys or key in visible_keys:
            raise SystemExit(f"repeated condition {key}")
        hidden_keys.add(key)
    visible_solvent = json.loads((ROOT / "environment" / "data" / "solvent.json").read_text())["name"]
    hidden_solvent = json.loads((ROOT / "tests" / "hidden_campaign" / "solvent.json").read_text())["name"]
    visible_metal = json.loads((ROOT / "environment" / "data" / "catalyst.json").read_text())["metal_location"]
    hidden_metal = json.loads((ROOT / "tests" / "hidden_campaign" / "catalyst.json").read_text())["metal_location"]
    if visible_solvent == hidden_solvent or visible_metal == hidden_metal:
        raise SystemExit("hidden campaign repeats the visible solvent or metal location")
    if not any(row["catalyst_state"] == "regenerated" for row in hidden):
        raise SystemExit("hidden campaign has no regenerated charge")
    hidden_runs = list(csv.DictReader((ROOT / "tests" / "hidden_campaign" / "rate_runs.csv").open()))
    if not any(row["status"] != "ok" for row in hidden_runs):
        raise SystemExit("hidden campaign has no aborted row")
    if any(row["dp_um"] in {"40", "40.0", "80", "80.0", "150", "150.0"} for row in hidden):
        raise SystemExit("hidden case uses a visible sieve cut")
    print(f"hidden cases {len(hidden)} visible cases {len(visible)} overlap 0")


def write_truth() -> None:
    campaign = forecast.load_campaign(ROOT / "tests" / "hidden_campaign")
    rows = forecast.load_cases(ROOT / "tests" / "hidden_campaign" / "predict_cases.csv")
    params = (0.0134, 64000.0, 0.31, -39000.0)
    print(f"latent truth for {len(rows)} cases at {TRUTH_NODES} nodes, {TRUTH_STEPS} steps")
    packed = forecast.forecast_cases(campaign, rows, params, nodes=TRUTH_NODES, steps=TRUTH_STEPS)
    # Spot-check a finer grid so the stored truth is not a coarse truncation.
    sample = [rows[0], rows[len(rows) // 2], rows[-4], rows[-1]]
    finer = forecast.forecast_cases(campaign, sample, params, nodes=112, steps=320)
    coarser = forecast.forecast_cases(campaign, sample, params, nodes=TRUTH_NODES, steps=TRUTH_STEPS)
    worst_rate = 0.0
    worst_consumed = 0.0
    for left, right in zip(coarser, finer):
        worst_rate = max(worst_rate, abs(left["rate_mol_m3_s"] - right["rate_mol_m3_s"]) / right["rate_mol_m3_s"])
        worst_consumed = max(worst_consumed, abs(left["hydrogen_consumed_mol_m3"] - right["hydrogen_consumed_mol_m3"]) / right["hydrogen_consumed_mol_m3"])
    print(f"truth vs finer grid max rel rate {worst_rate:.6f} consumed {worst_consumed:.6f}")
    if worst_rate > 0.008 or worst_consumed > 0.008:
        raise SystemExit("truth grid is not converged")
    for row in packed:
        fractions = (
            row["resistance_fraction_gas_liquid"]
            + row["resistance_fraction_liquid_solid"]
            + row["resistance_fraction_kinetic"]
        )
        if abs(fractions - 1.0) > 1e-8 or row["rate_mol_m3_s"] <= 0.0 or row["hydrogen_consumed_mol_m3"] <= 0.0:
            raise SystemExit(f"bad truth row {row['case_id']}")
    dump_json(ROOT / "tests" / "truth" / "expected.json", {"cases": packed})
    partials = []
    for row in rows:
        temperature = float(row["bath_C"]) + 273.15
        partial = float(row["P_total_barg"]) + campaign["p_atm"] - float(forecast.vapor_pressure(campaign["antoine"], temperature))
        partials.append(partial)
        if partial <= 0.04:
            raise SystemExit(f"{row['case_id']} hydrogen partial pressure {partial}")
    print(f"hydrogen partial pressure min {min(partials):.4f} max {max(partials):.4f}")


def write_provenance(hidden_cases: int) -> None:
    provenance = ROOT / "authoring" / "provenance"
    dump_json(provenance / "latent.json", {
        "visible": {
            "solvent": "toluene",
            "metal_location": "eggshell",
            "k_ref_mol_kg_s": 0.0125,
            "K_ref_m3_mol": 0.28,
            "E_J_mol": 62000.0,
            "dH_J_mol": -42000.0,
            "T_ref_K": T_REF,
            "noise_seed": VISIBLE_SEED,
            "noise_sigma_relative": 0.011,
        },
        "hidden": {
            "solvent": "anisole",
            "metal_location": "uniform",
            "effective_diffusivity_m2_per_s": 7.5e-10,
            "k_ref_mol_kg_s": 0.0134,
            "K_ref_m3_mol": 0.31,
            "E_J_mol": 64000.0,
            "dH_J_mol": -39000.0,
            "T_ref_K": T_REF,
            "truth_nodes": TRUTH_NODES,
            "truth_steps": TRUTH_STEPS,
            "noise_seed": HIDDEN_SEED,
            "noise_sigma_relative": 0.012,
            "case_count": hidden_cases,
            "shift": [
                "uniform palladium through large particles, so internal diffusion changes the rate",
                "anisole instead of toluene, with its own Antoine constants and solubility table",
                "temperatures, stirrer speeds, diameters, and catalyst masses outside the visible window",
                "regenerated catalyst state, absent from the visible sites table",
                "aborted chromatograph rows, a third response-factor window, and a curved thermocouple table",
                "near-boiling hydrogen partial pressure and holds long enough to empty a small headspace",
            ],
        },
    })
    dump_json(provenance / "seeds.json", {
        "visible_noise_seed": VISIBLE_SEED,
        "hidden_noise_seed": HIDDEN_SEED,
        "hidden_is_not_a_reseed": True,
        "hidden_design": "explicit factorial plus four edge charges in generate_campaign.py hidden_campaign()",
        "truth": "latent parameters on the fine finite-volume grid, not the reference least-squares fit",
    })
    paths = [
        ROOT / "environment" / "data_note.md",
        *sorted((ROOT / "environment" / "data").iterdir()),
        *sorted((ROOT / "tests" / "hidden_campaign").iterdir()),
        ROOT / "tests" / "truth" / "expected.json",
    ]
    lines = []
    for path in paths:
        if path.is_file():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            lines.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
    (provenance / "checksums.sha256").write_text("\n".join(lines) + "\n")


def main() -> None:
    visible_campaign()
    cases = hidden_campaign()
    assert_shift()
    write_truth()
    write_provenance(len(cases))
    print("wrote campaigns", len(cases))


if __name__ == "__main__":
    main()
