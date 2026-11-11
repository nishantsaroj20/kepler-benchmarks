#!/usr/bin/env python3
"""Seeded campaign generator.

Ground truth is this latent forward model. The reference solution is not
consulted. Noise, the blank tachometer row, and the baffles-removed dissolution
row are part of the shipped records; the latent rate has neither.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

SEED = 20261111
R_GAS = 8.314462618
P_ATM = 1.01325  # bar
V_SLURRY = 0.0015  # m3
RHO_P = 1900.0  # kg/m3
T_REF = 348.15  # K
WINDOW_S = 600.0
FRESH_CO = 186.0  # umol/g
EGGSHELL_UM = 3.0

ANTOINE_A = 4.07827
ANTOINE_B = 1343.943
ANTOINE_C = -53.773

K_REF = 0.0125  # mol/kg/s
KADS_REF = 0.28  # m3/mol
E_K = 62000.0  # J/mol
DH_ADS = -42000.0  # J/mol
CSTAR_REF = 2.85  # mol/m3 at 1 bar and T_REF
DH_SOL = 4800.0  # J/mol
KLA0 = 0.048
KLA_EXP = 0.82
KC0 = 6.8e-5
N_CHAR = 21.0
RATE_SIGMA = 0.011

RF_BEFORE = 2.62e-6
RF_AFTER = 2.02e-6
MAINT = "2024-04-16"

REL_TOL = 0.040
ABS_TOL = 0.0006
FRAC_TOL = 0.030

SPEEDS = (6.0, 9.0, 12.0, 16.0)
DIAMETERS = (40.0, 80.0, 150.0)

RUNS = [
    (60, 0.8, 12, 22, 80, "2024-04-18"),
    (60, 1.6, 12, 22, 80, "2024-04-09"),
    (60, 3.0, 12, 22, 80, "2024-04-19"),
    (60, 5.0, 12, 22, 80, "2024-04-08"),
    (60, 7.0, 12, 22, 80, "2024-04-11"),
    (75, 0.8, 12, 22, 80, "2024-04-17"),
    (75, 1.6, 12, 22, 80, "2024-04-10"),
    (75, 3.0, 12, 22, 80, "2024-04-18"),
    (75, 5.0, 12, 22, 80, "2024-04-07"),
    (75, 7.0, 12, 22, 80, "2024-04-12"),
    (90, 0.8, 12, 22, 80, "2024-04-20"),
    (90, 1.6, 12, 22, 80, "2024-04-08"),
    (90, 3.0, 12, 22, 80, "2024-04-22"),
    (90, 5.0, 12, 22, 80, "2024-04-09"),
    (90, 7.0, 12, 22, 80, "2024-04-14"),
    (75, 3.0, 12, 8, 80, "2024-04-11"),
    (75, 3.0, 12, 14, 80, "2024-04-23"),
    (75, 3.0, 12, 35, 80, "2024-04-06"),
    (75, 3.0, 12, 50, 80, "2024-04-24"),
    (75, 3.0, 6, 22, 80, "2024-04-13"),
    (75, 3.0, 9, 22, 80, "2024-04-16"),
    (75, 3.0, 16, 22, 80, "2024-04-05"),
    (75, 3.0, 12, 22, 40, "2024-04-15"),
    (75, 3.0, 12, 22, 150, "2024-04-25"),
    (75, 3.0, 6, 8, 40, "2024-04-04"),
    (75, 3.0, 16, 8, 40, "2024-04-26"),
    (75, 2.0, 6, 50, 150, "2024-04-03"),
    (60, 3.0, 6, 40, 150, "2024-04-27"),
    (90, 2.0, 16, 10, 40, "2024-04-02"),
    (75, 3.0, 12, 22, 80, "2024-04-28"),
]

CASES = [
    ("H1", 75.0, 11.5, 12.0, 40.0, 40.0, "conditioned"),
    ("L1", 75.0, 0.25, 12.0, 22.0, 80.0, "conditioned"),
    ("T1", 114.0, 1.55, 12.0, 10.0, 80.0, "conditioned"),
    ("W1", 75.0, 2.0, 6.0, 85.0, 80.0, "conditioned"),
    ("D1", 60.0, 5.5, 9.0, 36.0, 150.0, "conditioned"),
    ("F1", 75.0, 4.5, 16.0, 10.0, 40.0, "fresh"),
    ("X1", 108.0, 6.5, 6.0, 78.0, 150.0, "fresh"),
    ("K1", 60.0, 0.85, 16.0, 7.0, 40.0, "conditioned"),
]


def psat_bar(T_K: float) -> float:
    return 10 ** (ANTOINE_A - ANTOINE_B / (T_K + ANTOINE_C))


def k_of_T(T_K: float) -> float:
    return K_REF * math.exp(-E_K / R_GAS * (1.0 / T_K - 1.0 / T_REF))


def Kads_of_T(T_K: float) -> float:
    return KADS_REF * math.exp(-DH_ADS / R_GAS * (1.0 / T_K - 1.0 / T_REF))


def cstar_1bar(T_K: float) -> float:
    return CSTAR_REF * math.exp(-DH_SOL / R_GAS * (1.0 / T_K - 1.0 / T_REF))


def kla_of(N: float) -> float:
    return KLA0 * (N / 12.0) ** KLA_EXP


def kc_of(N: float, dp_um: float) -> float:
    return KC0 * math.sqrt(N / 12.0) * math.sqrt(80.0 / dp_um)


def activity(N: float, state: str) -> float:
    if state == "fresh":
        return 1.0
    return 1.0 / (1.0 + (N / N_CHAR) ** 2)


def tic_from_ref(T_C: float) -> float:
    return T_C + 4.6 + 0.018 * (T_C - 75.0)


def hydrogen_partial(T_K: float, P_barg: float) -> float:
    return P_barg + P_ATM - psat_bar(T_K)


def ks_vol(N: float, dp_um: float, w: float) -> float:
    area = 6.0 * w / (RHO_P * (dp_um * 1e-6))
    return kc_of(N, dp_um) * area


def rate_closed(Cstar: float, kLa: float, ks: float, w: float, a: float, k: float, K: float) -> float:
    Rmt = 1.0 / kLa + 1.0 / ks
    B = w * a * k
    Aq = K * Rmt
    Bq = -(1.0 + K * Cstar + B * K * Rmt)
    Cq = B * K * Cstar
    disc = Bq * Bq - 4.0 * Aq * Cq
    if disc < 0.0 and disc > -1e-18:
        disc = 0.0
    r = (2.0 * Cq) / (-Bq + math.sqrt(disc))
    return r


def evaluate(T_C: float, P_barg: float, N: float, m_g: float, dp_um: float, state: str):
    T_K = T_C + 273.15
    w = (m_g / 1000.0) / V_SLURRY
    a = activity(N, state)
    p_h2 = hydrogen_partial(T_K, P_barg)
    Cstar = cstar_1bar(T_K) * p_h2
    kLa = kla_of(N)
    ks = ks_vol(N, dp_um, w)
    r = rate_closed(Cstar, kLa, ks, w, a, k_of_T(T_K), Kads_of_T(T_K))
    d_gl = r / kLa
    d_ls = r / ks
    Cs = Cstar - d_gl - d_ls
    return {
        "rate_mol_m3_s": r,
        "resistance_fraction_gas_liquid": d_gl / Cstar,
        "resistance_fraction_liquid_solid": d_ls / Cstar,
        "resistance_fraction_kinetic": Cs / Cstar,
        "p_h2_bar": p_h2,
        "Cstar": Cstar,
        "Cs": Cs,
    }


def rf_for(date: str) -> float:
    return RF_BEFORE if date < MAINT else RF_AFTER


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def data_note() -> str:
    return """\
Files are UTF-8. Each CSV file has one header row and no index column. Each JSON file is one object. Dates are calendar dates in the form YYYY-MM-DD. Stirrer speeds are revolutions per second. Particle diameters are micrometers. A pressure written as gauge is bar gauge. The site atmospheric pressure in bar, the slurry volume in cubic meters, and the baffles flag are in `/app/data/vessel.json`.

`/app/data/rate_runs.csv` is the hydrogenation campaign. `tic201_C` is thermowell TIC-201 in degrees Celsius. `P_total_barg` is the headspace gauge pressure. `m_cat_g` is the dry catalyst charge in grams. `h2_peak_area` is the chromatograph area for the hydrogen consumed during `window_s` seconds. `sheet_rate_mol_m3_s` is the number exported by the acquisition spreadsheet.

Moles of hydrogen equal `h2_peak_area` times `response_factor_mol_per_count` from `/app/data/gc_standards.csv`. Use the standards row whose `start_date` and `end_date` both include `run_date`. The slurry-volume-specific rate, in mol m$ ^{-3} $ s$ ^{-1} $, is those moles divided by `window_s` and by the slurry volume. The substrate is in large excess and each window stays at low conversion, so this rate is the initial hydrogen rate. Catalyst mass concentration w, in kg m$ ^{-3} $, is the charge in kilograms divided by the slurry volume.

Absolute headspace pressure in bar is `P_total_barg` plus the site atmospheric pressure. `/app/data/solvent.json` gives the toluene vapor pressure through log10($ P_{sat} $ / 1 bar) = A - B / ($ T_K $ + C), with $ T_K $ the reference-probe temperature in kelvin and $ P_{sat} $ in bar. Hydrogen partial pressure in bar is the absolute headspace pressure minus that vapor pressure.

`/app/data/h2_solubility.csv` lists the dissolved-hydrogen concentration, in mol m$^{-3}$, at a hydrogen partial pressure of 1 bar, against reference-probe temperature in degrees Celsius. Between adjacent rows, interpolate the natural logarithm of that concentration linearly in reciprocal thermodynamic temperature. The interface concentration is the interpolated 1 bar value multiplied by the hydrogen partial pressure in bar.

`/app/data/thermocouple_calibration.csv` pairs `tic201_C` with `reference_C`. Interpolate `reference_C` linearly in `tic201_C`. Temperatures in `/app/data/predict_cases.csv` are already the reference probe, in degrees Celsius. Every historical rate was measured after the catalyst had been conditioned for 4 h at that run's stirrer speed.

`/app/data/kla_measurements.csv` gives the baffled-vessel gas-liquid coefficient `kla_per_s`, in s$ ^{-1} $. `/app/data/kc_measurements.csv` gives the liquid-solid coefficient `kc_m_per_s`, in m s$ ^{-1} $. A row with a blank stirrer speed is not usable. A row whose note is `baffles removed` is not the baffled campaign vessel. External catalyst area per slurry volume is $ 6 w / (\\rho_p d_p) $, where $ \\rho_p $ is the particle density in `/app/data/catalyst.json` and $ d_p $ is the particle diameter in meters. The liquid-solid volumetric coefficient is `kc_m_per_s` times that area.

The catalytic step scales in proportion to CO chemisorption uptake, in micromoles per gram. `/app/data/catalyst.json` gives the fresh uptake, the particle density, and the eggshell metal depth. `/app/data/catalyst_sites.csv` gives the uptake after conditioning at each stirrer speed. A charge with `catalyst_state` equal to `conditioned` uses the uptake at its own stirrer speed. A charge with `catalyst_state` equal to `fresh` uses the fresh uptake.

`/app/data/effectiveness_check.csv` is a crushed-versus-whole side test on one conditioned aliquot. It is not a row of `/app/data/rate_runs.csv`. `/app/data/operator_note.txt` is the operator's weekly note.
"""


def assert_baseline_gaps(expected: list[dict]) -> None:
    rates = np.array([row["rate_mol_m3_s"] for row in expected], dtype=float)
    fracs = {
        key: np.array([row[key] for row in expected], dtype=float)
        for key in (
            "resistance_fraction_gas_liquid",
            "resistance_fraction_liquid_solid",
            "resistance_fraction_kinetic",
        )
    }
    mean_r = float(rates.mean())
    med_r = float(np.median(rates))

    def rate_tol(value: float) -> float:
        return max(REL_TOL * abs(value), ABS_TOL)

    for row, rate in zip(expected, rates):
        if abs(rate - mean_r) <= rate_tol(rate):
            raise SystemExit(f"mean rate lands in {row['case_id']}")
        if abs(rate - med_r) <= rate_tol(rate):
            raise SystemExit(f"median rate lands in {row['case_id']}")
        if abs(rate - 0.0) <= rate_tol(rate) or abs(rate - 1.0) <= rate_tol(rate):
            raise SystemExit(f"constant 0 or 1 lands in {row['case_id']}")
        for key, values in fracs.items():
            if values.min() <= FRAC_TOL or values.max() >= 1.0 - FRAC_TOL:
                raise SystemExit(f"fraction {key} approaches 0 or 1")
            if abs(sum(row[k] for k in fracs) - 1.0) > 1e-12:
                raise SystemExit(f"fractions do not sum to 1 for {row['case_id']}")

    for i, left in enumerate(rates):
        for right in rates[i + 1 :]:
            if abs(left - right) <= max(rate_tol(left), rate_tol(right)):
                raise SystemExit("two truth rates sit inside one tolerance")

    inputs = {
        "bath_C": np.array([c[1] for c in CASES], dtype=float),
        "P_total_barg": np.array([c[2] for c in CASES], dtype=float),
        "N_rps": np.array([c[3] for c in CASES], dtype=float),
        "m_cat_g": np.array([c[4] for c in CASES], dtype=float),
        "dp_um": np.array([c[5] for c in CASES], dtype=float),
    }
    for name, column in inputs.items():
        for rate, copied in zip(rates, column):
            if abs(rate - copied) <= rate_tol(rate):
                raise SystemExit(f"copying {name} lands on a rate")


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    data_dir = root / "environment" / "data"
    truth_path = root / "tests" / "truth" / "expected.json"
    prov = Path(__file__).resolve().parent
    data_dir.mkdir(parents=True, exist_ok=True)
    truth_path.parent.mkdir(parents=True, exist_ok=True)

    for T_C, P, N, m, dp, date in RUNS:
        ev = evaluate(T_C, P, N, m, dp, "conditioned")
        if ev["p_h2_bar"] <= 0.25 or ev["Cs"] <= 0.2:
            raise SystemExit(f"training point not physical: {T_C} {P} {N} {m} {dp}")
    expected = []
    for case in CASES:
        cid, T_C, P, N, m, dp, state = case
        ev = evaluate(T_C, P, N, m, dp, state)
        if ev["p_h2_bar"] <= 0.25 or ev["Cs"] <= 0.05:
            raise SystemExit(f"case {cid} not physical: {ev}")
        expected.append(
            {
                "case_id": cid,
                "rate_mol_m3_s": ev["rate_mol_m3_s"],
                "resistance_fraction_gas_liquid": ev["resistance_fraction_gas_liquid"],
                "resistance_fraction_liquid_solid": ev["resistance_fraction_liquid_solid"],
                "resistance_fraction_kinetic": ev["resistance_fraction_kinetic"],
            }
        )
    assert_baseline_gaps(expected)

    rng = np.random.default_rng(SEED)
    rate_rows = []
    for index, (T_C, P, N, m, dp, date) in enumerate(RUNS, start=1):
        ev = evaluate(T_C, P, N, m, dp, "conditioned")
        observed = ev["rate_mol_m3_s"] * (1.0 + float(rng.normal(0.0, RATE_SIGMA)))
        factor = rf_for(date)
        area = observed * V_SLURRY * WINDOW_S / factor
        sheet = area * RF_BEFORE / (V_SLURRY * WINDOW_S)
        rate_rows.append(
            {
                "run_id": f"R{index:02d}",
                "run_date": date,
                "tic201_C": f"{tic_from_ref(T_C):.6f}",
                "P_total_barg": f"{P:.4f}",
                "N_rps": f"{N:.0f}",
                "m_cat_g": f"{m:.4f}",
                "dp_um": f"{dp:.1f}",
                "h2_peak_area": f"{area:.8e}",
                "window_s": f"{WINDOW_S:.1f}",
                "sheet_rate_mol_m3_s": f"{sheet:.8e}",
            }
        )

    kla_rows = []
    for speed in SPEEDS:
        kla_rows.append(
            {
                "N_rps": f"{speed:.0f}",
                "kla_per_s": f"{kla_of(speed):.8e}",
                "note": "",
            }
        )
    kla_rows.append({"N_rps": "", "kla_per_s": f"{0.041:.8e}", "note": "tachometer blank"})

    kc_rows = []
    for speed in SPEEDS:
        for diameter in DIAMETERS:
            kc_rows.append(
                {
                    "N_rps": f"{speed:.0f}",
                    "dp_um": f"{diameter:.1f}",
                    "kc_m_per_s": f"{kc_of(speed, diameter):.8e}",
                    "note": "",
                }
            )
    kc_rows.append(
        {
            "N_rps": "12",
            "dp_um": "80.0",
            "kc_m_per_s": f"{1.65 * kc_of(12.0, 80.0):.8e}",
            "note": "baffles removed",
        }
    )

    site_rows = []
    for speed in SPEEDS:
        site_rows.append(
            {
                "N_rps": f"{speed:.0f}",
                "co_uptake_umol_per_g": f"{FRESH_CO * activity(speed, 'conditioned'):.6f}",
            }
        )

    cal_rows = []
    for T_C in range(40, 131, 5):
        cal_rows.append(
            {
                "tic201_C": f"{tic_from_ref(float(T_C)):.6f}",
                "reference_C": f"{float(T_C):.4f}",
            }
        )

    sol_rows = []
    for T_C in range(40, 131, 2):
        sol_rows.append(
            {
                "reference_C": f"{float(T_C):.4f}",
                "cstar_mol_m3_at_1bar": f"{cstar_1bar(T_C + 273.15):.8e}",
            }
        )

    # Side test: whole 40 um versus crushed 8 um, kinetic corner, independent noise.
    side_rows = []
    for label, diameter in (("whole", 40.0), ("crushed", 8.0)):
        ev = evaluate(75.0, 3.0, 16.0, 8.0, diameter, "conditioned")
        observed = ev["rate_mol_m3_s"] * (1.0 + float(rng.normal(0.0, RATE_SIGMA)))
        factor = RF_AFTER
        area = observed * V_SLURRY * WINDOW_S / factor
        side_rows.append(
            {
                "aliquot": label,
                "run_date": "2024-04-27",
                "tic201_C": f"{tic_from_ref(75.0):.6f}",
                "P_total_barg": "3.0000",
                "N_rps": "16",
                "m_cat_g": "8.0000",
                "dp_um": f"{diameter:.1f}",
                "co_uptake_umol_per_g": f"{FRESH_CO * activity(16.0, 'conditioned'):.6f}",
                "h2_peak_area": f"{area:.8e}",
                "window_s": f"{WINDOW_S:.1f}",
                "note": "side aliquot conditioned at 16 rps; crushed sample was mortared",
            }
        )

    case_rows = []
    for cid, T_C, P, N, m, dp, state in CASES:
        case_rows.append(
            {
                "case_id": cid,
                "bath_C": f"{T_C:.4f}",
                "P_total_barg": f"{P:.4f}",
                "N_rps": f"{N:.0f}",
                "m_cat_g": f"{m:.4f}",
                "dp_um": f"{dp:.1f}",
                "catalyst_state": state,
            }
        )

    write_csv(
        data_dir / "rate_runs.csv",
        [
            "run_id",
            "run_date",
            "tic201_C",
            "P_total_barg",
            "N_rps",
            "m_cat_g",
            "dp_um",
            "h2_peak_area",
            "window_s",
            "sheet_rate_mol_m3_s",
        ],
        rate_rows,
    )
    write_csv(data_dir / "kla_measurements.csv", ["N_rps", "kla_per_s", "note"], kla_rows)
    write_csv(
        data_dir / "kc_measurements.csv",
        ["N_rps", "dp_um", "kc_m_per_s", "note"],
        kc_rows,
    )
    write_csv(
        data_dir / "catalyst_sites.csv",
        ["N_rps", "co_uptake_umol_per_g"],
        site_rows,
    )
    write_csv(
        data_dir / "thermocouple_calibration.csv",
        ["tic201_C", "reference_C"],
        cal_rows,
    )
    write_csv(
        data_dir / "h2_solubility.csv",
        ["reference_C", "cstar_mol_m3_at_1bar"],
        sol_rows,
    )
    write_csv(
        data_dir / "gc_standards.csv",
        ["start_date", "end_date", "response_factor_mol_per_count"],
        [
            {
                "start_date": "2024-04-01",
                "end_date": "2024-04-15",
                "response_factor_mol_per_count": f"{RF_BEFORE:.8e}",
            },
            {
                "start_date": "2024-04-16",
                "end_date": "2024-04-30",
                "response_factor_mol_per_count": f"{RF_AFTER:.8e}",
            },
        ],
    )
    write_csv(
        data_dir / "effectiveness_check.csv",
        [
            "aliquot",
            "run_date",
            "tic201_C",
            "P_total_barg",
            "N_rps",
            "m_cat_g",
            "dp_um",
            "co_uptake_umol_per_g",
            "h2_peak_area",
            "window_s",
            "note",
        ],
        side_rows,
    )
    write_csv(
        data_dir / "predict_cases.csv",
        ["case_id", "bath_C", "P_total_barg", "N_rps", "m_cat_g", "dp_um", "catalyst_state"],
        case_rows,
    )

    (data_dir / "vessel.json").write_text(
        json.dumps(
            {
                "slurry_volume_m3": V_SLURRY,
                "site_atmospheric_pressure_bar": P_ATM,
                "baffled": True,
                "working_volume_description": "1.5 L slurry, gas headspace above a single pitched-blade impeller",
            },
            indent=2,
        )
        + "\n"
    )
    (data_dir / "solvent.json").write_text(
        json.dumps(
            {
                "name": "toluene",
                "antoine_form": "log10(P_sat_bar) = A - B / (T_K + C)",
                "A": ANTOINE_A,
                "B": ANTOINE_B,
                "C": ANTOINE_C,
                "T_unit": "kelvin",
                "P_unit": "bar",
            },
            indent=2,
        )
        + "\n"
    )
    (data_dir / "catalyst.json").write_text(
        json.dumps(
            {
                "metal": "palladium",
                "support": "activated carbon",
                "metal_location": "eggshell",
                "eggshell_depth_um": EGGSHELL_UM,
                "particle_density_kg_m3": RHO_P,
                "fresh_co_uptake_umol_per_g": FRESH_CO,
            },
            indent=2,
        )
        + "\n"
    )
    (data_dir / "operator_note.txt").write_text(
        "Weekly note, 28 Apr 2024. Across the 12 rps and 9 rps sheet rates the "
        "uptake looked flat, so I treated the hydrogenation as kinetic and first "
        "order in hydrogen, same as the 2022 anisole charges. The sheet rates are "
        "what went into the report.\n"
    )
    note_path = root / "environment" / "data_note.md"
    note_path.write_text(data_note())

    payload = {"cases": expected}
    truth_path.write_text(json.dumps(payload, indent=2) + "\n")
    latent = {
        "seed": SEED,
        "rate_sigma": RATE_SIGMA,
        "gas_constant_J_per_mol_K": R_GAS,
        "T_ref_K": T_REF,
        "k_ref_mol_per_kg_s": K_REF,
        "K_ads_ref_m3_per_mol": KADS_REF,
        "E_J_per_mol": E_K,
        "delta_H_ads_J_per_mol": DH_ADS,
        "cstar_ref_mol_m3": CSTAR_REF,
        "delta_H_solution_J_per_mol": DH_SOL,
        "kla0_per_s": KLA0,
        "kla_exponent": KLA_EXP,
        "kc0_m_per_s": KC0,
        "N_char_rps": N_CHAR,
        "source": "Latent two-film slurry hydrogenation. Langmuir-Hinshelwood in surface hydrogen, Arrhenius surface coefficient, van't Hoff adsorption coefficient, activity proportional to CO uptake. Prediction truth is this forward model, not a fit to the noisy rates.",
        "cases": expected,
    }
    (prov / "latent.json").write_text(json.dumps(latent, indent=2) + "\n")
    (prov / "seeds.json").write_text(
        json.dumps(
            {
                "numpy_generator": "np.random.default_rng",
                "seed": SEED,
                "rate_sigma_relative": RATE_SIGMA,
                "draws": "one normal draw per rate run, then one per effectiveness aliquot, in file order",
            },
            indent=2,
        )
        + "\n"
    )

    checksum_lines = []
    files = sorted(data_dir.rglob("*")) + [note_path, truth_path]
    for path in files:
        if path.is_file():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            checksum_lines.append(f"{digest}  {path.relative_to(root)}")
    (prov / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n")
    print(f"wrote {len(rate_rows)} runs and {len(expected)} cases")
    for row in expected:
        print(
            f"  {row['case_id']} {row['rate_mol_m3_s']:.6f} "
            f"{row['resistance_fraction_gas_liquid']:.3f}/"
            f"{row['resistance_fraction_liquid_solid']:.3f}/"
            f"{row['resistance_fraction_kinetic']:.3f}"
        )


if __name__ == "__main__":
    main()
