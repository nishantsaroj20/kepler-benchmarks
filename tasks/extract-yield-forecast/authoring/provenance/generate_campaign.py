#!/usr/bin/env python3
"""Seed the extraction campaign and the latent forecast truth.

The truth is this forward model on the charge inputs. It is not the
least-squares fit used by the reference solution.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from datetime import date, timedelta
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "environment" / "data"
TRUTH = ROOT / "tests" / "truth" / "expected.json"
EVIDENCE = ROOT / "authoring" / "evidence" / "calibration_results.json"
SEEDS = ROOT / "authoring" / "provenance" / "seeds.json"
CHECKSUMS = ROOT / "authoring" / "provenance" / "checksums.sha256"
FORECAST = ROOT / "solution" / "forecast.py"

SEED = 20260318
PKA = 4.58
K_REF = 62.0
TEMPERATURE_SLOPE = 2450.0
KLA_REF = 0.055
REL_TOL = 0.025
ABS_RECOVERY = 0.006
ABS_ORGANIC = 0.20

CALIBRATION = {
    "solvent": "4-methyl-2-pentanone",
    "solute": "K-17 carboxylic acid",
    "temperature_reference_k": 298.15,
    "stirrer_reference_rps": 3.2,
    "stirrer_exponent": 0.85,
    "organic_assay_blank_mol_m3": 1.8,
    "thermocouple_celsius": {"a0": -5.5, "a1": 1.08, "a2": -0.0008},
    "ph_windows": [
        {
            "start": "2024-02-01",
            "end": "2024-05-18",
            "intercept": 7.06,
            "millivolt_offset": 4.0,
            "slope_mv_per_ph": 61.4,
            "sheet_slope_sign": -1,
        },
        {
            "start": "2024-05-19",
            "end": "2024-09-30",
            "intercept": 6.97,
            "millivolt_offset": -1.2,
            "slope_mv_per_ph": 58.6,
            "sheet_slope_sign": 1,
        },
    ],
}


def load_forecast():
    spec = importlib.util.spec_from_file_location("extract_forecast", FORECAST)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def broth_celsius(tic_c: float) -> float:
    coeff = CALIBRATION["thermocouple_celsius"]
    return coeff["a0"] + coeff["a1"] * tic_c + coeff["a2"] * tic_c * tic_c


def window_for(day: str) -> dict:
    for window in CALIBRATION["ph_windows"]:
        if window["start"] <= day <= window["end"]:
            return window
    raise KeyError(day)


def millivolts(ph: float, window: dict) -> float:
    return window["millivolt_offset"] + (window["intercept"] - ph) * window["slope_mv_per_ph"]


def sheet_ph(millivolt: float, window: dict) -> float:
    sign = window["sheet_slope_sign"]
    return window["intercept"] - sign * (millivolt - window["millivolt_offset"]) / window["slope_mv_per_ph"]


def recovery_of(ph, temperature_c, solvent_ratio, stirrer_rps, contact_s, stages) -> float:
    forecast = load_forecast()
    return forecast.recovery_of(
        ph, temperature_c, solvent_ratio, stirrer_rps, contact_s, stages,
        PKA, K_REF, TEMPERATURE_SLOPE, KLA_REF,
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def dated(start: date, index: int) -> str:
    return (start + timedelta(days=int(index))) .isoformat()


def build_runs(rng: np.random.Generator) -> list[dict]:
    rows = []
    early_start = date(2024, 2, 6)
    late_start = date(2024, 5, 20)
    for index in range(120):
        early = index < 54
        window = CALIBRATION["ph_windows"][0 if early else 1]
        day = dated(early_start, index if early else index - 54)
        if not early and day < "2024-05-19":
            day = dated(late_start, index - 54)
        ph = float(rng.uniform(3.25, 6.05))
        tic = float(rng.uniform(16.0, 42.0))
        solvent_ratio = float(rng.uniform(0.22, 1.15))
        stirrer = float(rng.uniform(1.5, 5.5))
        contact = float(rng.uniform(24.0, 200.0))
        stages = 2 if rng.random() < 0.34 else 1
        feed = float(rng.uniform(9.0, 38.0))
        temperature = broth_celsius(tic)
        recovery = recovery_of(ph, temperature, solvent_ratio, stirrer, contact, stages)
        status = "ok"
        if index in {7, 18, 29, 41, 58, 66, 77, 83, 96, 109}:
            status = "aborted"
        elif index in {11, 24, 37, 52, 71, 88, 102, 114}:
            status = "unmixed"
        if status == "ok":
            aqueous = feed * (1.0 - recovery) * (1.0 + 0.004 * float(rng.normal()))
            organic = recovery * feed / solvent_ratio + CALIBRATION["organic_assay_blank_mol_m3"]
            organic += 0.03 * float(rng.normal())
        elif status == "aborted":
            aqueous = feed * float(rng.uniform(0.15, 1.25))
            organic = float(rng.uniform(0.0, 20.0))
        else:
            aqueous = feed * (1.0 - 0.04 * recovery)
            organic = 0.05 * recovery * feed / solvent_ratio
        mv = millivolts(ph, window)
        rows.append({
            "case_id": f"R{index + 1:03d}",
            "date": day,
            "status": status,
            "tic_c": f"{tic:.4f}",
            "ph_mv": f"{mv:.6f}",
            "ph_sheet": f"{sheet_ph(mv, window):.6f}",
            "solvent_ratio": f"{solvent_ratio:.6f}",
            "stirrer_rps": f"{stirrer:.4f}",
            "contact_s": f"{contact:.2f}",
            "stage_count": str(stages),
            "c_feed_mol_m3": f"{feed:.6f}",
            "c_aq_mol_m3": f"{aqueous:.6f}",
            "c_org_assay_mol_m3": f"{organic:.6f}",
            "recovery_sheet": f"{(organic * solvent_ratio / feed):.6f}",
        })
    return rows


def build_charges(rng: np.random.Generator, runs: list[dict]) -> tuple[list[dict], list[dict]]:
    seen = {
        (
            row["date"][:7],
            row["tic_c"],
            row["solvent_ratio"],
            row["stirrer_rps"],
            row["contact_s"],
            row["stage_count"],
        )
        for row in runs
    }
    charges = []
    truth = []
    early_start = date(2024, 3, 2)
    late_start = date(2024, 6, 3)
    index = 0
    while len(charges) < 56:
        early = len(charges) < 24
        window = CALIBRATION["ph_windows"][0 if early else 1]
        day = dated(early_start if early else late_start, index)
        if early and day > "2024-05-18":
            day = "2024-05-18"
        if not early and day < "2024-05-19":
            day = "2024-05-19"
        ph = float(rng.uniform(3.1, 6.6))
        tic = float(rng.uniform(14.0, 50.0))
        solvent_ratio = float(rng.uniform(0.16, 1.85))
        stirrer = float(rng.uniform(1.2, 7.4))
        contact = float(rng.uniform(8.0, 280.0))
        stages = 2 if rng.random() < 0.4 else 1
        feed = float(rng.uniform(8.0, 42.0))
        key = (day[:7], f"{tic:.4f}", f"{solvent_ratio:.6f}", f"{stirrer:.4f}", f"{contact:.2f}", str(stages))
        index += 1
        if key in seen:
            continue
        seen.add(key)
        temperature = broth_celsius(tic)
        recovery = recovery_of(ph, temperature, solvent_ratio, stirrer, contact, stages)
        organic = recovery * feed / solvent_ratio
        mv = millivolts(ph, window)
        case_id = f"C{len(charges) + 1:03d}"
        charges.append({
            "case_id": case_id,
            "date": day,
            "status": "ok",
            "tic_c": f"{tic:.4f}",
            "ph_mv": f"{mv:.6f}",
            "ph_sheet": f"{sheet_ph(mv, window):.6f}",
            "solvent_ratio": f"{solvent_ratio:.6f}",
            "stirrer_rps": f"{stirrer:.4f}",
            "contact_s": f"{contact:.2f}",
            "stage_count": str(stages),
            "c_feed_mol_m3": f"{feed:.6f}",
        })
        truth.append({
            "case_id": case_id,
            "recovery": recovery,
            "organic_mol_m3": organic,
        })
    return charges, truth


def outside(predicted: float, truth_value: float, absolute: float) -> bool:
    return abs(predicted - truth_value) > max(REL_TOL * abs(truth_value), absolute)


def score(document: dict, truth_rows: list[dict]) -> dict:
    mapped = {row["case_id"]: row for row in document["charges"]}
    recovery_miss = []
    organic_miss = []
    max_recovery = 0.0
    max_organic = 0.0
    for truth_row in truth_rows:
        got = mapped[truth_row["case_id"]]
        recovery_gap = abs(got["recovery"] - truth_row["recovery"]) / max(abs(truth_row["recovery"]), 1e-12)
        organic_gap = abs(got["organic_mol_m3"] - truth_row["organic_mol_m3"]) / max(abs(truth_row["organic_mol_m3"]), 1e-12)
        max_recovery = max(max_recovery, recovery_gap)
        max_organic = max(max_organic, organic_gap)
        if outside(got["recovery"], truth_row["recovery"], ABS_RECOVERY):
            recovery_miss.append(truth_row["case_id"])
        if outside(got["organic_mol_m3"], truth_row["organic_mol_m3"], ABS_ORGANIC):
            organic_miss.append(truth_row["case_id"])
    return {
        "max_relative_recovery": max_recovery,
        "max_relative_organic": max_organic,
        "recovery_misses": len(recovery_miss),
        "organic_misses": len(organic_miss),
    }


def baseline_sheet(forecast, runs, charges) -> dict:
    """Fit the controller temperature and the spreadsheet pH, and keep bad rows."""
    design = []
    observed = []
    for row in runs:
        design.append([
            float(row["ph_sheet"]),
            float(row["tic_c"]),
            float(row["solvent_ratio"]),
            float(row["stirrer_rps"]),
            float(row["contact_s"]),
            int(row["stage_count"]),
        ])
        observed.append(float(row["recovery_sheet"]))
    matrix = np.asarray(design, dtype=float)
    target = np.asarray(observed, dtype=float)

    def residual(theta):
        predicted = np.empty(len(target))
        for index, row in enumerate(matrix):
            predicted[index] = forecast.recovery_of(
                row[0], row[1], row[2], row[3], row[4], int(row[5]),
                theta[0], math.exp(theta[1]), theta[2], math.exp(theta[3]),
            )
        return predicted - target

    start = np.array([4.4, math.log(40.0), 1800.0, math.log(0.04)])
    bounds = (
        [3.2, math.log(8.0), 400.0, math.log(0.008)],
        [6.2, math.log(300.0), 7000.0, math.log(0.4)],
    )
    from scipy.optimize import least_squares
    fitted = least_squares(residual, start, bounds=bounds)
    rows = []
    for row in charges:
        recovery = forecast.recovery_of(
            float(row["ph_sheet"]), float(row["tic_c"]), float(row["solvent_ratio"]),
            float(row["stirrer_rps"]), float(row["contact_s"]), int(row["stage_count"]),
            float(fitted.x[0]), math.exp(float(fitted.x[1])), float(fitted.x[2]), math.exp(float(fitted.x[3])),
        )
        rows.append({
            "case_id": row["case_id"],
            "recovery": recovery,
            "organic_mol_m3": recovery * float(row["c_feed_mol_m3"]) / float(row["solvent_ratio"]),
        })
    return {"charges": rows}


def baseline_mean(runs, charges) -> dict:
    recoveries = []
    for row in runs:
        if row["status"] != "ok":
            continue
        recoveries.append(1.0 - float(row["c_aq_mol_m3"]) / float(row["c_feed_mol_m3"]))
    value = float(np.mean(recoveries))
    return {"charges": [
        {
            "case_id": row["case_id"],
            "recovery": value,
            "organic_mol_m3": value * float(row["c_feed_mol_m3"]) / float(row["solvent_ratio"]),
        }
        for row in charges
    ]}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    rng = np.random.default_rng(SEED)
    runs = build_runs(rng)
    charges, truth_rows = build_charges(rng, runs)
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "calibration.json").write_text(json.dumps(CALIBRATION, indent=2) + "\n")
    write_csv(DATA / "runs.csv", runs)
    write_csv(DATA / "charges.csv", charges)
    TRUTH.parent.mkdir(parents=True, exist_ok=True)
    TRUTH.write_text(json.dumps({"charges": truth_rows}, indent=2) + "\n")
    SEEDS.write_text(json.dumps({
        "seed": SEED,
        "pKa": PKA,
        "K_ref": K_REF,
        "temperature_slope_k": TEMPERATURE_SLOPE,
        "kLa_ref_per_s": KLA_REF,
        "aqueous_relative_noise": 0.004,
    }, indent=2) + "\n")

    forecast = load_forecast()
    calibration = json.loads((DATA / "calibration.json").read_text())
    run_rows = forecast.read_csv(DATA / "runs.csv")
    charge_rows = forecast.read_csv(DATA / "charges.csv")
    reference = forecast.build_predictions(DATA)
    other_start = np.array([5.1, math.log(90.0), 3200.0, math.log(0.09)])
    other_theta = forecast.fit_constants(run_rows, calibration, start=other_start)
    second = {"charges": forecast.forecasts(charge_rows, calibration, other_theta)}
    reference_score = score(reference, truth_rows)
    second_score = score(second, truth_rows)
    sheet_score = score(baseline_sheet(forecast, runs, charges), truth_rows)
    mean_score = score(baseline_mean(runs, charges), truth_rows)
    gap = 0.0
    for left, right in zip(reference["charges"], second["charges"]):
        gap = max(gap, abs(left["recovery"] - right["recovery"]))
    results = {
        "reference": reference_score,
        "alternate_start": second_score,
        "spreadsheet_columns": sheet_score,
        "mean_recovery": mean_score,
        "route_gap_max_abs_recovery": gap,
        "tolerances": {
            "relative": REL_TOL,
            "absolute_recovery": ABS_RECOVERY,
            "absolute_organic_mol_m3": ABS_ORGANIC,
        },
    }
    EVIDENCE.parent.mkdir(parents=True, exist_ok=True)
    EVIDENCE.write_text(json.dumps(results, indent=2) + "\n")
    if reference_score["recovery_misses"] or reference_score["organic_misses"]:
        raise SystemExit(f"reference misses charges: {results}")
    if second_score["recovery_misses"] or second_score["organic_misses"]:
        raise SystemExit(f"alternate start misses charges: {results}")
    if sheet_score["recovery_misses"] < len(truth_rows) // 2:
        raise SystemExit(f"spreadsheet baseline is too close: {sheet_score}")
    if mean_score["recovery_misses"] < len(truth_rows) // 2:
        raise SystemExit(f"mean baseline is too close: {mean_score}")
    ok_runs = [row for row in runs if row["status"] == "ok"]

    def span(rows, key):
        values = [float(row[key]) for row in rows]
        return min(values), max(values)

    for key in ("tic_c", "solvent_ratio", "stirrer_rps", "contact_s"):
        train_lo, train_hi = span(ok_runs, key)
        outside_key = [row for row in charges if float(row[key]) < train_lo or float(row[key]) > train_hi]
        if not outside_key:
            raise SystemExit(f"charges stay inside the historical {key} window")
    lines = []
    for path in sorted(DATA.rglob("*")):
        if path.is_file():
            lines.append(f"{sha256(path)}  environment/data/{path.relative_to(DATA)}")
    lines.append(f"{sha256(TRUTH)}  tests/truth/expected.json")
    CHECKSUMS.write_text("\n".join(lines) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
