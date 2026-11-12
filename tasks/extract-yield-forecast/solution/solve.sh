#!/bin/sh
set -eu
mkdir -p /app /app/output
cat > /tmp/forecast.py << 'EXTRACT_FORECAST_PY_EOF'
#!/usr/bin/env python3
"""Forecast weak-acid extraction recoveries from the plant log."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

T_REF = 298.15
N_REF = 3.2
SPEED_POWER = 0.85


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def broth_celsius(tic_c: float, calibration: dict) -> float:
    coeff = calibration["thermocouple_celsius"]
    return coeff["a0"] + coeff["a1"] * tic_c + coeff["a2"] * tic_c * tic_c


def window_for(day: str, calibration: dict) -> dict:
    chosen = calibration["ph_windows"][0]
    for window in calibration["ph_windows"]:
        if window["start"] <= day <= window["end"]:
            chosen = window
    return chosen


def broth_ph(day: str, millivolts: float, calibration: dict) -> float:
    window = window_for(day, calibration)
    return window["intercept"] - (millivolts - window["millivolt_offset"]) / window["slope_mv_per_ph"]


def recovery_of(
    ph: float,
    temperature_c: float,
    solvent_ratio: float,
    stirrer_rps: float,
    contact_s: float,
    stages: int,
    pka: float,
    k_ref: float,
    temperature_slope: float,
    kla_ref: float,
) -> float:
    kelvin = temperature_c + 273.15
    neutral = 1.0 / (1.0 + 10.0 ** (ph - pka))
    distribution = k_ref * neutral * math.exp(temperature_slope * (1.0 / T_REF - 1.0 / kelvin))
    approach = 1.0 - math.exp(-kla_ref * (stirrer_rps / N_REF) ** SPEED_POWER * contact_s)
    portion = solvent_ratio / stages
    remaining = 1.0
    for _ in range(stages):
        step = approach * (distribution * portion) / (1.0 + distribution * portion)
        remaining *= 1.0 - step
    return 1.0 - remaining


def usable(row: dict[str, str]) -> bool:
    return row["status"] == "ok"


def fit_constants(runs: list[dict[str, str]], calibration: dict, start: np.ndarray | None = None) -> np.ndarray:
    design = []
    observed = []
    for row in runs:
        if not usable(row):
            continue
        ph = broth_ph(row["date"], float(row["ph_mv"]), calibration)
        temperature = broth_celsius(float(row["tic_c"]), calibration)
        solvent_ratio = float(row["solvent_ratio"])
        stirrer = float(row["stirrer_rps"])
        contact = float(row["contact_s"])
        stages = int(row["stage_count"])
        feed = float(row["c_feed_mol_m3"])
        aqueous = float(row["c_aq_mol_m3"])
        design.append([ph, temperature, solvent_ratio, stirrer, contact, stages])
        observed.append(1.0 - aqueous / feed)
    matrix = np.asarray(design, dtype=float)
    target = np.asarray(observed, dtype=float)

    def residual(theta: np.ndarray) -> np.ndarray:
        pka, log_k, slope, log_kla = theta
        predicted = np.empty(len(target))
        for index, row in enumerate(matrix):
            predicted[index] = recovery_of(
                row[0], row[1], row[2], row[3], row[4], int(row[5]),
                pka, math.exp(log_k), slope, math.exp(log_kla),
            )
        return predicted - target

    if start is None:
        start = np.array([4.4, math.log(40.0), 1800.0, math.log(0.04)])
    bounds = (
        [3.2, math.log(8.0), 400.0, math.log(0.008)],
        [6.2, math.log(300.0), 7000.0, math.log(0.4)],
    )
    fitted = least_squares(residual, start, bounds=bounds, xtol=1e-14, ftol=1e-14, gtol=1e-14)
    return fitted.x


def forecasts(charges: list[dict[str, str]], calibration: dict, theta: np.ndarray) -> list[dict]:
    pka, log_k, slope, log_kla = (float(theta[0]), float(theta[1]), float(theta[2]), float(theta[3]))
    k_ref = math.exp(log_k)
    kla_ref = math.exp(log_kla)
    rows = []
    for row in charges:
        ph = broth_ph(row["date"], float(row["ph_mv"]), calibration)
        temperature = broth_celsius(float(row["tic_c"]), calibration)
        solvent_ratio = float(row["solvent_ratio"])
        recovery = recovery_of(
            ph, temperature, solvent_ratio, float(row["stirrer_rps"]), float(row["contact_s"]),
            int(row["stage_count"]), pka, k_ref, slope, kla_ref,
        )
        organic = recovery * float(row["c_feed_mol_m3"]) / solvent_ratio
        rows.append({
            "case_id": row["case_id"],
            "recovery": recovery,
            "organic_mol_m3": organic,
        })
    return rows


def build_predictions(data_dir: Path) -> dict:
    calibration = json.loads((data_dir / "calibration.json").read_text())
    runs = read_csv(data_dir / "runs.csv")
    charges = read_csv(data_dir / "charges.csv")
    theta = fit_constants(runs, calibration)
    return {"charges": forecasts(charges, calibration, theta)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    document = build_predictions(args.data)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n")


if __name__ == "__main__":
    main()
EXTRACT_FORECAST_PY_EOF
python3 /tmp/forecast.py --data /app/data --output /app/output/predictions.json
test -s /app/output/predictions.json
