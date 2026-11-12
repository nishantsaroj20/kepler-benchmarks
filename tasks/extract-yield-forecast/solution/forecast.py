#!/usr/bin/env python3
"""Forecast K-17 extraction charges from the plant log.

The shake tests are equilibrium contacts. Their organic loading rises
faster than the aqueous neutral concentration, which is the signature of
dimerization of the extracted monomer in the ketone. Each pilot contact
then approaches that equilibrium transfer on the plant stirrer law.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

T_REF = 298.15


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


def _neutral_alpha(ph: float, temperature_c: float, pka: float, partition: float, slope: float) -> float:
    kelvin = temperature_c + 273.15
    neutral = 1.0 / (1.0 + 10.0 ** (ph - pka))
    monomer = partition * math.exp(slope * (1.0 / T_REF - 1.0 / kelvin))
    return monomer * neutral


def _aqueous_equilibrium(ca_in: float, solvent_portion: float, alpha: float, beta: float) -> float:
    quad = solvent_portion * beta
    linear = 1.0 + solvent_portion * alpha
    if quad <= 1e-18:
        return ca_in / linear
    discriminant = linear * linear + 4.0 * quad * ca_in
    return (-linear + math.sqrt(discriminant)) / (2.0 * quad)


def _approach(stirrer_rps: float, contact_s: float, kla_ref: float, calibration: dict) -> float:
    reference = calibration["stirrer_reference_rps"]
    exponent = calibration["stirrer_exponent"]
    kla = kla_ref * (stirrer_rps / reference) ** exponent
    return 1.0 - math.exp(-kla * contact_s)


def recovery_of(
    feed: float,
    solvent_ratio: float,
    stages: int,
    ph: float,
    temperature_c: float,
    stirrer_rps: float,
    contact_s: float,
    pka: float,
    partition: float,
    slope: float,
    dimer: float,
    kla_ref: float,
    calibration: dict,
    equilibrium: bool = False,
) -> tuple[float, float]:
    alpha = _neutral_alpha(ph, temperature_c, pka, partition, slope)
    beta = 2.0 * dimer * alpha * alpha
    portion = solvent_ratio / stages
    aqueous = feed
    for _ in range(stages):
        equilibrium_aqueous = _aqueous_equilibrium(aqueous, portion, alpha, beta)
        if equilibrium:
            aqueous = equilibrium_aqueous
        else:
            moved = _approach(stirrer_rps, contact_s, kla_ref, calibration)
            aqueous = aqueous - moved * (aqueous - equilibrium_aqueous)
    removed = 1.0 - aqueous / feed
    organic = (feed - aqueous) / solvent_ratio
    return removed, organic


def _observed(row: dict[str, str], calibration: dict) -> tuple[float, float, float, float, float, int, float, float]:
    ph = broth_ph(row["date"], float(row["ph_mv"]), calibration)
    temperature = broth_celsius(float(row["tic_c"]), calibration)
    solvent_ratio = float(row["solvent_ratio"])
    stirrer = float(row["stirrer_rps"])
    contact = float(row["contact_s"])
    stages = int(row["stage_count"])
    feed = float(row["c_feed_mol_m3"])
    aqueous = float(row["c_aq_mol_m3"])
    return ph, temperature, solvent_ratio, stirrer, contact, stages, feed, 1.0 - aqueous / feed


def fit_constants(runs: list[dict[str, str]], shakes: list[dict[str, str]], calibration: dict) -> np.ndarray:
    rows = [row for row in runs if row["status"] == "ok"] + list(shakes)
    observed = [_observed(row, calibration) for row in rows]

    def residual(theta: np.ndarray) -> np.ndarray:
        pka, log_partition, slope, log_dimer, log_kla = theta
        partition = math.exp(log_partition)
        dimer = math.exp(log_dimer)
        kla_ref = math.exp(log_kla)
        gap = np.empty(len(observed))
        for index, row in enumerate(observed):
            predicted, _organic = recovery_of(
                row[6], row[2], row[5], row[0], row[1], row[3], row[4],
                pka, partition, slope, dimer, kla_ref, calibration,
            )
            gap[index] = predicted - row[7]
        return gap

    start = np.array([4.4, math.log(8.0), 1200.0, math.log(0.002), math.log(0.03)])
    bounds = (
        [3.5, math.log(1.0), 200.0, math.log(1e-4), math.log(0.005)],
        [6.0, math.log(80.0), 5000.0, math.log(0.05), math.log(0.3)],
    )
    fitted = least_squares(residual, start, bounds=bounds, xtol=1e-14, ftol=1e-14, gtol=1e-14)
    return fitted.x


def forecasts(charges: list[dict[str, str]], calibration: dict, theta: np.ndarray) -> list[dict]:
    pka = float(theta[0])
    partition = math.exp(float(theta[1]))
    slope = float(theta[2])
    dimer = math.exp(float(theta[3]))
    kla_ref = math.exp(float(theta[4]))
    rows = []
    for row in charges:
        ph = broth_ph(row["date"], float(row["ph_mv"]), calibration)
        temperature = broth_celsius(float(row["tic_c"]), calibration)
        solvent_ratio = float(row["solvent_ratio"])
        recovery, organic = recovery_of(
            float(row["c_feed_mol_m3"]), solvent_ratio, int(row["stage_count"]),
            ph, temperature, float(row["stirrer_rps"]), float(row["contact_s"]),
            pka, partition, slope, dimer, kla_ref, calibration,
        )
        rows.append({
            "case_id": row["case_id"],
            "recovery": recovery,
            "organic_mol_m3": organic,
        })
    return rows


def build_predictions(data_dir: Path) -> dict:
    calibration = json.loads((data_dir / "calibration.json").read_text())
    runs = read_csv(data_dir / "runs.csv")
    shakes = read_csv(data_dir / "shakes.csv")
    charges = read_csv(data_dir / "charges.csv")
    theta = fit_constants(runs, shakes, calibration)
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
