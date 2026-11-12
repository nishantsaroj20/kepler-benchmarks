#!/usr/bin/env python3
"""Seed the extraction campaign and the latent forecast truth.

The truth is the forward model in this file, evaluated at the seeded
constants on the charge inputs. The reference solution is scored by
running solution/forecast.py. It is not used to build the truth.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
from datetime import date, timedelta
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "environment" / "data"
TRUTH = ROOT / "tests" / "truth" / "expected.json"
EVIDENCE = ROOT / "authoring" / "evidence" / "calibration_results.json"
SEEDS = ROOT / "authoring" / "provenance" / "seeds.json"
CHECKSUMS = ROOT / "authoring" / "provenance" / "checksums.sha256"
FORECAST = ROOT / "solution" / "forecast.py"

SEED = 20260318
PKA = 4.75
PARTITION = 12.0
TEMPERATURE_SLOPE = 1800.0
DIMER = 0.02
KLA_REF = 0.045
T_REF = 298.15
AQUEOUS_RELATIVE_NOISE = 0.004
REL_TOL = 0.02
ABS_RECOVERY = 0.012
ABS_ORGANIC = 0.60

CALIBRATION = {
    "solvent": "4-methyl-2-pentanone",
    "solute": "K-17 carboxylic acid",
    "temperature_reference_k": T_REF,
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

EQ_BOUNDS = (
    [3.5, math.log(1.0), 200.0, math.log(1e-4)],
    [6.0, math.log(80.0), 5000.0, math.log(0.05)],
)
RATE_BOUNDS = ([math.log(0.005)], [math.log(0.3)])
JOINT_BOUNDS = (
    [3.5, math.log(1.0), 200.0, math.log(0.005)],
    [6.0, math.log(80.0), 5000.0, math.log(0.3)],
)


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


def forward(
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
    equilibrium: bool = False,
) -> tuple[float, float]:
    """Closed-contact balance with organic-phase dimerization.

    alpha is the monomer distribution of the neutral acid. The organic
    loading is alpha*Ca + 2*dimer*(alpha*Ca)^2, and each stage conserves
    solute between the aqueous holdup and its own solvent portion.
    """
    kelvin = temperature_c + 273.15
    neutral = 1.0 / (1.0 + 10.0 ** (ph - pka))
    alpha = partition * math.exp(slope * (1.0 / T_REF - 1.0 / kelvin)) * neutral
    beta = 2.0 * dimer * alpha * alpha
    portion = solvent_ratio / stages
    aqueous = feed
    for _ in range(stages):
        quad = portion * beta
        linear = 1.0 + portion * alpha
        if quad <= 1e-18:
            equilibrium_aqueous = aqueous / linear
        else:
            discriminant = linear * linear + 4.0 * quad * aqueous
            equilibrium_aqueous = (-linear + math.sqrt(discriminant)) / (2.0 * quad)
        if equilibrium:
            aqueous = equilibrium_aqueous
            continue
        reference = CALIBRATION["stirrer_reference_rps"]
        exponent = CALIBRATION["stirrer_exponent"]
        kla = kla_ref * (stirrer_rps / reference) ** exponent
        moved = 1.0 - math.exp(-kla * contact_s)
        aqueous = aqueous - moved * (aqueous - equilibrium_aqueous)
    removed = 1.0 - aqueous / feed
    organic = (feed - aqueous) / solvent_ratio
    return removed, organic


def truth_pair(feed, solvent_ratio, stages, ph, temperature_c, stirrer_rps, contact_s, equilibrium=False):
    return forward(
        feed, solvent_ratio, stages, ph, temperature_c, stirrer_rps, contact_s,
        PKA, PARTITION, TEMPERATURE_SLOPE, DIMER, KLA_REF, equilibrium,
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def dated(start: date, index: int) -> str:
    return (start + timedelta(days=int(index))).isoformat()


def row_record(case_id, day, window, ph, tic, solvent_ratio, stirrer, contact, stages, feed, status, aqueous, organic_assay):
    mv = millivolts(ph, window)
    return {
        "case_id": case_id,
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
        "c_org_assay_mol_m3": f"{organic_assay:.6f}",
        "recovery_sheet": f"{(organic_assay * solvent_ratio / feed):.6f}",
    }


def build_shakes(rng: np.random.Generator) -> list[dict]:
    rows = []
    for index in range(48):
        early = index < 22
        window = CALIBRATION["ph_windows"][0 if early else 1]
        day = dated(date(2024, 2, 8) if early else date(2024, 5, 21), index if early else index - 22)
        ph = float(rng.uniform(3.15, 6.25))
        tic = float(rng.uniform(17.0, 41.0))
        solvent_ratio = float(rng.uniform(0.25, 1.35))
        stirrer = 5.8
        contact = 900.0
        stages = 2 if index % 6 == 0 else 1
        feed = float(rng.uniform(4.0, 78.0))
        temperature = broth_celsius(tic)
        recovery, organic = truth_pair(feed, solvent_ratio, stages, ph, temperature, stirrer, contact, True)
        aqueous = max(feed * (1.0 - recovery) * (1.0 + AQUEOUS_RELATIVE_NOISE * float(rng.normal())), 1e-8)
        assay = organic + CALIBRATION["organic_assay_blank_mol_m3"] + 0.04 * float(rng.normal())
        rows.append(row_record(
            f"S{index + 1:03d}", day, window, ph, tic, solvent_ratio, stirrer, contact,
            stages, feed, "ok", aqueous, assay,
        ))
    return rows


def build_runs(rng: np.random.Generator) -> list[dict]:
    aborted = {7, 18, 29, 41, 58, 66, 77, 91}
    unmixed = {11, 24, 37, 52, 71, 88, 102, 114}
    rows = []
    for index in range(120):
        early = index < 56
        window = CALIBRATION["ph_windows"][0 if early else 1]
        origin = date(2024, 2, 6) if early else date(2024, 5, 20)
        day = dated(origin, index if early else index - 56)
        ph = float(rng.uniform(3.30, 5.95))
        tic = float(rng.uniform(18.0, 40.0))
        solvent_ratio = float(rng.uniform(0.28, 1.10))
        stirrer = float(rng.uniform(1.6, 5.2))
        contact = float(rng.uniform(18.0, 170.0))
        stages = 2 if rng.random() < 0.35 else 1
        feed = float(rng.uniform(10.0, 36.0))
        temperature = broth_celsius(tic)
        recovery, organic = truth_pair(feed, solvent_ratio, stages, ph, temperature, stirrer, contact, False)
        if index in aborted:
            status = "aborted"
            aqueous = feed * float(rng.uniform(0.15, 1.25))
            assay = float(rng.uniform(0.0, 20.0))
        elif index in unmixed:
            status = "unmixed"
            aqueous = feed * (1.0 - 0.04 * recovery)
            assay = 0.05 * organic
        else:
            status = "ok"
            aqueous = max(feed * (1.0 - recovery) * (1.0 + AQUEOUS_RELATIVE_NOISE * float(rng.normal())), 1e-8)
            assay = organic + CALIBRATION["organic_assay_blank_mol_m3"] + 0.03 * float(rng.normal())
        rows.append(row_record(
            f"R{index + 1:03d}", day, window, ph, tic, solvent_ratio, stirrer, contact,
            stages, feed, status, aqueous, assay,
        ))
    return rows


def build_charges(rng: np.random.Generator, runs: list[dict]) -> tuple[list[dict], list[dict]]:
    seen = {
        (row["date"][:7], row["tic_c"], row["solvent_ratio"], row["stirrer_rps"], row["contact_s"], row["stage_count"])
        for row in runs
    }
    charges = []
    truth = []
    index = 0
    while len(charges) < 52:
        early = len(charges) < 22
        window = CALIBRATION["ph_windows"][0 if early else 1]
        day = dated(date(2024, 3, 4) if early else date(2024, 6, 4), index)
        if early and day > "2024-05-18":
            day = "2024-05-18"
        if not early and day < "2024-05-19":
            day = "2024-05-19"
        # Most charges sit where recovery is still moving. A fully extracted
        # contact hides the isotherm, so the tail is only there to leave the
        # pilot window.
        if rng.random() < 0.82:
            ph = float(rng.uniform(4.55, 6.45))
            tic = float(rng.uniform(16.0, 46.0))
            solvent_ratio = float(rng.uniform(0.18, 0.85))
            stirrer = float(rng.uniform(1.3, 6.4))
            contact = float(rng.uniform(8.0, 95.0))
            feed = float(rng.uniform(18.0, 70.0))
        else:
            ph = float(rng.uniform(3.05, 6.60))
            tic = float(rng.uniform(11.0, 54.0))
            solvent_ratio = float(rng.uniform(0.12, 1.95))
            stirrer = float(rng.uniform(1.0, 8.0))
            contact = float(rng.uniform(5.0, 320.0))
            feed = float(rng.uniform(8.0, 70.0))
        stages = 2 if rng.random() < 0.40 else 1
        key = (day[:7], f"{tic:.4f}", f"{solvent_ratio:.6f}", f"{stirrer:.4f}", f"{contact:.2f}", str(stages))
        index += 1
        if key in seen or index > 5000:
            continue
        seen.add(key)
        temperature = broth_celsius(tic)
        recovery, organic = truth_pair(feed, solvent_ratio, stages, ph, temperature, stirrer, contact, False)
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
        truth.append({"case_id": case_id, "recovery": recovery, "organic_mol_m3": organic})
    return charges, truth


def _observed(row: dict) -> dict:
    window = window_for(row["date"])
    ph = window["intercept"] - (float(row["ph_mv"]) - window["millivolt_offset"]) / window["slope_mv_per_ph"]
    return {
        "ph": ph,
        "temperature": broth_celsius(float(row["tic_c"])),
        "solvent_ratio": float(row["solvent_ratio"]),
        "stirrer": float(row["stirrer_rps"]),
        "contact": float(row["contact_s"]),
        "stages": int(row["stage_count"]),
        "feed": float(row["c_feed_mol_m3"]),
        "aqueous": float(row["c_aq_mol_m3"]),
        "recovery": 1.0 - float(row["c_aq_mol_m3"]) / float(row["c_feed_mol_m3"]),
    }


def predict_rows(theta: dict, rows: list[dict], equilibrium: bool = False) -> list[dict]:
    out = []
    for row in rows:
        window = window_for(row["date"])
        ph = window["intercept"] - (float(row["ph_mv"]) - window["millivolt_offset"]) / window["slope_mv_per_ph"]
        temperature = broth_celsius(float(row["tic_c"]))
        feed = float(row["c_feed_mol_m3"])
        solvent_ratio = float(row["solvent_ratio"])
        recovery, organic = forward(
            feed, solvent_ratio, int(row["stage_count"]), ph, temperature,
            float(row["stirrer_rps"]), float(row["contact_s"]),
            theta["pka"], theta["partition"], theta["slope"], theta["dimer"], theta["kla"],
            equilibrium,
        )
        out.append({"case_id": row["case_id"], "recovery": recovery, "organic_mol_m3": organic})
    return out


def fit_equilibrium_then_rate(shakes: list[dict], runs: list[dict]) -> dict:
    """Identify the equilibrium on the shake study, then the rate on the pilot log."""
    shake_obs = [_observed(row) for row in shakes]
    run_obs = [_observed(row) for row in runs if row["status"] == "ok"]

    def equilibrium_residual(packed: np.ndarray) -> np.ndarray:
        pka, log_partition, slope, log_dimer = packed
        gap = np.empty(len(shake_obs))
        for index, row in enumerate(shake_obs):
            predicted, _organic = forward(
                row["feed"], row["solvent_ratio"], row["stages"], row["ph"], row["temperature"],
                row["stirrer"], row["contact"], pka, math.exp(log_partition), slope, math.exp(log_dimer),
                1.0, True,
            )
            gap[index] = predicted - row["recovery"]
        return gap

    equilibrium = least_squares(
        equilibrium_residual,
        np.array([5.05, math.log(20.0), 900.0, math.log(0.0015)]),
        bounds=EQ_BOUNDS,
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
    )
    if not equilibrium.success:
        raise RuntimeError(f"equilibrium route failed: {equilibrium.message}")
    pka = float(equilibrium.x[0])
    partition = math.exp(float(equilibrium.x[1]))
    slope = float(equilibrium.x[2])
    dimer = math.exp(float(equilibrium.x[3]))

    def rate_residual(packed: np.ndarray) -> np.ndarray:
        kla = math.exp(float(packed[0]))
        gap = np.empty(len(run_obs))
        for index, row in enumerate(run_obs):
            predicted, _organic = forward(
                row["feed"], row["solvent_ratio"], row["stages"], row["ph"], row["temperature"],
                row["stirrer"], row["contact"], pka, partition, slope, dimer, kla, False,
            )
            gap[index] = predicted - row["recovery"]
        return gap

    rate = least_squares(
        rate_residual,
        np.array([math.log(0.02)]),
        bounds=RATE_BOUNDS,
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
    )
    if not rate.success:
        raise RuntimeError(f"rate route failed: {rate.message}")
    return {
        "pka": pka,
        "partition": partition,
        "slope": slope,
        "dimer": dimer,
        "kla": math.exp(float(rate.x[0])),
        "equilibrium_nfev": int(equilibrium.nfev),
        "rate_nfev": int(rate.nfev),
    }


def fit_monomer_only(shakes: list[dict], runs: list[dict]) -> dict:
    """Reconcile the meters, then force the dimerization constant to zero."""
    observed = [_observed(row) for row in shakes]
    observed.extend(_observed(row) for row in runs if row["status"] == "ok")

    def residual(packed: np.ndarray) -> np.ndarray:
        pka, log_partition, slope, log_kla = packed
        gap = np.empty(len(observed))
        for index, row in enumerate(observed):
            predicted, _organic = forward(
                row["feed"], row["solvent_ratio"], row["stages"], row["ph"], row["temperature"],
                row["stirrer"], row["contact"], pka, math.exp(log_partition), slope, 0.0,
                math.exp(log_kla), False,
            )
            gap[index] = predicted - row["recovery"]
        return gap

    fitted = least_squares(
        residual,
        np.array([4.5, math.log(18.0), 1400.0, math.log(0.04)]),
        bounds=JOINT_BOUNDS,
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
    )
    return {
        "pka": float(fitted.x[0]),
        "partition": math.exp(float(fitted.x[1])),
        "slope": float(fitted.x[2]),
        "dimer": 0.0,
        "kla": math.exp(float(fitted.x[3])),
    }


def baseline_sheet(runs: list[dict], charges: list[dict]) -> dict:
    """Fit the controller temperature, the spreadsheet pH, and the blanked recovery."""
    design = []
    for row in runs:
        design.append({
            "ph": float(row["ph_sheet"]),
            "temperature": float(row["tic_c"]),
            "solvent_ratio": float(row["solvent_ratio"]),
            "stirrer": float(row["stirrer_rps"]),
            "contact": float(row["contact_s"]),
            "stages": int(row["stage_count"]),
            "feed": float(row["c_feed_mol_m3"]),
            "recovery": float(row["recovery_sheet"]),
        })

    def residual(packed: np.ndarray) -> np.ndarray:
        pka, log_partition, slope, log_kla = packed
        gap = np.empty(len(design))
        for index, row in enumerate(design):
            predicted, _organic = forward(
                row["feed"], row["solvent_ratio"], row["stages"], row["ph"], row["temperature"],
                row["stirrer"], row["contact"], pka, math.exp(log_partition), slope, 0.0,
                math.exp(log_kla), False,
            )
            gap[index] = predicted - row["recovery"]
        return gap

    fitted = least_squares(
        residual,
        np.array([4.5, math.log(18.0), 1400.0, math.log(0.04)]),
        bounds=JOINT_BOUNDS,
        xtol=1e-12,
        ftol=1e-12,
    )
    theta = {
        "pka": float(fitted.x[0]),
        "partition": math.exp(float(fitted.x[1])),
        "slope": float(fitted.x[2]),
        "dimer": 0.0,
        "kla": math.exp(float(fitted.x[3])),
    }
    rows = []
    for row in charges:
        feed = float(row["c_feed_mol_m3"])
        solvent_ratio = float(row["solvent_ratio"])
        recovery, organic = forward(
            feed, solvent_ratio, int(row["stage_count"]), float(row["ph_sheet"]), float(row["tic_c"]),
            float(row["stirrer_rps"]), float(row["contact_s"]),
            theta["pka"], theta["partition"], theta["slope"], 0.0, theta["kla"], False,
        )
        rows.append({"case_id": row["case_id"], "recovery": recovery, "organic_mol_m3": organic})
    return {"charges": rows}


def baseline_mean(runs: list[dict], charges: list[dict]) -> dict:
    recoveries = [
        1.0 - float(row["c_aq_mol_m3"]) / float(row["c_feed_mol_m3"])
        for row in runs
        if row["status"] == "ok"
    ]
    value = float(np.mean(recoveries))
    return {"charges": [
        {
            "case_id": row["case_id"],
            "recovery": value,
            "organic_mol_m3": value * float(row["c_feed_mol_m3"]) / float(row["solvent_ratio"]),
        }
        for row in charges
    ]}


def outside(predicted: float, truth_value: float, absolute: float) -> bool:
    return abs(predicted - truth_value) > max(REL_TOL * abs(truth_value), absolute)


def score(document: dict, truth_rows: list[dict]) -> dict:
    mapped = {row["case_id"]: row for row in document["charges"]}
    recovery_miss = 0
    organic_miss = 0
    charges_outside = 0
    max_rel_recovery = 0.0
    max_rel_organic = 0.0
    max_abs_recovery = 0.0
    max_abs_organic = 0.0
    for truth_row in truth_rows:
        got = mapped[truth_row["case_id"]]
        abs_recovery = abs(got["recovery"] - truth_row["recovery"])
        abs_organic = abs(got["organic_mol_m3"] - truth_row["organic_mol_m3"])
        max_rel_recovery = max(max_rel_recovery, abs_recovery / max(abs(truth_row["recovery"]), 1e-12))
        max_rel_organic = max(max_rel_organic, abs_organic / max(abs(truth_row["organic_mol_m3"]), 1e-12))
        max_abs_recovery = max(max_abs_recovery, abs_recovery)
        max_abs_organic = max(max_abs_organic, abs_organic)
        recovery_bad = outside(got["recovery"], truth_row["recovery"], ABS_RECOVERY)
        organic_bad = outside(got["organic_mol_m3"], truth_row["organic_mol_m3"], ABS_ORGANIC)
        recovery_miss += int(recovery_bad)
        organic_miss += int(organic_bad)
        charges_outside += int(recovery_bad or organic_bad)
    return {
        "max_relative_recovery": max_rel_recovery,
        "max_relative_organic": max_rel_organic,
        "max_absolute_recovery": max_abs_recovery,
        "max_absolute_organic_mol_m3": max_abs_organic,
        "recovery_misses": recovery_miss,
        "organic_misses": organic_miss,
        "charges_outside": charges_outside,
    }


def route_gap(left: dict, right: dict) -> dict:
    abs_recovery = 0.0
    abs_organic = 0.0
    rel_recovery = 0.0
    for a, b in zip(left["charges"], right["charges"]):
        abs_recovery = max(abs_recovery, abs(a["recovery"] - b["recovery"]))
        abs_organic = max(abs_organic, abs(a["organic_mol_m3"] - b["organic_mol_m3"]))
        rel_recovery = max(rel_recovery, abs(a["recovery"] - b["recovery"]) / max(abs(a["recovery"]), 1e-12))
    return {
        "max_absolute_recovery": abs_recovery,
        "max_absolute_organic_mol_m3": abs_organic,
        "max_relative_recovery": rel_recovery,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def reference_predictions() -> dict:
    destination = DATA / "_reference_predictions.json"
    subprocess.run(
        [sys.executable, str(FORECAST), "--data", str(DATA), "--output", str(destination)],
        check=True,
    )
    document = json.loads(destination.read_text())
    destination.unlink()
    return document


def main() -> None:
    rng = np.random.default_rng(SEED)
    shakes = build_shakes(rng)
    runs = build_runs(rng)
    charges, truth_rows = build_charges(rng, runs)
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "calibration.json").write_text(json.dumps(CALIBRATION, indent=2) + "\n")
    write_csv(DATA / "shakes.csv", shakes)
    write_csv(DATA / "runs.csv", runs)
    write_csv(DATA / "charges.csv", charges)
    TRUTH.parent.mkdir(parents=True, exist_ok=True)
    TRUTH.write_text(json.dumps({"charges": truth_rows}, indent=2) + "\n")
    SEEDS.write_text(json.dumps({
        "seed": SEED,
        "pKa": PKA,
        "partition_ref": PARTITION,
        "temperature_slope_k": TEMPERATURE_SLOPE,
        "dimerization_m3_per_mol": DIMER,
        "kLa_ref_per_s": KLA_REF,
        "aqueous_relative_noise": AQUEOUS_RELATIVE_NOISE,
        "equilibrium_shakes": len(shakes),
        "pilot_rows": len(runs),
        "charges": len(charges),
    }, indent=2) + "\n")

    reference = reference_predictions()
    sequential = fit_equilibrium_then_rate(shakes, runs)
    sequential_doc = {"charges": predict_rows(sequential, charges)}
    monomer = fit_monomer_only(shakes, runs)
    monomer_doc = {"charges": predict_rows(monomer, charges)}
    reference_score = score(reference, truth_rows)
    sequential_score = score(sequential_doc, truth_rows)
    monomer_score = score(monomer_doc, truth_rows)
    sheet_score = score(baseline_sheet(runs, charges), truth_rows)
    mean_score = score(baseline_mean(runs, charges), truth_rows)
    gap = route_gap(reference, sequential_doc)
    smallest_recovery = min(row["recovery"] for row in truth_rows)
    assay_recovery_sigma = AQUEOUS_RELATIVE_NOISE * (1.0 - smallest_recovery)
    results = {
        "reference_joint_least_squares": reference_score,
        "equilibrium_then_rate": {
            **sequential_score,
            "fitted": {key: sequential[key] for key in ("pka", "partition", "slope", "dimer", "kla")},
        },
        "route_gap": gap,
        "monomer_only_reconciled": monomer_score,
        "spreadsheet_columns": sheet_score,
        "mean_recovery": mean_score,
        "assay": {
            "aqueous_relative_noise": AQUEOUS_RELATIVE_NOISE,
            "smallest_true_recovery": smallest_recovery,
            "recovery_sigma_at_smallest_recovery": assay_recovery_sigma,
            "three_sigma_recovery": 3.0 * assay_recovery_sigma,
        },
        "tolerances": {
            "relative": REL_TOL,
            "absolute_recovery": ABS_RECOVERY,
            "absolute_organic_mol_m3": ABS_ORGANIC,
        },
    }
    EVIDENCE.parent.mkdir(parents=True, exist_ok=True)
    EVIDENCE.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
    if reference_score["recovery_misses"] or reference_score["organic_misses"]:
        raise SystemExit("reference misses charges")
    if sequential_score["recovery_misses"] or sequential_score["organic_misses"]:
        raise SystemExit("equilibrium-then-rate route misses charges")
    half = len(truth_rows) // 2
    if monomer_score["charges_outside"] < half:
        raise SystemExit(f"monomer-only baseline is too close: {monomer_score}")
    if sheet_score["charges_outside"] < half:
        raise SystemExit(f"spreadsheet baseline is too close: {sheet_score}")
    if mean_score["charges_outside"] < half:
        raise SystemExit(f"mean baseline is too close: {mean_score}")
    worse_relative = max(
        reference_score["max_relative_recovery"],
        reference_score["max_relative_organic"],
        sequential_score["max_relative_recovery"],
        sequential_score["max_relative_organic"],
        gap["max_relative_recovery"],
    )
    if REL_TOL <= 3.0 * worse_relative:
        raise SystemExit("relative tolerance does not clear three times the worse correct route")
    if ABS_RECOVERY <= max(gap["max_absolute_recovery"], 3.0 * assay_recovery_sigma):
        raise SystemExit("absolute recovery floor does not clear the route gap and the assay floor")
    if ABS_ORGANIC <= 3.0 * gap["max_absolute_organic_mol_m3"]:
        raise SystemExit("absolute organic floor does not clear three times the route gap")

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


if __name__ == "__main__":
    main()
