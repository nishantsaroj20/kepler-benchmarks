"""One fresh pilot entry point.

Creates a new output directory, writes metadata before the baseline search,
and does not read previous run directories.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crc.baselines import run_campaign
from crc.controllers import diagonal_pi
from crc.evaluate import evaluate_blocks
from crc.norms import bounded_real_certificate, hamiltonian_norm, slicot_linf
from crc.params import (
    A_NOISE_PU,
    A_REF_PU,
    A_V_PU,
    D_SIGMA_MAX,
    F_RES,
    GAMMA_PASS,
    GAMMA_WITNESS,
    PLANT_SET_ID,
    SCHEMA_VERSION,
    STAB_PASS,
    STAB_WITNESS,
    W_E_EPS,
    W_E_HZ,
    W_E_M,
    W_U_KINF,
    W_U_POLE,
    W_U_ZERO,
)
from crc.plant import build_generalized, compare_paths, plant_public_dict
from crc.schema import dump_controller, load_controller
def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def _git_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT.parents[1], text=True).strip()
    except Exception as exc:
        return f"unavailable:{exc}"


def norm_selftest() -> list[dict]:
    rows = []

    def add(name, A, B, C, D, expect):
        A = np.atleast_2d(np.asarray(A, float))
        B = np.atleast_2d(np.asarray(B, float))
        C = np.atleast_2d(np.asarray(C, float))
        D = np.atleast_2d(np.asarray(D, float))
        ham = hamiltonian_norm(A, B, C, D)
        sl = slicot_linf(A, B, C, D)
        brl = bounded_real_certificate(A, B, C, D, expect * 1.02 + 1e-4)
        rows.append(
            {
                "name": name,
                "expect": expect,
                "slicot": sl["value"],
                "slicot_fpeak": sl["fpeak"],
                "hamiltonian_status": ham["status"],
                "hamiltonian_lower": ham["lower"],
                "hamiltonian_upper": ham["upper"],
                "brl_certified_near_expect": brl["certified"],
                "brl_status": brl["status"],
                "slicot_abs_error": abs(sl["value"] - expect),
            }
        )

    add("lowpass", [[-1.0]], [[1.0]], [[1.0]], [[0.0]], 1.0)
    add("feedthrough", [[-1.0]], [[0.0]], [[0.0]], [[2.5]], 2.5)
    wn, zeta = 10.0, 0.05
    exact = 1.0 / (2.0 * zeta * np.sqrt(1.0 - zeta**2))
    add(
        "second_order",
        [[0.0, 1.0], [-wn**2, -2 * zeta * wn]],
        [[0.0], [wn**2]],
        [[1.0, 0.0]],
        [[0.0]],
        exact,
    )
    wn, zeta = 200.0, 1e-4
    exact = 1.0 / (2.0 * zeta * np.sqrt(1.0 - zeta**2))
    add(
        "narrow_resonance",
        [[0.0, 1.0], [-wn**2, -2 * zeta * wn]],
        [[0.0], [wn**2]],
        [[1.0, 0.0]],
        [[0.0]],
        exact,
    )
    add(
        "cancellation",
        [[-3.0, -2.0], [1.0, 0.0]],
        [[1.0], [0.0]],
        [[1.0, 1.0]],
        [[0.0]],
        0.5,
    )
    # Unstable: L-infinity is finite, H-infinity is not.
    A = np.array([[1.0]])
    B = np.array([[1.0]])
    C = np.array([[1.0]])
    D = np.array([[0.0]])
    sl = slicot_linf(A, B, C, D)
    rows.append(
        {
            "name": "unstable_pole",
            "expect_linf": 1.0,
            "slicot_linf": sl["value"],
            "hinf": "infinite",
            "note": "finite L-infinity must not be read as H-infinity",
        }
    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    out = args.output
    if out.exists():
        raise SystemExit(f"refusing to reuse existing output directory {out}")
    out.mkdir(parents=True)
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    t0 = time.perf_counter()
    frozen = {
        "plant_set_id": PLANT_SET_ID,
        "schema_version": SCHEMA_VERSION,
        "gamma_pass": GAMMA_PASS,
        "gamma_witness": GAMMA_WITNESS,
        "stab_pass_physical": STAB_PASS,
        "stab_witness_physical": STAB_WITNESS,
        "disturbance_pu": {"voltage": A_V_PU, "reference": A_REF_PU, "noise": A_NOISE_PU},
        "W_e": {"hz": W_E_HZ, "M": W_E_M, "eps": W_E_EPS},
        "W_u": {"k_inf": W_U_KINF, "zero_rad_s": W_U_ZERO, "pole_rad_s": W_U_POLE},
        "d_sigma_max": D_SIGMA_MAX,
        "stability_scope": "physical plant plus controller states; weight poles are excluded from the margin and included only in the Hurwitz gate for H-infinity",
        "f_res_hz": F_RES,
    }
    (out / "FROZEN_REQUIREMENTS.json").write_text(json.dumps(frozen, indent=2) + "\n")
    meta = {
        "run_id": out.name,
        "git_head": _git_hash(),
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "threads_requested": 1,
        "started_unix": time.time(),
    }
    (out / "run_metadata_started.json").write_text(json.dumps(meta, indent=2) + "\n")

    selftest = norm_selftest()
    (out / "norm_selftest.json").write_text(json.dumps(selftest, indent=2) + "\n")

    from crc.params import CASES

    plants = []
    publics = []
    comparisons = []
    for case in CASES:
        plant = build_generalized(case)
        plants.append(plant)
        pub = plant_public_dict(plant)
        publics.append(pub)
        comparisons.append({"case_id": case.case_id, **compare_paths(case, plant)})
        (out / "plants").mkdir(exist_ok=True)
        (out / "plants" / f"{case.case_id}.json").write_text(json.dumps(pub) + "\n")
    (out / "path_comparison.json").write_text(json.dumps(comparisons, indent=2) + "\n")

    campaign = run_campaign(plants)
    # Full independent replay of the best tied-PI grid controller and the zero law.
    pi_blocks = campaign["pi_grid_best_blocks"]
    full_pi = evaluate_blocks(pi_blocks, plants, check_formula=True, certify=True)
    full_zero = evaluate_blocks(campaign["zero_blocks"], plants, check_formula=True, certify=True)
    # Explicit hand design inside the grid: kp=0.4, ki=4.
    hand = diagonal_pi(0.4, 4.0, 0.4, 4.0)
    full_hand = evaluate_blocks(hand, plants, check_formula=False, certify=True)
    # Rounding / neighborhood.
    rng = np.random.default_rng(1401)
    neighbors = []
    for scale in (1e-6, 1e-4, 1e-3):
        for draw in range(4):
            noise = rng.normal(scale=scale, size=4)
            gains = np.array([0.4, 4.0, 0.4, 4.0]) * (1.0 + noise)
            ev = evaluate_blocks(diagonal_pi(*gains.tolist()), plants, certify=False)
            neighbors.append(
                {
                    "scale": scale,
                    "draw": draw,
                    "gains": gains.tolist(),
                    "pass_all": ev["pass_all"],
                    "worst_upper": max(r["hinf_upper_used"] or 1e9 for r in ev["per_case"]),
                    "worst_abscissa": max(r["abscissa"] for r in ev["per_case"]),
                }
            )
    rounded = [round(v, 6) for v in (0.4, 4.0, 0.4, 4.0)]
    serialized = out / "hand_pi_controller.json"
    dump_controller(serialized, hand)
    reloaded = load_controller(serialized)
    full_reloaded = evaluate_blocks(reloaded, plants, certify=False)

    def strip(ev):
        return {k: v for k, v in ev.items() if k != "blocks"}

    replay = {
        "grid_best": strip(full_pi),
        "hand_kp0p4_ki4": strip(full_hand),
        "zero": strip(full_zero),
        "reloaded_json": strip(full_reloaded),
        "neighbors": neighbors,
        "controller_sha256": _sha256_file(serialized),
    }
    (out / "independent_replay.json").write_text(json.dumps(replay, indent=2, default=float) + "\n")
    (out / "baseline_records.json").write_text(json.dumps(campaign["records"], indent=2, default=float) + "\n")
    (out / "centralized.json").write_text(json.dumps(campaign["centralized"], indent=2, default=float) + "\n")

    # Flat CSV of fast baseline rows.
    csv_path = out / "per_case_baselines.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "family",
                "start_id",
                "seed",
                "case_id",
                "pass_stability",
                "pass_hinf",
                "abscissa",
                "slicot_linf",
                "worst_slicot",
                "n_evaluations",
                "wall_s",
                "note",
            ],
        )
        writer.writeheader()
        for rec in campaign["records"]:
            fast = rec.get("fast") or {}
            for row in fast.get("per_case") or []:
                writer.writerow(
                    {
                        "family": rec["family"],
                        "start_id": rec["start_id"],
                        "seed": rec["seed"],
                        "case_id": row.get("case_id"),
                        "pass_stability": row.get("pass_stability"),
                        "pass_hinf": row.get("pass_hinf"),
                        "abscissa": row.get("abscissa"),
                        "slicot_linf": row.get("slicot_linf"),
                        "worst_slicot": fast.get("worst_slicot"),
                        "n_evaluations": rec["n_evaluations"],
                        "wall_s": rec["wall_s"],
                        "note": rec["note"],
                    }
                )

    simple_pass = bool(full_hand["pass_all"] and full_hand["witness_margin"])
    zero_pass = bool(full_zero["pass_all"])
    if simple_pass or zero_pass:
        verdict = "STOP"
        reason = (
            "A tied diagonal PI with gains kp=0.4 and ki=4 on both power channels "
            "meets internal stability and the stacked weighted H-infinity bound on every "
            "public case under SLICOT and the refined frequency lower bound, and the "
            "bounded-real lemma certifies the pass threshold. Ordinary tuning solves the "
            "frozen specification. Thresholds were not tightened."
        )
    else:
        verdict = "RECONSIDER"
        reason = (
            "The pre-registered simple PI did not clear both evaluation paths with the "
            "bounded-real certificate. No PROCEED claim is made."
        )
    verdict_doc = {
        "verdict": verdict,
        "reason": reason,
        "simple_pi_pass_all": full_hand["pass_all"],
        "simple_pi_witness": full_hand["witness_margin"],
        "zero_pass_all": zero_pass,
        "future_harbor_harness": "UNVERIFIED",
        "future_instruction_md": "not authored",
        "matlab_hinfstruct": "UNVERIFIED",
        "python_control_import": "UNVERIFIED_PYPI_QUARANTINE",
        "difficulty_band_current_screen": "UNVERIFIED",
        "wall_s": time.perf_counter() - t0,
        "campaign_wall_s": campaign["wall_s"],
    }
    (out / "stage_gate_verdicts.json").write_text(json.dumps(verdict_doc, indent=2) + "\n")
    versions = {
        "numpy": np.__version__,
        "python": platform.python_version(),
    }
    for name in ("scipy", "cvxpy", "clarabel", "slycot"):
        module = __import__(name)
        versions[name] = getattr(module, "__version__", "unknown")
    meta["versions"] = versions
    meta["finished_unix"] = time.time()
    meta["verdict"] = verdict
    (out / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(verdict_doc, indent=2))


if __name__ == "__main__":
    main()
