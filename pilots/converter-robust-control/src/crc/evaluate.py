"""Evaluate one admissible or comparator controller on the public plants."""

from __future__ import annotations

import numpy as np

from .closed_loop import closed_loop_formula, formula_matches_interconnection, stack_controller
from .norms import (
    bounded_real_certificate,
    classify_stability,
    evaluate_transfer,
    refined_frequency_lower_bound,
    stability_report,
)
from .schema import d_sigma_max, gain_ok
from .params import (
    GAMMA_PASS,
    GAMMA_WITNESS,
    METHOD_AGREE_ABS,
    N_X_PHYS,
    STAB_PASS,
    STAB_WITNESS,
)

# Weight poles sit near -0.126 rad/s by construction. They are analysis
# states, not equipment modes, and they do not feed back into the loop.
# H-infinity only requires the augmented realization to be clearly Hurwitz.
HINF_HURWITZ = -1.0e-3


def _agree(a: float, b: float) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= max(METHOD_AGREE_ABS, 0.01 * scale)


def physical_closed_loop(plant: dict, K: dict) -> np.ndarray:
    """Internal dynamics of the scaled physical plant and the controller."""
    phys = {
        "A": plant["physical"]["A_scaled_physical"],
        "Bw": np.zeros((N_X_PHYS, 1)),
        "Bu": plant["physical"]["Bu_pu"],
        "Cy": plant["physical"]["Cy_scaled"],
        "Dyw": np.zeros((4, 1)),
        "Dyu": np.zeros((4, 4)),
        "Cz": plant["physical"]["Cy_scaled"],
        "Dzw": np.zeros((4, 1)),
        "Dzu": np.zeros((4, 4)),
    }
    return closed_loop_formula(phys, K)["A"]


def evaluate_blocks(
    blocks: list[dict],
    plants: list[dict],
    check_formula: bool = False,
    certify: bool = False,
) -> dict:
    sig = d_sigma_max(blocks)
    K = stack_controller(blocks)
    per_case = []
    for plant in plants:
        if check_formula:
            formula_matches_interconnection(plant, K)
        A_phys = physical_closed_loop(plant, K)
        phys_report = stability_report(A_phys)
        phys_class = classify_stability(phys_report, STAB_PASS)
        cl = closed_loop_formula(plant, K)
        report = evaluate_transfer(cl["A"], cl["B"], cl["C"], cl["D"], HINF_HURWITZ)
        stab = report["stability_class"]
        ham = report.get("hamiltonian") or {}
        slicot = (report.get("slicot_linf") or {}).get("value")
        refined = None
        agrees = False
        upper = None
        brl = None
        if stab == "pass" and slicot is not None:
            refined_report = refined_frequency_lower_bound(cl["A"], cl["B"], cl["C"], cl["D"])
            refined = float(refined_report["lower"])
            agrees = _agree(refined, float(slicot))
            upper = float(max(refined, slicot))
        hinf_pass = bool(
            phys_class == "pass"
            and stab == "pass"
            and agrees
            and upper is not None
            and upper <= GAMMA_PASS
        )
        if certify and stab == "pass":
            brl = bounded_real_certificate(cl["A"], cl["B"], cl["C"], cl["D"], GAMMA_PASS)
        witness = False
        if hinf_pass and certify:
            witness = bool(
                upper <= GAMMA_WITNESS
                and phys_report["abscissa"] + phys_report["uncertainty"] <= STAB_WITNESS
                and brl
                and brl.get("certified")
                and float(slicot) <= GAMMA_WITNESS
                and float(refined) <= GAMMA_WITNESS
            )
        per_case.append(
            {
                "case_id": plant["case_id"],
                "stability_class": phys_class,
                "augmented_hurwitz": stab,
                "abscissa": phys_report["abscissa"],
                "uncertainty": phys_report["uncertainty"],
                "eig_cond_max": phys_report.get("eig_cond_max"),
                "hamiltonian_status": ham.get("status"),
                "hamiltonian_lower": ham.get("lower"),
                "hamiltonian_upper": ham.get("upper"),
                "refined_lower": refined,
                "slicot_linf": slicot,
                "slicot_fpeak": (report.get("slicot_linf") or {}).get("fpeak"),
                "method_agree": agrees,
                "hinf_upper_used": upper,
                "bounded_real_at_1": None if brl is None else brl.get("certified"),
                "bounded_real_lmi_max_eig": None if brl is None else brl.get("lmi_max_eig"),
                "bounded_real_status": None if brl is None else brl.get("status"),
                "pass_stability": phys_class == "pass" and stab == "pass",
                "pass_hinf": hinf_pass,
                "witness_margin": witness,
                "hinf_status": (report.get("hinf") or {}).get("status"),
            }
        )
    all_pass = all(row["pass_stability"] and row["pass_hinf"] for row in per_case) and gain_ok(blocks)
    witness = all(row["witness_margin"] for row in per_case) and gain_ok(blocks)
    return {
        "d_sigma_max": sig,
        "gain_bound_pass": gain_ok(blocks),
        "pass_all": all_pass,
        "witness_margin": witness,
        "per_case": per_case,
        "blocks": blocks,
    }
