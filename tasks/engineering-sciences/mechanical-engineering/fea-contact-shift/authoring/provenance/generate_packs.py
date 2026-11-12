#!/usr/bin/env python3
"""Generate practice and held-out contact packs (synthetic, seeded)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def gap_profile(x: np.ndarray, patches: list[dict]) -> np.ndarray:
    """Multi-patch initial gap: pointwise minimum of parabolic asperities."""
    g = np.full_like(x, 1e3, dtype=float)
    for p in patches:
        g = np.minimum(g, p["clearance"] + (x - p["center"]) ** 2 / (2.0 * p["radius"]))
    return g


def load_history(kind: str, Fn: float, mu: float, n_steps: int, rng: np.random.Generator):
    fn = np.zeros(n_steps)
    ft = np.zeros(n_steps)
    if kind == "normal_then_shear":
        n1 = max(2, n_steps // 2)
        fn[:n1] = np.linspace(Fn / n1, Fn, n1)
        fn[n1:] = Fn
        ft[:n1] = 0.0
        ft_max = float(rng.uniform(0.25, 0.75)) * mu * Fn
        ft[n1:] = np.linspace(0.0, ft_max, n_steps - n1)
    elif kind == "proportional":
        ft_max = float(rng.uniform(0.2, 0.65)) * mu * Fn
        fn[:] = np.linspace(Fn / n_steps, Fn, n_steps)
        ft[:] = np.linspace(ft_max / n_steps, ft_max, n_steps)
    elif kind == "shear_reverse":
        n1 = max(2, n_steps // 3)
        n2 = max(n1 + 1, 2 * n_steps // 3)
        fn[:] = Fn
        fn[:n1] = np.linspace(Fn / n1, Fn, n1)
        ft_max = float(rng.uniform(0.3, 0.7)) * mu * Fn
        ft[:n1] = 0.0
        ft[n1:n2] = np.linspace(0.0, ft_max, n2 - n1)
        ft[n2:] = np.linspace(ft_max, -0.5 * ft_max, n_steps - n2)
    else:
        raise ValueError(kind)
    return fn, ft


def make_pack(spec: dict, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(int(spec["seed"]))
    x_min, x_max = spec["x_range"]
    n = int(spec["n_nodes"])
    x = np.linspace(x_min, x_max, n)
    g0 = gap_profile(x, spec["patches"])
    E = float(spec["E"])
    nu = float(spec["nu"])
    mu = float(spec["mu"])
    rref = float(spec.get("rref", 1.0))
    # Target Fn from first-patch Hertz estimate, then scale
    p0 = spec["patches"][0]
    Estar = E / (1.0 - nu * nu)
    a = float(spec.get("a_target", 0.15 * (x_max - x_min)))
    R = float(p0["radius"])
    Fn_h = np.pi * Estar * a**2 / (4.0 * R)
    Fn = float(spec.get("Fn_scale", 1.0)) * Fn_h
    # For multi-patch, increase load so both patches engage when intended
    if len(spec["patches"]) > 1:
        Fn *= float(spec.get("multi_load_factor", 1.8))
    fn, ft = load_history(spec["load_kind"], Fn, mu, int(spec["n_steps"]), rng)

    meta = {
        "pack_id": spec["pack_id"],
        "E": E,
        "nu": nu,
        "mu": mu,
        "rref": rref,
        "plane_strain": True,
        "kernel": "plane_strain_log",
        "n_nodes": n,
        "n_steps": int(spec["n_steps"]),
        "notes": "Synthetic multi-patch frictional contact pack",
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    with (out_dir / "interface.csv").open("w") as f:
        f.write("x,gap0\n")
        for xi, gi in zip(x, g0):
            f.write(f"{xi:.16e},{gi:.16e}\n")
    with (out_dir / "loads.csv").open("w") as f:
        f.write("step,Fn,Ft\n")
        for i, (a, b) in enumerate(zip(fn, ft), start=1):
            f.write(f"{i},{a:.16e},{b:.16e}\n")


def practice_specs() -> list[dict]:
    return [
        {
            "pack_id": "practice_01",
            "seed": 101,
            "x_range": (-0.006, 0.006),
            "n_nodes": 121,
            "n_steps": 8,
            "E": 210e9,
            "nu": 0.30,
            "mu": 0.25,
            "rref": 1.0,
            "a_target": 0.0011,
            "load_kind": "normal_then_shear",
            "patches": [{"center": 0.0, "radius": 0.05, "clearance": 0.0}],
        },
        {
            "pack_id": "practice_02",
            "seed": 102,
            "x_range": (-0.008, 0.008),
            "n_nodes": 141,
            "n_steps": 10,
            "E": 70e9,
            "nu": 0.33,
            "mu": 0.40,
            "rref": 1.0,
            "a_target": 0.0014,
            "load_kind": "proportional",
            "patches": [
                {"center": -0.0025, "radius": 0.04, "clearance": 0.0},
                {"center": 0.0028, "radius": 0.035, "clearance": 1.2e-6},
            ],
        },
        {
            "pack_id": "practice_03",
            "seed": 103,
            "x_range": (-0.007, 0.007),
            "n_nodes": 131,
            "n_steps": 12,
            "E": 110e9,
            "nu": 0.28,
            "mu": 0.18,
            "rref": 1.0,
            "a_target": 0.0012,
            "load_kind": "shear_reverse",
            "Fn_scale": 1.1,
            "patches": [{"center": 0.0004, "radius": 0.06, "clearance": 0.0}],
        },
        {
            "pack_id": "practice_04",
            "seed": 104,
            "x_range": (-0.009, 0.009),
            "n_nodes": 161,
            "n_steps": 10,
            "E": 200e9,
            "nu": 0.29,
            "mu": 0.35,
            "rref": 1.0,
            "a_target": 0.0013,
            "load_kind": "normal_then_shear",
            "multi_load_factor": 2.2,
            "patches": [
                {"center": -0.0035, "radius": 0.045, "clearance": 0.0},
                {"center": 0.0, "radius": 0.08, "clearance": 2.5e-6},
                {"center": 0.0032, "radius": 0.04, "clearance": 0.8e-6},
            ],
        },
    ]


def hidden_specs(seed: int = 9001, n_packs: int = 8) -> list[dict]:
    """Deliberately shifted families (geometry / mu / load path), not reseeds of practice."""
    rng = np.random.default_rng(seed)
    specs = []
    load_kinds = ["normal_then_shear", "proportional", "shear_reverse"]
    for i in range(n_packs):
        n_patches = int(rng.integers(1, 4))
        span = float(rng.uniform(0.0055, 0.010))
        patches = []
        for k in range(n_patches):
            patches.append(
                {
                    "center": float(rng.uniform(-0.45 * span, 0.45 * span)),
                    "radius": float(rng.uniform(0.03, 0.09)),
                    "clearance": float(rng.uniform(0.0, 3.5e-6) if k > 0 else 0.0),
                }
            )
        specs.append(
            {
                "pack_id": f"held_{i+1:02d}",
                "seed": int(rng.integers(1_000, 1_000_000)),
                "x_range": (-span, span),
                "n_nodes": int(rng.choice([111, 131, 151, 171])),
                "n_steps": int(rng.choice([8, 10, 12])),
                "E": float(rng.uniform(60e9, 230e9)),
                "nu": float(rng.uniform(0.25, 0.34)),
                "mu": float(rng.uniform(0.12, 0.45)),
                "rref": 1.0,
                "a_target": float(rng.uniform(0.0009, 0.0016)),
                "load_kind": str(load_kinds[i % len(load_kinds)]),
                "Fn_scale": float(rng.uniform(0.9, 1.25)),
                "multi_load_factor": float(rng.uniform(1.6, 2.4)),
                "patches": patches,
            }
        )
    return specs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--practice-out", type=Path, required=True)
    ap.add_argument("--hidden-out", type=Path, default=None)
    ap.add_argument("--hidden-seed", type=int, default=9001)
    args = ap.parse_args()
    for spec in practice_specs():
        make_pack(spec, args.practice_out / spec["pack_id"])
    if args.hidden_out is not None:
        for spec in hidden_specs(seed=args.hidden_seed):
            make_pack(spec, args.hidden_out / spec["pack_id"])
    print("done")


if __name__ == "__main__":
    main()
