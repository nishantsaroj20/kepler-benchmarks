"""Inverse that uses only the direct Newton plate.

No import of the reference condensation, its residual, or its least-squares fit.
Centres are scanned on the direct-model frequency and shape residual. Bounded
least squares then polishes axes, angles, severities and both spring multipliers
on that same direct residual. The JSON this writes is what the verifier scores.
"""
import itertools
import json
import os
import subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares, linear_sum_assignment

from independent import DirectPlate

ROOT = Path(__file__).resolve().parents[2]
KEYS = ("xc", "yc", "a", "b", "theta_deg", "d")


def sensor_basis(spec, sensors):
    length, width = float(spec["Lx_m"]), float(spec["Ly_m"])
    modes = [(m, n) for m in range(1, 4) for n in range(1, 4)]
    basis = []
    for x, y in sensors:
        basis.append([np.sin(m * np.pi * x / length) * np.sin(n * np.pi * y / width) for m, n in modes])
    return np.asarray(basis, float)


def load():
    spec = json.loads((ROOT / "environment/data/plate_specs.json").read_text())
    data = ROOT / "environment/data"
    sensors = np.loadtxt(data / "sensors.csv", delimiter=",", skiprows=1)[:, 1:]
    obs = []
    for rec in spec["records"]:
        folder = data / rec["record_id"]
        freq = np.loadtxt(folder / "frequencies.csv", delimiter=",", skiprows=1)[:, 1]
        shapes = np.loadtxt(folder / "mode_shapes.csv", delimiter=",", skiprows=1)[:, 1:]
        shapes = shapes / np.linalg.norm(shapes, axis=0)
        obs.append((float(rec["force_N"]), freq, shapes))
    bounds = spec["damage_model"]["bounds"]
    n_zones = int(spec["damage_model"]["zones_per_plate"])
    low = np.array([bounds[k][0] for k in KEYS] * n_zones + [spec["fixture_model"]["bounds"][0]] * 2, float)
    high = np.array([bounds[k][1] for k in KEYS] * n_zones + [spec["fixture_model"]["bounds"][1]] * 2, float)
    return spec, obs, sensors, low, high


def assign_block(predicted_f, predicted_y, observed_f, observed_y):
    """Frequency-and-shape assignment. Predicted columns are the direct-model eigenvectors."""
    cost = 10.0 * (predicted_f[:, None] / observed_f - 1.0) ** 2 + (1.0 - np.abs(predicted_y.T @ observed_y))
    rows, cols = linear_sum_assignment(cost)
    order = np.argsort(cols)
    rows = rows[order]
    sign = np.sign(np.sum(predicted_y[:, rows] * observed_y, axis=0) + 1e-15)
    return np.r_[predicted_f[rows] / observed_f - 1.0, (predicted_y[:, rows] * sign - observed_y).ravel()]


class DirectSearch:
    def __init__(self, spec, obs, sensors):
        self.spec = spec
        self.obs = obs
        self.model = DirectPlate(spec)
        self.basis = sensor_basis(spec, sensors)

    def modes(self, force, p, branch):
        self.model.prepare(p, branch > 0)
        _q, hessian = self.model.equilibrium(force, branch)
        values, vectors = eigh(hessian)
        freq = np.sqrt(values) * self.model.omega / (2.0 * np.pi)
        shapes = self.basis @ vectors[:, :7]
        shapes = shapes / np.linalg.norm(shapes, axis=0)
        return freq, shapes

    def residual(self, p):
        blocks = []
        signs = []
        width = len(self.obs) * (len(self.obs[0][1]) + self.obs[0][2].size)
        n_zones = (len(p) - 2) // 6
        for force, freq, shapes in self.obs:
            best = None
            for branch in (-1, 1):
                try:
                    predicted_f, predicted_y = self.modes(force, p, branch)
                except (RuntimeError, ValueError, np.linalg.LinAlgError):
                    continue
                block = assign_block(predicted_f[:7], predicted_y, freq, shapes)
                err = float(block @ block)
                if best is None or err < best[0]:
                    best = (err, block, branch)
            if best is None:
                return np.ones(width + n_zones), [-1] * len(self.obs)
            blocks.append(best[1])
            signs.append(best[2])
        # One-sided so an axis swap (a <= b) is visible to the Jacobian and a valid ellipse is not.
        gap = [5.0 * max(0.0, p[6 * zone + 3] + 1e-4 - p[6 * zone + 2]) for zone in range(n_zones)]
        return np.concatenate(blocks + [np.asarray(gap, float)]), signs

    def cost(self, p):
        block, signs = self.residual(p)
        return 0.5 * float(block @ block), signs

    def place(self, found, zone, xs, ys):
        n_zones = (len(found) - 2) // 6
        best = None
        for yc in ys:
            for xc in xs:
                if any(
                    z != zone and found[6 * z + 5] > 0.5 and np.hypot(xc - found[6 * z], yc - found[6 * z + 1]) < 0.07
                    for z in range(n_zones)
                ):
                    continue
                trial = found.copy()
                trial[6 * zone] = xc
                trial[6 * zone + 1] = yc
                trial[6 * zone + 5] = 0.85
                value = self.cost(trial)[0]
                if best is None or value < best[0]:
                    best = (value, float(xc), float(yc))
        return best

    def polish(self, found, low, high, free, max_nfev, diff_step):
        free = np.asarray(free, int)
        scale = (high - low)[free]
        start = np.clip((found[free] - low[free]) / scale, 0.0, 1.0)

        def residual(t):
            trial = found.copy()
            trial[free] = low[free] + np.asarray(t, float) * scale
            return self.residual(trial)[0]

        sol = least_squares(
            residual,
            start,
            bounds=(np.zeros(len(free)), np.ones(len(free))),
            max_nfev=max_nfev,
            diff_step=diff_step,
            ftol=1e-14,
            xtol=1e-14,
        )
        found = found.copy()
        found[free] = low[free] + sol.x * scale
        return found, float(sol.cost)


def main():
    spec, obs, sensors, low, high = load()
    search = DirectSearch(spec, obs, sensors)
    n_zones = int(spec["damage_model"]["zones_per_plate"])
    found = np.zeros(6 * n_zones + 2)
    for zone in range(n_zones):
        found[6 * zone : 6 * zone + 6] = (0.06, 0.055, 0.024, 0.013, 0.0, 0.4)
    found[-2:] = (1.0, 1.0)
    xs = np.arange(0.08, 0.44, 0.012)
    ys = np.arange(0.08, 0.33, 0.012)
    for zone in range(n_zones):
        _cost, xc, yc = search.place(found, zone, xs, ys)
        found[6 * zone] = xc
        found[6 * zone + 1] = yc
        found[6 * zone + 5] = 0.85
        free = [6 * z + k for z in range(zone + 1) for k in (2, 3, 4, 5)]
        free += [len(found) - 2, len(found) - 1]
        found, cost = search.polish(found, low, high, free, max_nfev=18, diff_step=0.02)
        print(
            "centre",
            zone,
            round(float(found[6 * zone]), 4),
            round(float(found[6 * zone + 1]), 4),
            "springs",
            round(float(found[-2]), 4),
            round(float(found[-1]), 4),
            "cost",
            cost,
            flush=True,
        )
    for zone in range(n_zones):
        wx = np.arange(found[6 * zone] - 0.018, found[6 * zone] + 0.019, 0.004)
        wy = np.arange(found[6 * zone + 1] - 0.018, found[6 * zone + 1] + 0.019, 0.004)
        _cost, xc, yc = search.place(found, zone, wx, wy)
        found[6 * zone] = xc
        found[6 * zone + 1] = yc
        found[6 * zone + 5] = max(found[6 * zone + 5], 0.7)
        print("refine", zone, round(xc, 4), round(yc, 4), flush=True)
    seeds = [found]
    shifted = found.copy()
    for zone in range(n_zones):
        shifted[6 * zone + 4] = float(np.clip(shifted[6 * zone + 4] + 30.0, -90.0, 90.0))
    seeds.append(shifted)
    best = None
    for guess in seeds:
        polished, cost = search.polish(guess, low, high, list(range(len(guess))), 48, 0.01)
        if best is None or cost < best[0]:
            best = (cost, polished)
        print("polish", cost, flush=True)
    cost, found = best
    _score, signs = search.cost(found)
    snap = float(search.model.snap(found))
    forecast = search.model.forecast(found)
    zones = []
    for zone in range(n_zones):
        zones.append({key: float(found[6 * zone + i]) for i, key in enumerate(KEYS)})
    artifact = {
        "plates": [
            {
                "plate_id": spec["specimen_id"],
                "zones": zones,
                "fixture": {"x0_multiplier": float(found[-2]), "yL_multiplier": float(found[-1])},
                "branches": {rec["record_id"]: int(sign) for rec, sign in zip(spec["records"], signs)},
                "snap_load_N": snap,
                "open_branch_forecast_Hz": forecast,
            }
        ]
    }
    publish(cost, artifact)
    print("completed", cost, "signs", signs, "snap", snap, flush=True)


def angle_gap(left, right):
    return abs((left - right + 90.0) % 180.0 - 90.0)


def matched_errors(artifact):
    """Largest miss of each field after the same one-to-one zone pairing the verifier uses."""
    truth = json.loads((ROOT / "tests/ground_truth.json").read_text())
    plate = truth["plates"][0]
    tol = truth["tolerances"]
    got = artifact["plates"][0]
    best = None
    for order in itertools.permutations(range(len(plate["zones"]))):
        rows = []
        for i, j in enumerate(order):
            left, right = got["zones"][i], plate["zones"][j]
            rows.append(
                {
                    "dx": abs(left["xc"] - right["xc"]),
                    "dy": abs(left["yc"] - right["yc"]),
                    "a": abs(left["a"] - right["a"]),
                    "b": abs(left["b"] - right["b"]),
                    "theta": angle_gap(left["theta_deg"], right["theta_deg"]),
                    "d": abs(left["d"] - right["d"]),
                }
            )
        ok = all(
            row["dx"] <= tol["zone_centre_m"]
            and row["dy"] <= tol["zone_centre_m"]
            and row["a"] <= tol["zone_a_m"]
            and row["b"] <= tol["zone_b_m"]
            and row["theta"] <= tol["zone_theta_deg"]
            and row["d"] <= tol["zone_d"]
            for row in rows
        )
        score = sum(row["dx"] + row["dy"] for row in rows)
        if best is None or (ok and not best[0]) or (ok == best[0] and score < best[1]):
            best = (ok, score, rows)
    rows = best[2]
    forecast = [
        abs(a / b - 1.0)
        for got_row, true_row in zip(got["open_branch_forecast_Hz"], plate["open_branch_forecast_Hz"])
        for a, b in zip(got_row["frequencies_Hz"], true_row["frequencies_Hz"])
    ]
    return {
        "inside_caps": bool(best[0]),
        "centre_x_m": max(row["dx"] for row in rows),
        "centre_y_m": max(row["dy"] for row in rows),
        "a_m": max(row["a"] for row in rows),
        "b_m": max(row["b"] for row in rows),
        "theta_deg": max(row["theta"] for row in rows),
        "d": max(row["d"] for row in rows),
        "x0": abs(got["fixture"]["x0_multiplier"] - plate["fixture"]["x0_multiplier"]),
        "yL": abs(got["fixture"]["yL_multiplier"] - plate["fixture"]["yL_multiplier"]),
        "snap_relative": abs(got["snap_load_N"] / plate["snap_load_N"] - 1.0),
        "forecast_max_relative": max(forecast),
    }


def score_artifact(artifact, path):
    path.write_text(json.dumps(artifact, indent=2) + "\n")
    proc = subprocess.run(
        ["pytest", "tests/test_outputs.py", "-q", "--tb=line"],
        cwd=ROOT,
        env={**os.environ, "LDI_ARTIFACT": str(path)},
        capture_output=True,
        text=True,
    )
    tail = "\n".join((proc.stdout + proc.stderr).strip().splitlines()[-3:])
    return {"returncode": proc.returncode, "summary": tail}


def publish(cost, artifact):
    """Write the direct inverse and the real-verifier scores of that JSON."""
    passed = score_artifact(artifact, Path("/tmp/damage-independent.json"))
    midpoint = json.loads(json.dumps(artifact))
    for zone in midpoint["plates"][0]["zones"]:
        zone.update({"xc": 0.25, "yc": 0.20, "a": 0.027, "b": 0.014, "theta_deg": 0.0, "d": 0.67})
    rejected = score_artifact(midpoint, Path("/tmp/damage-midpoint.json"))
    out = ROOT / "authoring/evidence/independent-inverse.json"
    out.write_text(
        json.dumps(
            {
                "method": (
                    "Full inverse on the direct Newton plate only. Strains are integrated on the "
                    "30 by 26 Gauss rule and equilibria are continued by Newton. The residual is "
                    "that model's frequencies and shapes. The script does not import the reference "
                    "condensation, its residual, or its least-squares fit, and it does not start "
                    "from the authored parameters."
                ),
                "modal_cost": cost,
                "matched_errors": matched_errors(artifact),
                "verifier": passed,
                "rejected_interval_midpoint": rejected,
                "artifact": artifact,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
