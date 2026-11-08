"""Inverse that uses only the direct Newton plate.

No import of the reference condensation, its residual, or its least-squares fit.
Centres are scanned on the direct-model frequency and shape residual. Bounded
least squares then polishes axes, angles, severities and both spring multipliers
on that same direct residual. The JSON this writes is what the verifier scores.
"""
import json
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
    out = ROOT / "authoring/evidence/independent-inverse.json"
    out.write_text(json.dumps({"modal_cost": cost, "artifact": artifact}, indent=2) + "\n")
    print("completed", cost, "signs", signs, "snap", snap, flush=True)


if __name__ == "__main__":
    main()
