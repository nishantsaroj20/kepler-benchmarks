"""Identify four breathing delamination zones and forecast the closed-branch snap."""
import os
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, linear_sum_assignment
from nonlinear import Nonlinear
from split_plate import SplitPlate

KEYS = ("xc", "yc", "a", "b", "theta_deg", "d")


def membrane_parameter(fix):
    z = np.zeros(16)
    z[-2:] = fix
    z[0], z[1], z[7], z[8] = 0.2, 0.2, 0.3, 0.2
    z[2], z[3], z[9], z[10] = 0.08, 0.03, 0.08, 0.03
    return z


class Campaign:
    def __init__(self, spec, sensors):
        self.spec = spec
        self.plate = SplitPlate(spec, sensors)
        self.model = Nonlinear(spec, sensors)

    def apply(self, p, branch):
        p = np.asarray(p, float)
        open_branch = branch > 0
        loss = float(self.spec["damage_model"]["open_membrane_loss"]) if open_branch else 0.0
        self.model.damage_membrane(self.plate.separated(p), loss)
        Kb = self.plate.bending(p, open_branch)
        m = self.model
        m.K = (Kb * m.mass * m.h**2 - m.N * m.G) / m.E0

    def observe(self, force, p, branch):
        self.apply(p, branch)
        return self.model.observe(force, branch)

    def stable(self, force, p, branch):
        self.apply(p, branch)
        try:
            q, H = self.model.state(force, branch)
        except ValueError:
            return False
        return float(np.linalg.eigvalsh(H)[0]) > 0 and np.sign(q[0]) == branch

    def snap_load(self, p, branch=-1, lo=30.0, hi=120.0):
        if not self.stable(lo, p, branch):
            raise RuntimeError("closed branch is already unstable at the lower force")
        for _ in range(24):
            mid = 0.5 * (lo + hi)
            if self.stable(mid, p, branch):
                lo = mid
            else:
                hi = mid
        return float(lo)

    def _block(self, force, freq, shapes, p, branch):
        self.apply(p, branch)
        _q, f, Y = self.model.observe(force, branch)
        f = f[:7]
        Y = Y[:, :7]
        C = 10.0 * (f[:, None] / freq - 1.0) ** 2 + (1.0 - np.abs(Y.T @ shapes))
        ii, jj = linear_sum_assignment(C)
        ii = ii[np.argsort(jj)]
        sign = np.sign(np.sum(Y[:, ii] * shapes, axis=0) + 1e-15)
        block = np.r_[f[ii] / freq - 1.0, (Y[:, ii] * sign - shapes).ravel()]
        return block, float(block @ block)

    def match(self, p, obs):
        """Return the residual and the better preparation sign of each record."""
        blocks = []
        signs = []
        width = len(obs) * (len(obs[0][1]) + obs[0][2].size)
        for force, freq, shapes in obs:
            best = None
            for branch in (-1, 1):
                try:
                    block, cost = self._block(force, freq, shapes, p, branch)
                except ValueError:
                    continue
                if best is None or cost < best[0]:
                    best = (cost, block, branch)
            if best is None:
                return np.ones(width), [-1] * len(obs)
            blocks.append(best[1])
            signs.append(best[2])
        return np.concatenate(blocks), signs

    def residual(self, p, obs):
        return self.match(p, obs)[0]

    def _cost_of(self, p, obs):
        residual, signs = self.match(p, obs)
        return 0.5 * float(residual @ residual), signs

    def _place_centres(self, found, obs, zone, xs, ys):
        """Best centre for one zone. Parked zones (severity at its lower bound) do not repel."""
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
                cost = self._cost_of(trial, obs)[0]
                if best is None or cost < best[0]:
                    best = (cost, float(xc), float(yc))
        return best

    def _polish_free(self, found, obs, free_idx, max_nfev):
        """Bounded least squares on a subset of the parameter vector."""
        low, high = bounds(self.spec)
        free_idx = np.asarray(free_idx, int)
        scale = (high - low)[free_idx]
        start = np.clip((found[free_idx] - low[free_idx]) / scale, 0.0, 1.0)

        def residual(t):
            trial = found.copy()
            trial[free_idx] = low[free_idx] + t * scale
            return self.residual(trial, obs)

        sol = least_squares(
            residual,
            start,
            bounds=(np.zeros(len(free_idx)), np.ones(len(free_idx))),
            max_nfev=max_nfev,
            diff_step=0.02,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-11,
        )
        polished = found.copy()
        polished[free_idx] = low[free_idx] + sol.x * scale
        return polished, float(sol.cost)

    def restart_ellipses(self, found, obs):
        """Try other major-axis, severity and angle combinations on each patch.

        A continuous polish stays on the angle it was given. Area, severity and
        angle trade, so the lower residual can sit on a different angle.
        """
        best_cost, _signs = self._cost_of(found, obs)
        best = np.array(found, float)
        n_zones = (len(best) - 2) // 6
        for zone in range(n_zones):
            local_cost, local = best_cost, best
            for theta in np.arange(-90.0, 91.0, 15.0):
                for severity in (0.55, 0.70, 0.85, 0.95):
                    for major in (0.018, 0.024, 0.030, 0.036):
                        trial = best.copy()
                        minor = min(float(trial[6 * zone + 3]), major - 0.001)
                        if minor < 0.008:
                            continue
                        trial[6 * zone + 2] = major
                        trial[6 * zone + 3] = minor
                        trial[6 * zone + 4] = float(theta)
                        trial[6 * zone + 5] = severity
                        cost = self._cost_of(trial, obs)[0]
                        if cost < local_cost:
                            local_cost, local = cost, trial
            best_cost, best = local_cost, local
            print("ellipse", zone, round(best_cost, 8), flush=True)
        return best

    def fit(self, obs):
        """Place one ellipse at a time, updating the shared springs before the next search.

        Closed-branch bending still carries the split, so the springs are not a separate
        intact-plate fit. A centre scan at the wrong springs prefers a compromise blob;
        refitting the springs on the patches found so far moves the next scan onto a real patch.
        """
        low, high = bounds(self.spec)
        scale = high - low
        n_zones = int(self.spec["damage_model"]["zones_per_plate"])
        xs = np.arange(0.08, 0.44, 0.012)
        ys = np.arange(0.08, 0.33, 0.012)
        found = np.zeros(6 * n_zones + 2)
        for zone in range(n_zones):
            found[6 * zone] = 0.06
            found[6 * zone + 1] = 0.055
            found[6 * zone + 2] = 0.024
            found[6 * zone + 3] = 0.013
            found[6 * zone + 5] = 0.4
        found[-2], found[-1] = 1.0, 1.0
        for zone in range(n_zones):
            _cost, xc, yc = self._place_centres(found, obs, zone, xs, ys)
            found[6 * zone] = xc
            found[6 * zone + 1] = yc
            found[6 * zone + 5] = 0.85
            print("placed", zone, round(xc, 4), round(yc, 4), flush=True)
            free = [6 * z + k for z in range(zone + 1) for k in (2, 3, 4, 5)]
            free += [len(found) - 2, len(found) - 1]
            found, _cost = self._polish_free(found, obs, free, max_nfev=max(60, 6 * (len(free) + 1)))
        for zone in range(n_zones):
            wx = np.arange(found[6 * zone] - 0.018, found[6 * zone] + 0.019, 0.004)
            wy = np.arange(found[6 * zone + 1] - 0.018, found[6 * zone + 1] + 0.019, 0.004)
            _cost, xc, yc = self._place_centres(found, obs, zone, wx, wy)
            found[6 * zone] = xc
            found[6 * zone + 1] = yc
            found[6 * zone + 5] = max(found[6 * zone + 5], 0.7)
        centre_free = [6 * z + k for z in range(n_zones) for k in (0, 1)]
        centre_free += [len(found) - 2, len(found) - 1]
        found, _cost = self._polish_free(found, obs, centre_free, max_nfev=8 * (len(centre_free) + 1))
        seeds = [found]
        shifted = found.copy()
        for zone in range(n_zones):
            shifted[6 * zone + 4] = float(np.clip(shifted[6 * zone + 4] + 30.0, -90.0, 90.0))
        seeds.append(shifted)
        best = None
        for guess in seeds:
            v0 = np.clip((guess - low) / scale, 0.0, 1.0)
            sol = least_squares(
                lambda t: self.residual(low + t * scale, obs),
                v0,
                bounds=(np.zeros(len(low)), np.ones(len(low))),
                max_nfev=14 * (len(low) + 1),
                diff_step=0.004,
                ftol=1e-12,
                xtol=1e-12,
            )
            polished = low + sol.x * scale
            if best is None or sol.cost < best[0]:
                best = (float(sol.cost), polished, int(sol.nfev))
        _cost, p, nfev = best
        p = self.restart_ellipses(p, obs)
        v0 = np.clip((p - low) / scale, 0.0, 1.0)
        sol = least_squares(
            lambda t: self.residual(low + t * scale, obs),
            v0,
            bounds=(np.zeros(len(low)), np.ones(len(low))),
            max_nfev=10 * (len(low) + 1),
            diff_step=0.004,
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        p = low + sol.x * scale
        _residual, signs = self.match(p, obs)
        return float(sol.cost), p, signs, int(sol.nfev)

    def forecast(self, p):
        """Lowest four frequencies of the positive preparation at the requested forces."""
        rows = []
        sign = int(self.spec["forecast"]["preparation_sign"])
        count = int(self.spec["forecast"]["ascending_modes"])
        for force in self.spec["forecast"]["forces_N"]:
            freq = self.observe(float(force), p, sign)[1][:count]
            rows.append({"force_N": float(force), "frequencies_Hz": [float(v) for v in freq]})
        return rows


def bounds(spec):
    b = spec["damage_model"]["bounds"]
    n_zones = int(spec["damage_model"]["zones_per_plate"])
    low = [b[k][0] for k in KEYS] * n_zones + [spec["fixture_model"]["bounds"][0]] * 2
    high = [b[k][1] for k in KEYS] * n_zones + [spec["fixture_model"]["bounds"][1]] * 2
    return np.array(low, float), np.array(high, float)


def load_observations(data, spec):
    pts = np.loadtxt(data / "sensors.csv", delimiter=",", skiprows=1)[:, 1:]
    obs = []
    for rec in spec["records"]:
        folder = data / rec["record_id"]
        freq = np.loadtxt(folder / "frequencies.csv", delimiter=",", skiprows=1)[:, 1]
        shapes = np.loadtxt(folder / "mode_shapes.csv", delimiter=",", skiprows=1)[:, 1:]
        shapes = shapes / np.linalg.norm(shapes, axis=0)
        obs.append((float(rec["force_N"]), freq, shapes))
    return pts, obs


def artifact(spec, p, signs, snap, forecast):
    n_zones = int(spec["damage_model"]["zones_per_plate"])
    zones = []
    for j in range(n_zones):
        zones.append({k: float(p[6 * j + i]) for i, k in enumerate(KEYS)})
    return {
        "plates": [
            {
                "plate_id": spec["specimen_id"],
                "zones": zones,
                "fixture": {
                    "x0_multiplier": float(p[-2]),
                    "yL_multiplier": float(p[-1]),
                },
                "branches": {
                    rec["record_id"]: int(sign) for rec, sign in zip(spec["records"], signs)
                },
                "snap_load_N": float(snap),
                "open_branch_forecast_Hz": forecast,
            }
        ]
    }


def main():
    data = Path(os.environ.get("LDI_DATA_DIR", "/app/data"))
    out = Path(os.environ.get("LDI_OUT_FILE", "/app/output/damage.json"))
    spec = json.loads((data / "plate_specs.json").read_text())
    pts, obs = load_observations(data, spec)
    camp = Campaign(spec, pts)
    cost, p, signs, nfev = camp.fit(obs)
    if cost > 2.0e-5:
        raise RuntimeError(f"modal records unresolved (cost {cost})")
    snap = camp.snap_load(p, branch=-1)
    forecast = camp.forecast(p)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact(spec, p, signs, snap, forecast), indent=2) + "\n")
    print("completed", cost, "signs", signs, "snap", snap, "nfev", nfev, flush=True)


if __name__ == "__main__":
    main()
