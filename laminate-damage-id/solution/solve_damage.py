"""Identify two breathing delamination zones and forecast the closed-branch snap."""
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
        self.model.assemble(membrane_parameter(p[-2:]))
        Kb = self.plate.bending(np.asarray(p, float), branch > 0)
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

    def fit(self, obs):
        """Scan one elliptical support at a time, then polish the best pair.

        A trial zone does not change the plate until its ellipse overlaps the split,
        so random multistart stays in decoys whose residual is several times the noise.
        """
        low, high = bounds(self.spec)
        scale = high - low
        xs = np.linspace(0.09, 0.41, 9)
        ys = np.linspace(0.09, 0.31, 7)
        base = np.array([0.20, 0.16, 0.034, 0.015, 0.0, 0.78, 0.30, 0.22, 0.030, 0.014, 0.0, 0.0, 1.05, 0.90])
        ranked = []
        for yc in ys:
            for xc in xs:
                trial = base.copy()
                trial[0], trial[1] = xc, yc
                cost, _signs = self._cost_of(trial, obs)
                ranked.append((cost, xc, yc))
        ranked.sort()
        anchors = []
        for _cost, xc, yc in ranked:
            if all(np.hypot(xc - ax, yc - ay) > 0.06 for ax, ay in anchors):
                anchors.append((xc, yc))
            if len(anchors) == 2:
                break
        (x1, y1), (x2, y2) = anchors
        seeds = []
        for xc1, yc1, xc2, yc2 in ((x1, y1, x2, y2), (x2, y2, x1, y1)):
            for th1, th2 in ((15.0, -40.0), (40.0, -15.0), (-20.0, 30.0)):
                guess = base.copy()
                guess[11] = 0.72
                guess[0], guess[1], guess[6], guess[7] = xc1, yc1, xc2, yc2
                guess[4], guess[10] = th1, th2
                seeds.append(guess)
        best = None
        for guess in seeds:
            v0 = np.clip((guess - low) / scale, 0.0, 1.0)
            sol = least_squares(
                lambda t: self.residual(low + t * scale, obs),
                v0,
                bounds=(np.zeros(len(low)), np.ones(len(low))),
                max_nfev=16,
                diff_step=0.012,
                ftol=1e-14,
                xtol=1e-14,
            )
            polished = low + sol.x * scale
            if best is None or sol.cost < best[0]:
                best = (float(sol.cost), polished, int(sol.nfev))
        _cost, p, nfev = best
        _residual, signs = self.match(p, obs)
        return best[0], p, signs, nfev

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
    low = [b[k][0] for k in KEYS] * 2 + [spec["fixture_model"]["bounds"][0]] * 2
    high = [b[k][1] for k in KEYS] * 2 + [spec["fixture_model"]["bounds"][1]] * 2
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
    zones = []
    for j in range(2):
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
    if cost > 8e-6:
        raise RuntimeError(f"modal records unresolved (cost {cost})")
    snap = camp.snap_load(p, branch=-1)
    forecast = camp.forecast(p)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact(spec, p, signs, snap, forecast), indent=2) + "\n")
    print("completed", cost, "signs", signs, "snap", snap, "nfev", nfev, flush=True)


if __name__ == "__main__":
    main()
