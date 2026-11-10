"""Calibrate the verifier with a noise-weighted inverse and direct Newton forecasts.

Starts use the recorded independent inverse and seeded perturbations, not truth.
The verifier is consulted only after each complete fit has been reported.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.linalg import solve
from scipy.optimize import linear_sum_assignment

from independent_inverse import DirectSearch, KEYS, ROOT, load, matched_errors

sys.path.insert(0, str(ROOT / "solution"))
from solve_damage import Campaign


class LikelihoodSearch(DirectSearch):
    def residual(self, p):
        noise = self.spec["observation"]["noise"]
        frequency_std = noise["frequency_relative_std"]
        shape_std = noise["unit_shape_component_std"]
        blocks = []
        signs = []
        for force, frequency, shape in self.obs:
            best = None
            for branch in (-1, 1):
                try:
                    predicted_f, predicted_y = self.modes(force, p, branch)
                except (RuntimeError, ValueError, np.linalg.LinAlgError):
                    continue
                relative = (frequency / predicted_f[:7, None] - 1) / frequency_std
                alignment = np.clip(predicted_y.T @ shape, -1, 1)
                # Profile the arbitrary signed scale by projecting the clean
                # unit vector onto the line of the observed shape.
                cost = relative**2 + (1 - alignment**2) / shape_std**2
                rows, cols = linear_sum_assignment(cost)
                indices = rows[np.argsort(cols)]
                projection = alignment[indices, np.arange(len(frequency))]
                mismatch = predicted_y[:, indices] - shape * projection
                block = np.r_[relative[indices, np.arange(len(frequency))],
                              mismatch.ravel() / shape_std]
                score = float(block @ block)
                if best is None or score < best[0]:
                    best = score, branch, block
            if best is None:
                width = sum(len(f) + y.size for _, f, y in self.obs) + (len(p) - 2) // 6
                return np.full(width, 1e6), [-1] * len(self.obs)
            _, branch, block = best
            blocks.append(block)
            signs.append(branch)
        axis_penalty = 1e5 * np.maximum(0.0, p[3:-2:6] + 1e-6 - p[2:-2:6])
        return np.concatenate([*blocks, axis_penalty]), signs


class CachedLikelihood(LikelihoodSearch):
    def __init__(self, spec, obs, sensors):
        super().__init__(spec, obs, sensors)
        self.campaign = Campaign(spec, sensors)
        m, b = self.campaign.model, self.campaign.plate
        strain = np.concatenate([m.B, m.R], axis=2)
        mat = m.A.sum(axis=0)
        gram = np.einsum("qai,ab,qbj,q->qij", strain, mat, strain, m.w, optimize=True)
        self.membrane = gram.reshape(len(gram), -1)
        self.intact = gram.sum(axis=0)
        self.bend = np.einsum("qai,ab,qbj,q->qij", b.B, b.delta, b.B, b.wt, optimize=True).reshape(len(gram), -1)
        self.intact_bend = np.einsum("qai,ab,qbj,q->ij", b.B, b.D_intact, b.B, b.wt, optimize=True)
        self.i, self.j = np.array(m.pairs).T
        self.eye = np.eye(9)
        m.features = self.features
        self.previous = None

    def features(self, q):
        i, j, eye = self.i, self.j, self.eye
        return np.r_[q, q[i] * q[j]], np.vstack([eye, eye[i] * q[j, None] + eye[j] * q[i, None]])

    def operators_for(self, p):
        m, b = self.campaign.model, self.campaign.plate
        if self.previous is None or not np.array_equal(self.previous, p):
            sep = b.separated(p)
            springs = sum(factor * stiffness for factor, stiffness in zip([p[-2], 1, 1, p[-1]], b.lin.edge))
            bend = (sep @ self.bend).reshape(9, 9)
            loss = self.spec["damage_model"]["open_membrane_loss"]
            retain = self.spec["damage_model"]["closed_bending_retain"]
            opened = self.intact - loss * (sep @ self.membrane).reshape(self.intact.shape)
            self.operators = {}
            for s, matrix, factor in [(-1, self.intact, retain), (1, opened, 1)]:
                c, d = matrix[:18, :18], matrix[:18, 18:]
                S = (matrix[18:, 18:] - d.T @ solve(c, d, assume_a="pos")) / m.E0
                K = ((self.intact_bend + springs + factor * bend) * m.mass * m.h**2 - m.N*m.G) / m.E0
                self.operators[s] = K, S
            self.previous = p.copy()
        return self.operators

    def modes(self, force, p, branch):
        m = self.campaign.model
        m.K, m.S = self.operators_for(p)[branch]
        _, f, y = m.observe(force, branch)
        return f, y[:, :7]


def main():
    spec, obs, sensors, low, high = load()
    search = CachedLikelihood(spec, obs, sensors)
    evidence = ROOT / "authoring/evidence"
    initial = json.loads((evidence / "independent-inverse.json").read_text())["artifact"]["plates"][0]
    start = np.array([zone[key] for zone in initial["zones"] for key in KEYS]
                     + [initial["fixture"][key] for key in ("x0_multiplier", "yL_multiplier")])
    rng = np.random.default_rng(831)
    guesses = [start] + [np.clip(start + (high - low) * rng.normal(0, 0.035, len(start)), low, high)
                         for _ in range(4)]
    results = []
    for index, guess in enumerate(guesses):
        fitted, cost = search.polish(guess, low, high, np.arange(len(start)), 100, 0.002)
        print(index, "chi_square", 2 * cost, flush=True)
        results.append((cost, fitted))
    # Select solely by observation loss, before reading verifier outcomes.
    chosen = min(range(len(results)), key=lambda i: results[i][0])
    report = {"method": "Condensed forward, noise-scaled relative frequency and profiled shape projection, "
                        "with direct Newton snap and forecasts. "
                        "One public-data inverse warm start and four seeded perturbations.",
              "selected_by_lowest_residual": chosen, "fits": []}
    for index, (cost, fitted) in enumerate(results):
        _, signs = search.residual(fitted)
        snap = float(search.model.snap(fitted))
        forecast = search.model.forecast(fitted)
        plate = {
            "plate_id": spec["specimen_id"],
            "zones": [dict(zip(KEYS, map(float, row))) for row in fitted[:-2].reshape(-1, 6)],
            "fixture": dict(zip(("x0_multiplier", "yL_multiplier"), map(float, fitted[-2:]))),
            "branches": {rec["record_id"]: sign for rec, sign in zip(spec["records"], signs)},
            "snap_load_N": snap,
            "open_branch_forecast_Hz": forecast,
        }
        artifact = {"plates": [plate]}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "damage.json"
            path.write_text(json.dumps(artifact))
            test = subprocess.run(
                [sys.executable, "-m", "pytest", str(ROOT / "tests/test_outputs.py"), "-q", "--tb=short"],
                env={**os.environ, "LDI_ARTIFACT": str(path)}, capture_output=True, text=True,
                check=False,
            )
        print(index, "verifier exit", test.returncode, flush=True)
        report["fits"].append({"chi_square": 2 * cost, "verifier_exit": test.returncode,
                               "verifier_output": test.stdout, "matched_errors": matched_errors(artifact),
                               "artifact": artifact})
    (evidence / "likelihood-calibration.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
