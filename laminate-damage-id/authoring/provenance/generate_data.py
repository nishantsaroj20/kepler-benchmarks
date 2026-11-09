"""Forward records from authored parameters. Does not read an inverse result."""
import csv
import json
import sys
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "solution"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from independent import graded_forward
from solve_damage import KEYS, Campaign, artifact

spec = json.loads((root / "environment/data/plate_specs.json").read_text())
authored = json.loads(Path(__file__).with_name("parameters.json").read_text())
plate = authored["plate"]
p = np.array([z[k] for z in plate["zones"] for k in KEYS] + list(plate["fixture"].values()), float)
pts = np.array(authored["sensor_coordinates_m"], float)
data = root / "environment" / "data"
camp = Campaign(spec, pts)
rng = np.random.default_rng(authored["generation_seed"])
with (data / "sensors.csv").open("w") as fh:
    w = csv.writer(fh)
    w.writerow(["sensor_id", "x_m", "y_m"])
    w.writerows([(i + 1, *pt) for i, pt in enumerate(pts)])
for rec in spec["records"]:
    branch = int(plate["branches"][rec["record_id"]])
    _q, freq, shapes = camp.observe(rec["force_N"], p, branch)
    idx = rng.choice(np.arange(7), 4, replace=False)
    rng.shuffle(idx)
    noisy_f = freq[idx] * (1 + rng.normal(0, spec["observation"]["noise"]["frequency_relative_std"], 4))
    noisy_y = shapes[:, idx] + rng.normal(0, spec["observation"]["noise"]["unit_shape_component_std"], (len(pts), 4))
    noisy_y *= rng.choice([-1.0, 1.0], 4) * rng.uniform(0.6, 1.8, 4)
    folder = data / rec["record_id"]
    folder.mkdir(exist_ok=True)
    with (folder / "frequencies.csv").open("w") as fh:
        w = csv.writer(fh)
        w.writerow(["resonance_id", "frequency_Hz"])
        w.writerows([(i + 1, format(v, ".14g")) for i, v in enumerate(noisy_f)])
    with (folder / "mode_shapes.csv").open("w") as fh:
        w = csv.writer(fh)
        w.writerow(["sensor_id"] + [f"resonance_{i+1}" for i in range(4)])
        w.writerows([(i + 1, *[format(v, ".14g") for v in row]) for i, row in enumerate(noisy_y)])
    print(rec["record_id"], "branch", branch, "modes", idx.tolist(), flush=True)
snap, forecast = graded_forward(spec, p, quad=tuple(spec["model"]["gauss_orders"]))
signs = [plate["branches"][rec["record_id"]] for rec in spec["records"]]
body = artifact(spec, p, signs, snap, forecast)
truth = {
    "plates": body["plates"],
    "tolerances": {
        "snap_relative": 0.0024,
        "forecast_relative": 0.0024,
        "zone_centre_m": 0.026,
        "zone_a_m": 0.022,
        "zone_b_m": 0.010,
        "zone_theta_deg": 90.0,
        "zone_d": 0.40,
        "fixture_absolute": 0.07,
    },
}
(root / "tests" / "ground_truth.json").write_text(json.dumps(truth, indent=2) + "\n")
print("snap", snap, flush=True)
