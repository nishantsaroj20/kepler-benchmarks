#!/bin/sh
set -eu
mkdir -p /app /app/output
cat > /app/forecast.py << 'SLURRY_FORECAST_PY_EOF'
#!/usr/bin/env python3
"""Forecast slurry hydrogen rates for a campaign directory.

The historical rates identify a saturating surface rate. Uniform particles
integrate that rate against steady diffusion; a closed headspace is then
integrated over the requested hold.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
from scipy.optimize import least_squares

R_GAS = 8.314462618
T_REF = 348.15
NODES = 64
HOLD_STEPS = 120
INITIAL_OUTER = 8
INITIAL_INNER = 4


def thomas(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    width = diag.shape[1]
    cp = np.empty_like(diag)
    dp = np.empty_like(rhs)
    inv = 1.0 / diag[:, 0]
    cp[:, 0] = upper[:, 0] * inv
    dp[:, 0] = rhs[:, 0] * inv
    for index in range(1, width):
        den = diag[:, index] - lower[:, index] * cp[:, index - 1]
        inv = 1.0 / den
        cp[:, index] = upper[:, index] * inv
        dp[:, index] = (rhs[:, index] - lower[:, index] * dp[:, index - 1]) * inv
    solution = np.empty_like(rhs)
    solution[:, -1] = dp[:, -1]
    for index in range(width - 2, -1, -1):
        solution[:, index] = dp[:, index] - cp[:, index] * solution[:, index + 1]
    return solution


def sphere_geometry(radius: float, nodes: int, shell_thickness: float | None):
    if shell_thickness is None:
        inner_radius = 0.0
        span = radius
    else:
        inner_radius = radius - shell_thickness
        span = shell_thickness
    dr = span / nodes
    index = np.arange(nodes)
    r_in = inner_radius + index * dr
    r_out = inner_radius + (index + 1) * dr
    volume = (4.0 / 3.0) * np.pi * (r_out**3 - r_in**3)
    area_out = 4.0 * np.pi * r_out**2
    area_in = 4.0 * np.pi * r_in**2
    node = 0.5 * (r_in + r_out)
    gap = np.empty(nodes + 1)
    gap[0] = node[0] - inner_radius
    gap[1:-1] = np.diff(node)
    gap[-1] = radius - node[-1]
    g_out = area_out / gap[1:]
    g_in = np.zeros(nodes)
    if shell_thickness is None and nodes > 1:
        g_in[1:] = area_in[1:] / gap[1:-1]
    elif shell_thickness is not None and nodes > 1:
        g_in[1:] = area_in[1:] / gap[1:-1]
    particle = (4.0 / 3.0) * np.pi * radius**3
    if shell_thickness is None:
        activity_scale = 1.0
    else:
        fraction = 1.0 - (inner_radius / radius) ** 3
        activity_scale = 1.0 / fraction
    return volume, g_in, g_out, particle, activity_scale


def solve_profiles(conc, surface, diffusivity, rho, activity, k_value, k_ads, volume, g_in, g_out, scale, rounds):
    g_out_b = diffusivity * g_out
    g_in_b = diffusivity * g_in
    for _ in range(rounds):
        theta = k_ads[:, None] * conc
        rate = rho * scale * activity[:, None] * k_value[:, None] * theta / (1.0 + theta)
        slope = rho * scale * activity[:, None] * k_value[:, None] * k_ads[:, None] / (1.0 + theta) ** 2
        outer = np.empty_like(conc)
        outer[:, :-1] = conc[:, 1:]
        outer[:, -1] = surface
        inner = np.zeros_like(conc)
        inner[:, 1:] = conc[:, :-1]
        supply = g_in_b * (inner - conc) + g_out_b * (outer - conc)
        residual = supply - rate * volume
        diag = -(g_in_b + g_out_b) - slope * volume
        lower = np.broadcast_to(g_in_b, conc.shape).copy()
        upper = np.broadcast_to(g_out_b, conc.shape).copy()
        upper[:, -1] = 0.0
        conc = np.maximum(conc + thomas(lower, diag, upper, -residual), 0.0)
    return conc


def catalytic_rate(conc, surface, diffusivity, g_out, particle, rho, loading):
    influx = diffusivity * g_out[-1] * (surface - conc[:, -1])
    return loading * influx / (particle * rho)


def initial_rates(pressure, temperature, kla, ks, loading, activity, k_value, k_ads, c_one, radius, diffusivity, rho, shell, nodes):
    volume, g_in, g_out, particle, scale = sphere_geometry(radius, nodes, shell)
    cstar = c_one * pressure
    resistance = 1.0 / kla + 1.0 / np.maximum(ks, 1e-12)
    surface = np.maximum(0.35 * cstar, 1e-10)
    conc = np.broadcast_to(surface[:, None], (surface.size, nodes)).copy()
    for _ in range(INITIAL_OUTER):
        conc = solve_profiles(conc, surface, diffusivity, rho, activity, k_value, k_ads, volume, g_in, g_out, scale, INITIAL_INNER)
        cat = catalytic_rate(conc, surface, diffusivity, g_out, particle, rho, loading)
        delta = np.maximum(2e-5 * np.maximum(surface, 1e-8), 1e-12)
        bumped = solve_profiles(conc, surface + delta, diffusivity, rho, activity, k_value, k_ads, volume, g_in, g_out, scale, 2)
        cat_up = catalytic_rate(bumped, surface + delta, diffusivity, g_out, particle, rho, loading)
        slope = (cat_up - cat) / delta
        film = (cstar - surface) / resistance
        surface = np.clip(surface - (film - cat) / (-1.0 / resistance - slope), 0.0, np.maximum(cstar, 0.0))
    conc = solve_profiles(conc, surface, diffusivity, rho, activity, k_value, k_ads, volume, g_in, g_out, scale, INITIAL_INNER + 1)
    rate = catalytic_rate(conc, surface, diffusivity, g_out, particle, rho, loading)
    bulk = cstar - rate / kla
    return rate, (cstar - bulk) / cstar, (bulk - surface) / cstar, surface / cstar, conc, surface, resistance, volume, g_in, g_out, particle, scale


def integrate_hold(pressure, temperature, kla, ks, loading, activity, k_value, k_ads, c_one, radius, diffusivity, rho, shell, slurry_volume, hold, head, nodes, steps):
    rate0, gas_liquid, liquid_solid, kinetic, conc, surface, resistance, volume, g_in, g_out, particle, scale = initial_rates(
        pressure, temperature, kla, ks, loading, activity, k_value, k_ads, c_one, radius, diffusivity, rho, shell, nodes
    )
    dt = hold / steps
    partial = np.array(pressure, dtype=float, copy=True)
    consumed = np.zeros_like(partial)
    for _ in range(steps):
        cstar = c_one * np.maximum(partial, 0.0)
        for _correction in range(2):
            conc = solve_profiles(conc, surface, diffusivity, rho, activity, k_value, k_ads, volume, g_in, g_out, scale, 2)
            cat = catalytic_rate(conc, surface, diffusivity, g_out, particle, rho, loading)
            delta = np.maximum(2e-5 * np.maximum(surface, 1e-8), 1e-12)
            bumped = solve_profiles(conc, surface + delta, diffusivity, rho, activity, k_value, k_ads, volume, g_in, g_out, scale, 2)
            cat_up = catalytic_rate(bumped, surface + delta, diffusivity, g_out, particle, rho, loading)
            slope = (cat_up - cat) / delta
            film = (cstar - surface) / resistance
            surface = np.clip(surface - (film - cat) / (-1.0 / resistance - slope), 0.0, np.maximum(cstar, 0.0))
        conc = solve_profiles(conc, surface, diffusivity, rho, activity, k_value, k_ads, volume, g_in, g_out, scale, 2)
        rate = catalytic_rate(conc, surface, diffusivity, g_out, particle, rho, loading)
        rate = np.where(partial <= 0.0, 0.0, rate)
        consumed += rate * dt
        partial = np.maximum(partial - rate * slurry_volume * R_GAS * temperature / (head * 1e5) * dt, 0.0)
    return rate0, gas_liquid, liquid_solid, kinetic, consumed


def _read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _interp_linear(x_values, y_values, x_query):
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    query = np.asarray(x_query, dtype=float)
    order = np.argsort(x_values)
    x_values = x_values[order]
    y_values = y_values[order]
    index = np.searchsorted(x_values, query, side="right") - 1
    index = np.clip(index, 0, len(x_values) - 2)
    x0 = x_values[index]
    x1 = x_values[index + 1]
    y0 = y_values[index]
    y1 = y_values[index + 1]
    weight = (query - x0) / (x1 - x0)
    return y0 + weight * (y1 - y0)


def load_campaign(campaign: Path):
    vessel = json.loads((campaign / "vessel.json").read_text())
    solvent = json.loads((campaign / "solvent.json").read_text())
    catalyst = json.loads((campaign / "catalyst.json").read_text())
    solubility = _read_csv(campaign / "h2_solubility.csv")
    calibration = _read_csv(campaign / "thermocouple_calibration.csv")
    standards = _read_csv(campaign / "gc_standards.csv")
    kla_rows = _read_csv(campaign / "kla_measurements.csv")
    kc_rows = _read_csv(campaign / "kc_measurements.csv")
    site_rows = _read_csv(campaign / "catalyst_sites.csv")
    runs = _read_csv(campaign / "rate_runs.csv")

    kla = {}
    for row in kla_rows:
        if not row["N_rps"].strip():
            continue
        if "baffles removed" in row.get("note", ""):
            continue
        kla[round(float(row["N_rps"]), 6)] = float(row["kla_per_s"])
    kc = {}
    for row in kc_rows:
        if not row["N_rps"].strip():
            continue
        if "baffles removed" in row.get("note", ""):
            continue
        kc[(round(float(row["N_rps"]), 6), round(float(row["dp_um"]), 6))] = float(row["kc_m_per_s"])
    sites: dict[str, list[tuple[float | None, float]]] = {}
    for row in site_rows:
        state = row["catalyst_state"].strip()
        speed = float(row["N_rps"]) if row["N_rps"].strip() else None
        sites.setdefault(state, []).append((speed, float(row["co_uptake_umol_per_g"])))
    return {
        "slurry_volume": float(vessel["slurry_volume_m3"]),
        "p_atm": float(vessel["site_atmospheric_pressure_bar"]),
        "antoine": (float(solvent["A"]), float(solvent["B"]), float(solvent["C"])),
        "rho": float(catalyst["particle_density_kg_m3"]),
        "fresh": float(catalyst["fresh_co_uptake_umol_per_g"]),
        "distribution": catalyst["metal_location"].strip(),
        "shell_um": float(catalyst["eggshell_depth_um"]) if str(catalyst.get("eggshell_depth_um", "")).strip() not in {"", "None"} else None,
        "diffusivity": float(catalyst["effective_diffusivity_m2_per_s"]),
        "sol_t": np.array([float(row["reference_C"]) for row in solubility]),
        "sol_c": np.array([float(row["cstar_mol_m3_at_1bar"]) for row in solubility]),
        "cal_tic": np.array([float(row["tic201_C"]) for row in calibration]),
        "cal_ref": np.array([float(row["reference_C"]) for row in calibration]),
        "standards": [(row["start_date"], row["end_date"], float(row["response_factor_mol_per_count"])) for row in standards],
        "kla": kla,
        "kc": kc,
        "sites": sites,
        "runs": runs,
    }


def vapor_pressure(antoine, temperature_k):
    a_coef, b_coef, c_coef = antoine
    return 10 ** (a_coef - b_coef / (temperature_k + c_coef))


def solubility_1bar(campaign, temperature_c):
    inverse = 1.0 / (np.asarray(temperature_c, dtype=float) + 273.15)
    known = 1.0 / (campaign["sol_t"] + 273.15)
    log_c = np.log(campaign["sol_c"])
    order = np.argsort(known)
    known = known[order]
    log_c = log_c[order]
    # Reciprocal temperature falls as Celsius rises, so search the ascending inverse table.
    ascending = known[0] < known[-1]
    if not ascending:
        known = known[::-1]
        log_c = log_c[::-1]
    query = np.asarray(inverse, dtype=float)
    index = np.searchsorted(known, query, side="right") - 1
    index = np.clip(index, 0, len(known) - 2)
    weight = (query - known[index]) / (known[index + 1] - known[index])
    return np.exp(log_c[index] + weight * (log_c[index + 1] - log_c[index]))


def reference_temperature(campaign, tic):
    return _interp_linear(campaign["cal_tic"], campaign["cal_ref"], tic)


def response_factor(campaign, run_date: str) -> float:
    for start, end, factor in campaign["standards"]:
        if start <= run_date <= end:
            return factor
    raise KeyError(run_date)


def uptake_for(campaign, state: str, speed: float) -> float:
    rows = campaign["sites"][state]
    for listed_speed, uptake in rows:
        if listed_speed is not None and abs(listed_speed - speed) < 1e-6:
            return uptake
    for listed_speed, uptake in rows:
        if listed_speed is None:
            return uptake
    raise KeyError(state)


def observed_runs(campaign):
    temperatures = []
    pressures = []
    speeds = []
    loadings = []
    activities = []
    diameters = []
    kla = []
    ks = []
    rates = []
    for row in campaign["runs"]:
        if row.get("status", "ok").strip() != "ok":
            continue
        tic = float(row["tic201_C"])
        temperature_c = float(reference_temperature(campaign, np.array([tic]))[0])
        temperature_k = temperature_c + 273.15
        gauge = float(row["P_total_barg"])
        speed = float(row["N_rps"])
        mass_g = float(row["m_cat_g"])
        diameter = float(row["dp_um"])
        partial = gauge + campaign["p_atm"] - float(vapor_pressure(campaign["antoine"], temperature_k))
        loading = (mass_g / 1000.0) / campaign["slurry_volume"]
        activity = uptake_for(campaign, "conditioned", speed) / campaign["fresh"]
        area = 6.0 * loading / (campaign["rho"] * diameter * 1e-6)
        factor = response_factor(campaign, row["run_date"])
        rate = float(row["h2_peak_area"]) * factor / (float(row["window_s"]) * campaign["slurry_volume"])
        temperatures.append(temperature_k)
        pressures.append(partial)
        speeds.append(speed)
        loadings.append(loading)
        activities.append(activity)
        diameters.append(diameter)
        kla.append(campaign["kla"][round(speed, 6)])
        ks.append(campaign["kc"][(round(speed, 6), round(diameter, 6))] * area)
        rates.append(rate)
    return {key: np.array(value, dtype=float) for key, value in {
        "temperature": temperatures,
        "pressure": pressures,
        "speed": speeds,
        "loading": loadings,
        "activity": activities,
        "diameter": diameters,
        "kla": kla,
        "ks": ks,
        "rate": rates,
    }.items()}


def k_of(temperature, k_ref, energy):
    return k_ref * np.exp(-energy / R_GAS * (1.0 / temperature - 1.0 / T_REF))


def kads_of(temperature, k_ref, enthalpy):
    return k_ref * np.exp(-enthalpy / R_GAS * (1.0 / temperature - 1.0 / T_REF))


def rates_for_params(campaign, bundle, k_ref, energy, kads_ref, enthalpy, nodes):
    shell = None
    if campaign["distribution"] == "eggshell":
        shell = campaign["shell_um"] * 1e-6
    c_one = solubility_1bar(campaign, bundle["temperature"] - 273.15)
    k_value = k_of(bundle["temperature"], k_ref, energy)
    k_ads = kads_of(bundle["temperature"], kads_ref, enthalpy)
    predicted = np.zeros(bundle["rate"].shape)
    for diameter in sorted(set(np.round(bundle["diameter"], 6))):
        chosen = np.where(np.abs(bundle["diameter"] - diameter) < 1e-6)[0]
        rate, *_ = initial_rates(
            bundle["pressure"][chosen],
            bundle["temperature"][chosen],
            bundle["kla"][chosen],
            bundle["ks"][chosen],
            bundle["loading"][chosen],
            bundle["activity"][chosen],
            k_value[chosen],
            k_ads[chosen],
            c_one[chosen],
            float(diameter) * 0.5e-6,
            campaign["diffusivity"],
            campaign["rho"],
            shell,
            nodes,
        )
        predicted[chosen] = rate
    return predicted


def fit_parameters(campaign, bundle):
    measured = bundle["rate"]

    def residuals(theta):
        log_k, energy_kj, log_kads, enthalpy_kj = theta
        predicted = rates_for_params(
            campaign,
            bundle,
            math.exp(log_k),
            energy_kj * 1000.0,
            math.exp(log_kads),
            enthalpy_kj * 1000.0,
            NODES,
        )
        return (predicted - measured) / np.maximum(measured, 1e-8)

    fit = least_squares(
        residuals,
        np.array([math.log(0.01), 55.0, math.log(0.3), -35.0]),
        bounds=([-9.2, 15.0, -4.6, -140.0], [0.0, 140.0, 3.0, -2.0]),
        method="trf",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
        max_nfev=50,
    )
    log_k, energy_kj, log_kads, enthalpy_kj = fit.x
    return math.exp(float(log_k)), float(energy_kj) * 1000.0, math.exp(float(log_kads)), float(enthalpy_kj) * 1000.0


def load_cases(path: Path):
    rows = _read_csv(path)
    return rows


def forecast_cases(campaign, rows, params, nodes=NODES, steps=HOLD_STEPS):
    k_ref, energy, kads_ref, enthalpy = params
    shell = None
    if campaign["distribution"] == "eggshell":
        shell = campaign["shell_um"] * 1e-6
    temperature_c = np.array([float(row["bath_C"]) for row in rows])
    temperature = temperature_c + 273.15
    gauge = np.array([float(row["P_total_barg"]) for row in rows])
    speed = np.array([float(row["N_rps"]) for row in rows])
    mass_g = np.array([float(row["m_cat_g"]) for row in rows])
    diameter = np.array([float(row["dp_um"]) for row in rows])
    hold = np.array([float(row["hold_s"]) for row in rows])
    head = np.array([float(row["headspace_m3"]) for row in rows])
    state = [row["catalyst_state"].strip() for row in rows]
    partial = gauge + campaign["p_atm"] - vapor_pressure(campaign["antoine"], temperature)
    loading = (mass_g / 1000.0) / campaign["slurry_volume"]
    activity = np.array([uptake_for(campaign, name, float(n_value)) / campaign["fresh"] for name, n_value in zip(state, speed)])
    kla = np.array([campaign["kla"][round(float(n_value), 6)] for n_value in speed])
    ks = np.array([
        campaign["kc"][(round(float(n_value), 6), round(float(dp_value), 6))] * 6.0 * load / (campaign["rho"] * dp_value * 1e-6)
        for n_value, dp_value, load in zip(speed, diameter, loading)
    ])
    c_one = solubility_1bar(campaign, temperature_c)
    k_value = k_of(temperature, k_ref, energy)
    k_ads = kads_of(temperature, kads_ref, enthalpy)
    output = {
        "rate": np.zeros(len(rows)),
        "gas_liquid": np.zeros(len(rows)),
        "liquid_solid": np.zeros(len(rows)),
        "kinetic": np.zeros(len(rows)),
        "consumed": np.zeros(len(rows)),
    }
    for size in sorted(set(np.round(diameter, 6))):
        chosen = np.where(np.abs(diameter - size) < 1e-6)[0]
        rate, gas_liquid, liquid_solid, kinetic, consumed = integrate_hold(
            partial[chosen],
            temperature[chosen],
            kla[chosen],
            ks[chosen],
            loading[chosen],
            activity[chosen],
            k_value[chosen],
            k_ads[chosen],
            c_one[chosen],
            float(size) * 0.5e-6,
            campaign["diffusivity"],
            campaign["rho"],
            shell,
            campaign["slurry_volume"],
            hold[chosen],
            head[chosen],
            nodes,
            steps,
        )
        output["rate"][chosen] = rate
        output["gas_liquid"][chosen] = gas_liquid
        output["liquid_solid"][chosen] = liquid_solid
        output["kinetic"][chosen] = kinetic
        output["consumed"][chosen] = consumed
    packed = []
    for index, row in enumerate(rows):
        packed.append({
            "case_id": row["case_id"],
            "rate_mol_m3_s": float(output["rate"][index]),
            "hydrogen_consumed_mol_m3": float(output["consumed"][index]),
            "resistance_fraction_gas_liquid": float(output["gas_liquid"][index]),
            "resistance_fraction_liquid_solid": float(output["liquid_solid"][index]),
            "resistance_fraction_kinetic": float(output["kinetic"][index]),
        })
    return packed


def write_predictions(path: Path, cases):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"cases": cases}, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--cases", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    campaign = load_campaign(Path(args.campaign))
    bundle = observed_runs(campaign)
    params = fit_parameters(campaign, bundle)
    rows = load_cases(Path(args.cases))
    write_predictions(Path(args.output), forecast_cases(campaign, rows, params))


if __name__ == "__main__":
    main()
SLURRY_FORECAST_PY_EOF
chmod a+r /app/forecast.py
test -s /app/forecast.py
