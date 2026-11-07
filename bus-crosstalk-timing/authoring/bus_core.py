"""Coupled RC bus driven by level-1 inverters."""

from __future__ import annotations

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import splu


def _mos(vgs, vds, drive, threshold, lam):
    """Drain-to-source current and dI/dVds with the gate held fixed."""
    if vgs <= threshold:
        return 0.0, 0.0
    overdrive = vgs - threshold
    if vds < 0.0:
        current, conductance = _mos(vgs - vds, -vds, drive, threshold, lam)
        return -current, conductance
    if vds < overdrive:
        base = drive * (overdrive * vds - 0.5 * vds * vds)
        current = base * (1.0 + lam * vds)
        conductance = drive * (overdrive - vds) * (1.0 + lam * vds) + base * lam
        return current, conductance
    base = 0.5 * drive * overdrive * overdrive
    current = base * (1.0 + lam * vds)
    return current, base * lam


def _driver_current(output_v, gate_v, spec):
    vdd = spec["vdd_v"]
    n_current, n_cond = _mos(gate_v, output_v, spec["kn_a_per_v2"], spec["vt_n_v"], spec["lambda_per_v"])
    p_current, p_cond = _mos(vdd - gate_v, vdd - output_v, spec["kp_a_per_v2"], spec["vt_p_v"], spec["lambda_per_v"])
    # Current leaving the output node, and its derivative with respect to that voltage.
    return n_current - p_current, n_cond + p_cond


def _gate_voltage(wire, spec, time_s):
    start = wire["skew_s"]
    ramp = spec["ramp_s"]
    if time_s <= start:
        return wire["gate_from_v"]
    if time_s >= start + ramp:
        return wire["gate_to_v"]
    fraction = (time_s - start) / ramp
    return wire["gate_from_v"] + fraction * (wire["gate_to_v"] - wire["gate_from_v"])


def _rest_voltage(wire, spec):
    return 0.0 if wire["gate_from_v"] > 0.5 * spec["vdd_v"] else spec["vdd_v"]


def build_system(spec, n_seg, coupling_scale=1.0):
    wires = spec["wire_count"]
    nodes = n_seg + 1
    count = wires * nodes
    resistance = spec["resistance_ohm_per_m"] * spec["length_m"] / n_seg
    ground = spec["ground_cap_f_per_m"] * spec["length_m"] / n_seg
    couple = coupling_scale * spec["coupling_cap_f_per_m"] * spec["length_m"] / n_seg
    conductance = lil_matrix((count, count))
    capacitance = lil_matrix((count, count))

    def index(wire, node):
        return wire * nodes + node

    for wire in range(wires):
        for node in range(n_seg):
            left = index(wire, node)
            right = index(wire, node + 1)
            value = 1.0 / resistance
            conductance[left, left] += value
            conductance[right, right] += value
            conductance[left, right] -= value
            conductance[right, left] -= value
        for node in range(nodes):
            share = 0.5 if node in (0, n_seg) else 1.0
            here = index(wire, node)
            capacitance[here, here] += share * ground
            if wire + 1 < wires:
                other = index(wire + 1, node)
                capacitance[here, here] += share * couple
                capacitance[other, other] += share * couple
                capacitance[here, other] -= share * couple
                capacitance[other, here] -= share * couple
        capacitance[index(wire, n_seg), index(wire, n_seg)] += spec["receiver_cap_f"]
    return conductance.tocsc(), capacitance.tocsc(), nodes


def simulate(spec, n_seg, dt_s, coupling_scale=1.0):
    conductance, capacitance, nodes = build_system(spec, n_seg, coupling_scale)
    wires = spec["wire_count"]
    count = wires * nodes
    voltage = np.zeros(count)
    for wire, item in enumerate(spec["wires"]):
        voltage[wire * nodes : (wire + 1) * nodes] = _rest_voltage(item, spec)
    over_dt = capacitance / dt_s
    horizon = spec["horizon_s"]
    steps = int(np.ceil(horizon / dt_s))
    times = np.zeros(steps + 1)
    far = np.zeros((steps + 1, wires))
    for wire in range(wires):
        far[0, wire] = voltage[wire * nodes + n_seg]
    previous = voltage.copy()
    for step in range(1, steps + 1):
        time_s = step * dt_s
        times[step] = time_s
        guess = previous.copy()
        converged = False
        for _ in range(16):
            residual_driver = np.zeros(count)
            jacobian = (over_dt + conductance).tolil()
            for wire, item in enumerate(spec["wires"]):
                node = wire * nodes
                current, slope = _driver_current(guess[node], _gate_voltage(item, spec, time_s), spec)
                residual_driver[node] = current
                jacobian[node, node] += slope
            residual = over_dt @ (guess - previous) + conductance @ guess + residual_driver
            if np.max(np.abs(residual)) < 1e-9:
                converged = True
                break
            delta = splu(jacobian.tocsc()).solve(-residual)
            guess = guess + delta
            if np.max(np.abs(delta)) < 1e-9 and np.max(np.abs(residual)) < 1e-8:
                converged = True
                break
        if not converged:
            raise RuntimeError(f"Newton did not settle at t = {time_s}")
        previous = guess
        for wire in range(wires):
            far[step, wire] = previous[wire * nodes + n_seg]
    return times, far


def crossing_times_s(spec, times, far):
    threshold = 0.5 * spec["vdd_v"]
    crossings = []
    for wire, item in enumerate(spec["wires"]):
        wave = far[:, wire]
        falling = _rest_voltage(item, spec) > threshold
        crossed = None
        for index in range(1, len(times)):
            before, after = wave[index - 1], wave[index]
            hit = before >= threshold > after if falling else before <= threshold < after
            if hit:
                span = after - before
                fraction = 0.0 if span == 0.0 else (threshold - before) / span
                crossed = times[index - 1] + fraction * (times[index] - times[index - 1])
                break
        if crossed is None:
            raise RuntimeError(f"wire {wire} did not cross half Vdd")
        crossings.append(float(crossed))
    return crossings


def limit_clock_mhz(spec, crossings_s):
    period = max(crossings_s) + spec["setup_time_s"]
    return 1e-6 / period


def evaluate(spec, n_seg, dt_s, coupling_scale=1.0):
    times, far = simulate(spec, n_seg, dt_s, coupling_scale)
    crossings = crossing_times_s(spec, times, far)
    return {
        "far_end_cross_ps": [value * 1e12 for value in crossings],
        "limit_clock_mhz": limit_clock_mhz(spec, crossings),
        "n_seg": n_seg,
        "dt_ps": dt_s * 1e12,
    }
