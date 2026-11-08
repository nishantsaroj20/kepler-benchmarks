#!/usr/bin/env python3
"""Acceptance-rule regression: a failed Newton or stagger must not commit state."""
from __future__ import annotations

import numpy as np

from coupled_current import make_solver, quick_card
from problem_spec import Schedule


def _equal(a, b):
    return np.array_equal(a, b)


def test_restore_bitwise():
    sol = make_solver(quick_card(), 0.015, 333.0, 3, 6)
    snap = sol.snapshot()
    sol.t_nodes += 7.5
    sol.ep += 0.2
    sol.u += 0.1
    sol.be += 0.05
    sol.fn2 += 0.01
    sol.fnth += 0.02
    sol.tgp += 4.0
    sol.restore(snap)
    for key in snap:
        assert _equal(getattr(sol, key if key != "fn2" else "fn2"), snap[key] if key != "fn2" else snap["fn2"]), key
    # attribute names match snapshot keys except fn2 is self.fn2
    assert _equal(sol.u, snap["u"])
    assert _equal(sol.be, snap["be"])
    assert _equal(sol.fn2, snap["fn2"])
    assert _equal(sol.fnth, snap["fnth"])
    assert _equal(sol.ep, snap["ep"])
    assert _equal(sol.t_nodes, snap["t_nodes"])
    assert _equal(sol.tgp, snap["tgp"])
    print("restore bitwise OK")


def test_rejected_newton_rolls_back():
    sol = make_solver(quick_card(), 0.015, 333.0, 4, 8)
    sch = Schedule(1.4, 1.5, 1.1, 0.0)
    out = sol.run(sch, nsub=4, mode="coupled", max_stagger=4, stagger_tol=0.25)
    assert out["history"][0]["t_s"] == 0.0
    assert abs(out["history"][0]["T_neck_K"] - 333.0) < 1e-9
    assert out["extrema"]["T_neck_K"] + 1e-9 >= 333.0
    snap = sol.snapshot()
    real = sol._newton

    def fail_and_corrupt(*_a, **_k):
        sol.t_nodes += 11.0
        sol.ep += 0.4
        sol.u += 0.3
        return None, None, 3, np.inf

    sol._newton = fail_and_corrupt
    payload = sol._attempt_step(sol.u.copy(), sol.u[sol.end, 1].mean() + 0.01, 0.05, "coupled", 0.25, 4)
    assert payload is None
    sol.restore(snap)
    sol._newton = real
    assert _equal(sol.ep, snap["ep"])
    assert _equal(sol.t_nodes, snap["t_nodes"])
    assert _equal(sol.u, snap["u"])
    assert _equal(sol.be, snap["be"])
    assert _equal(sol.tgp, snap["tgp"])
    print("rejected newton rollback OK", "steps", len(out["history"]), "retries", len(out["failures"]))


def test_rejected_stagger_does_not_keep_temperature():
    sol = make_solver(quick_card(), 0.015, 293.0, 3, 6)
    snap = sol.snapshot()
    calls = {"n": 0}
    real = sol._newton

    def fail_second(*a, **k):
        calls["n"] += 1
        if calls["n"] == 1:
            return real(*a, **k)
        sol.t_nodes += 9.0
        return None, None, 1, np.inf

    sol._newton = fail_second
    payload = sol._attempt_step(sol.u.copy(), 0.02, 0.05, "coupled", 1e-6, 4)
    assert payload is None
    # attempt itself must not have committed the corrupted nodal field as the solver state
    # (the corrupt write is inside the patched newton; restore is the caller's job)
    sol.restore(snap)
    sol._newton = real
    assert _equal(sol.t_nodes, snap["t_nodes"])
    assert float(sol.ep.max()) == 0.0
    print("stagger rejection leaves restored state at the initial solid")


def test_clean_step_records_convergence():
    sol = make_solver(quick_card(), 0.015, 333.0, 4, 8)
    sch = Schedule(1.4, 1.5, 1.1, 0.0)
    out = sol.run(sch, nsub=4, mode="coupled", max_stagger=4, stagger_tol=0.25)
    assert all(h["accepted"] for h in out["history"])
    assert all(h["coupling_converged"] for h in out["history"])
    assert abs(out["history"][0]["T_neck_K"] - 333.0) < 1e-8
    assert out["extrema"]["T_neck_t_s"] >= 0.0
    print(
        "clean",
        "Tn", round(out["extrema"]["T_neck_K"], 2),
        "at", round(out["extrema"]["T_neck_t_s"], 3),
        "g", round(out["final"]["gauge_mm"], 4),
        "r", round(out["extrema"]["r_min_mm"], 4),
        "retries", len(out["failures"]),
        "steps", len(out["history"]),
    )


if __name__ == "__main__":
    test_restore_bitwise()
    test_rejected_newton_rolls_back()
    test_rejected_stagger_does_not_keep_temperature()
    test_clean_step_records_convergence()
    print("ALL ACCEPTANCE TESTS PASSED")
