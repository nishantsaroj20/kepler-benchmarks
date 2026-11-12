#!/usr/bin/env python3
"""Verifier: run the submitted column_solver.py on held-out hidden cases.

Agent code is launched in an isolated temp workdir with only that case's inputs
copied in. Hidden truth under /tests/hidden stays root-only (700) so a naive
manifest-echo cheat cannot read it. Children run in a new process group with a
minimal env and, when available, as an unprivileged user.
"""

from __future__ import annotations

import json
import math
import os
import pwd
import shutil
import signal
import site
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

PROGRAM = Path("/workspace/output/column_solver.py")
HIDDEN_ROOT = Path("/tests/hidden")
MANIFEST = HIDDEN_ROOT / "manifest.json"

_REPO = Path(__file__).resolve().parents[1]
if not PROGRAM.exists():
    local_prog = _REPO / "output" / "column_solver.py"
    if local_prog.exists():
        PROGRAM = local_prog
if not MANIFEST.exists():
    local_hidden = _REPO / "tests" / "hidden"
    if (local_hidden / "manifest.json").exists():
        HIDDEN_ROOT = local_hidden
        MANIFEST = local_hidden / "manifest.json"

REQUIRED_KEYS = ("qm", "K", "k_LDF", "D_ax", "L_min_m")


def _harden_hidden_tree() -> None:
    """Best-effort: agent children must not read sealed hidden truth."""
    if not HIDDEN_ROOT.exists():
        return
    try:
        for root, dirs, files in os.walk(HIDDEN_ROOT):
            os.chmod(root, 0o700)
            for name in files:
                os.chmod(os.path.join(root, name), 0o600)
    except OSError:
        pass


def _system_site_pythonpath() -> str:
    return os.pathsep.join(site.getsitepackages())


def _agent_user_kwargs(env: dict) -> dict:
    """Drop privileges when nobody can still import the scientific stack."""
    if os.geteuid() != 0:
        return {}
    try:
        pw = pwd.getpwnam("nobody")
    except KeyError:
        return {}
    probe = subprocess.run(
        [sys.executable, "-s", "-c", "import numpy,scipy,pandas"],
        env=env,
        user=pw.pw_uid,
        group=pw.pw_gid,
        capture_output=True,
        text=True,
        timeout=30,
    )
    if probe.returncode != 0:
        return {}
    return {"user": pw.pw_uid, "group": pw.pw_gid}


def _minimal_env(work: Path) -> dict:
    path = os.environ.get("PATH", "/usr/local/bin:/usr/bin:/bin")
    env = {
        "PATH": path,
        "HOME": str(work),
        "TMPDIR": str(work / "tmp"),
        "TMP": str(work / "tmp"),
        "TEMP": str(work / "tmp"),
        "LANG": "C.UTF-8",
        "LC_ALL": "C.UTF-8",
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
        "PYTHONPATH": _system_site_pythonpath(),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
    }
    return env


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def max_outlet_under_truth(case: dict, params: dict, L: float, n_nodes: int = 64) -> float:
    """Independent MOL check used for L_min feasibility (not the agent's code)."""
    duty = case["design_duty"]
    eps = float(case["bed_voidage"])
    rho = float(case["particle_density_kg_m3"])
    c_in = float(duty["c_in_mol_m3"])
    u = float(duty["interstitial_velocity_m_s"])
    t_service = float(duty["service_time_s"])
    qm = float(params["qm"])
    K = float(params["K"])
    k_ldf = float(params["k_LDF"])
    D_ax = float(params["D_ax"])
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho
    y0 = np.zeros(2 * n_nodes)

    def rhs(_t, y):
        c = y[:n_nodes]
        q = y[n_nodes:]
        dqdt = k_ldf * (q_star(c, qm, K) - q)
        dc = np.empty(n_nodes)
        dc[0] = u * (c[0] - c_in) / dz
        dc[1:] = u * (c[1:] - c[:-1]) / dz
        c_im1 = np.empty(n_nodes)
        c_ip1 = np.empty(n_nodes)
        c_im1[0] = c_in
        c_im1[1:] = c[:-1]
        c_ip1[:-1] = c[1:]
        c_ip1[-1] = c[-1]
        d2 = (c_ip1 - 2.0 * c + c_im1) / (dz * dz)
        return np.concatenate([-dc + D_ax * d2 - factor * dqdt, dqdt])

    sol = solve_ivp(
        rhs,
        (0.0, t_service),
        y0,
        method="BDF",
        t_eval=np.linspace(0.0, t_service, 201),
        rtol=1e-6,
        atol=1e-9,
    )
    assert sol.success, sol.message
    return float(np.max(sol.y[n_nodes - 1, :]))


def run_program_isolated(program_path: Path, case_src: Path, timeout_s: float) -> dict:
    """Copy only one case into a temp dir and execute the agent program there."""
    _harden_hidden_tree()
    with tempfile.TemporaryDirectory(prefix="adsorb_case_") as tmp:
        work = Path(tmp)
        case_dir = work / "case"
        tmp_dir = work / "tmp"
        out_path = work / "result.json"
        prog_path = work / "column_solver.py"
        shutil.copytree(case_src, case_dir)
        shutil.copy2(program_path, prog_path)
        tmp_dir.mkdir()
        # Ensure the unprivileged child can read the sandboxed case + program.
        for p in (work, case_dir, tmp_dir):
            os.chmod(p, 0o755)
        for f in case_dir.rglob("*"):
            if f.is_file():
                os.chmod(f, 0o644)
        os.chmod(prog_path, 0o755)

        cmd = [
            sys.executable,
            "-s",  # no user site-packages; system site via PYTHONPATH only
            str(prog_path),
            "--case",
            str(case_dir / "case.json"),
            "--out",
            str(out_path),
        ]
        env = _minimal_env(work)
        # Allow the current interpreter's user site only when we cannot drop
        # privileges (local host installs often live under ~/.local).
        user_kwargs = _agent_user_kwargs(env)
        if not user_kwargs:
            user_site = site.getusersitepackages()
            if user_site:
                # Prefer the invoking user's scientific stack for local runs.
                env["PYTHONPATH"] = os.pathsep.join(
                    [user_site, env.get("PYTHONPATH", "")]
                )
        else:
            try:
                uid = user_kwargs["user"]
                gid = user_kwargs["group"]
                for root, dirs, files in os.walk(work):
                    os.chown(root, uid, gid)
                    for name in files:
                        os.chown(os.path.join(root, name), uid, gid)
            except OSError:
                user_kwargs = {}

        t0 = time.monotonic()
        try:
            proc = subprocess.run(
                cmd,
                cwd=str(work),
                env=env,
                capture_output=True,
                text=True,
                timeout=timeout_s,
                start_new_session=True,
                **user_kwargs,
            )
        except subprocess.TimeoutExpired as exc:
            # Kill the whole process group if still alive.
            try:
                if exc.process is not None and exc.process.pid:
                    os.killpg(exc.process.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError, AttributeError):
                pass
            pytest.fail(
                f"case {case_src.name}: program exceeded {timeout_s}s per-run limit "
                f"(timeout fails the task); stdout={exc.stdout!r} stderr={exc.stderr!r}"
            )
        elapsed = time.monotonic() - t0
        if proc.returncode != 0:
            pytest.fail(
                f"case {case_src.name}: program exited {proc.returncode} in {elapsed:.2f}s\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )
        assert out_path.is_file(), f"case {case_src.name}: no output written"
        assert not out_path.is_symlink(), f"case {case_src.name}: symlink output rejected"
        payload = json.loads(out_path.read_text())
        assert isinstance(payload, dict), f"case {case_src.name}: output must be a JSON object"
        # Copy payload out before temp dir is destroyed.
        return {"payload": payload, "elapsed_s": elapsed}


@pytest.fixture(scope="module")
def manifest():
    _harden_hidden_tree()
    assert MANIFEST.is_file(), f"missing manifest {MANIFEST}"
    data = json.loads(MANIFEST.read_text())
    assert isinstance(data, dict)
    return data


@pytest.fixture(scope="module")
def program_path():
    assert PROGRAM.is_file(), f"missing program artifact {PROGRAM}"
    assert not PROGRAM.is_symlink(), "program must not be a symlink"
    text = PROGRAM.read_text(errors="ignore")
    assert text.strip(), "program file is empty"
    return PROGRAM


@pytest.fixture(scope="module")
def hidden_results(program_path, manifest):
    timeout_s = float(manifest["per_run_timeout_s"])
    results = {}
    for entry in manifest["hidden_cases"]:
        case_id = entry["case_id"]
        case_src = HIDDEN_ROOT / case_id
        assert (case_src / "case.json").is_file(), f"missing hidden case {case_src}"
        case = json.loads((case_src / "case.json").read_text())
        run = run_program_isolated(program_path, case_src, timeout_s)
        results[case_id] = {
            "payload": run["payload"],
            "elapsed_s": run["elapsed_s"],
            "truth": entry,
            "case": case,
        }
    return results


def _finite_positive(value, name):
    assert isinstance(value, (int, float)) and not isinstance(value, bool), (
        f"{name} must be numeric"
    )
    value = float(value)
    assert math.isfinite(value), f"{name} must be finite"
    assert value > 0.0, f"{name} must be positive"
    return value


def test_program_exists(program_path):
    assert program_path.name.endswith(".py")


def test_manifest_echo_cheat_fails(program_path, manifest):
    """A short manifest-echo cheat must not recover sealed truth."""
    timeout_s = min(30.0, float(manifest["per_run_timeout_s"]))
    cheat = '''#!/usr/bin/env python3
import json, sys
from pathlib import Path
for p in ("/tests/hidden/manifest.json", "tests/hidden/manifest.json", "../hidden/manifest.json"):
    try:
        data = json.loads(Path(p).read_text())
        entry = data["hidden_cases"][0]
        out = {"qm": entry["parameters"]["qm"], "K": entry["parameters"]["K"],
               "k_LDF": entry["parameters"]["k_LDF"], "D_ax": entry["parameters"]["D_ax"],
               "L_min_m": entry["L_min_m"]}
        Path(sys.argv[sys.argv.index("--out")+1]).write_text(json.dumps(out))
        raise SystemExit(0)
    except Exception:
        pass
raise SystemExit("no manifest")
'''
    with tempfile.TemporaryDirectory() as tmp:
        cheat_path = Path(tmp) / "cheat.py"
        cheat_path.write_text(cheat)
        case_src = HIDDEN_ROOT / manifest["hidden_cases"][0]["case_id"]
        cheated = False
        try:
            run_program_isolated(cheat_path, case_src, timeout_s)
            cheated = True
        except BaseException:
            cheated = False
        assert not cheated, "manifest-echo cheat must fail under isolation"


def test_all_hidden_cases_run(hidden_results, manifest):
    expected = {e["case_id"] for e in manifest["hidden_cases"]}
    assert set(hidden_results) == expected


def test_constant_baseline_fails_every_case(manifest):
    """Hardcoded constants (incl. visible V4-like values) must fail all hidden cases."""
    # Deliberately close to the old near-duplicate V4/H3 cluster.
    constant = {
        "qm": 0.43,
        "K": 2.7,
        "k_LDF": 0.027,
        "D_ax": 2.3e-4,
        "L_min_m": 0.12,
    }
    failures = []
    for entry in manifest["hidden_cases"]:
        case_id = entry["case_id"]
        tol = entry["tolerances"]
        params = entry["parameters"]
        ok_params = True
        for key, rel in (
            ("qm", tol["qm_rel"]),
            ("K", tol["K_rel"]),
            ("k_LDF", tol["k_LDF_rel"]),
            ("D_ax", tol["D_ax_rel"]),
        ):
            err = abs(constant[key] - params[key]) / abs(params[key])
            if err > rel:
                ok_params = False
                break
        case = json.loads((HIDDEN_ROOT / case_id / "case.json").read_text())
        c_in = float(case["design_duty"]["c_in_mol_m3"])
        c_limit = float(case["design_duty"]["outlet_mole_fraction_limit"]) * c_in
        L_ref = float(entry["L_min_m"])
        L = constant["L_min_m"]
        oversized = L > float(tol["L_min_oversize_ceil"]) * L_ref
        try:
            max_out = max_outlet_under_truth(case, params, L, n_nodes=64)
            feasible = max_out <= c_limit * float(tol["outlet_feasibility_slack"])
        except Exception:
            feasible = False
        passes = ok_params and feasible and not oversized
        if passes:
            failures.append(case_id)
    assert not failures, (
        f"constant baseline incorrectly passes cases: {failures}"
    )


@pytest.mark.parametrize(
    "case_id",
    [
        "H1_new_material",
        "H2_dispersion_regime",
        "H3_sensor_faults",
        "H4_high_cin",
        "H5_near_linear",
    ],
)
def test_hidden_case_accuracy(case_id, hidden_results, manifest):
    assert case_id in hidden_results, f"missing result for {case_id}"
    blob = hidden_results[case_id]
    payload = blob["payload"]
    truth = blob["truth"]
    case = blob["case"]
    tol = truth["tolerances"]
    for key in REQUIRED_KEYS:
        assert key in payload, f"{case_id}: missing key {key}"
        _finite_positive(payload[key], f"{case_id}.{key}")

    params = truth["parameters"]
    for key, rel in (
        ("qm", tol["qm_rel"]),
        ("K", tol["K_rel"]),
        ("k_LDF", tol["k_LDF_rel"]),
        ("D_ax", tol["D_ax_rel"]),
    ):
        est = float(payload[key])
        ref = float(params[key])
        assert abs(est - ref) / abs(ref) <= rel + 1e-12, (
            f"{case_id}: {key} rel err too large est={est} ref={ref} tol={rel}"
        )

    L_est = float(payload["L_min_m"])
    L_ref = float(truth["L_min_m"])
    assert L_est <= float(tol["L_min_oversize_ceil"]) * L_ref + 1e-12, (
        f"{case_id}: L_min trivially oversized est={L_est} ref={L_ref}"
    )
    c_in = float(case["design_duty"]["c_in_mol_m3"])
    c_limit = float(case["design_duty"]["outlet_mole_fraction_limit"]) * c_in
    slack = float(tol["outlet_feasibility_slack"])
    max_out = max_outlet_under_truth(case, params, L_est, n_nodes=64)
    assert max_out <= c_limit * slack + 1e-12, (
        f"{case_id}: true-model outlet {max_out} exceeds limit {c_limit} at L={L_est}"
    )
