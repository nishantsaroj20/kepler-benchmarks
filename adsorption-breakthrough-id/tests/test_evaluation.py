#!/usr/bin/env python3
"""Verifier: run the submitted column_solver.py on held-out hidden cases.

Agent code is launched in an isolated temp workdir with only that case's inputs
copied in. Hidden truth under /tests/hidden stays root-only (700) so a naive
manifest-echo cheat cannot read it. Children run in a new process group with a
minimal env and, when root, as an unprivileged user (no silent root fallback).
Process groups are always cleaned up after each run.
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
VISIBLE_TRUTH = Path("/tests/../authoring/provenance/visible_truth.json")

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
_LOCAL_VISIBLE = _REPO / "authoring" / "provenance" / "visible_truth.json"
if _LOCAL_VISIBLE.exists():
    VISIBLE_TRUTH = _LOCAL_VISIBLE

REQUIRED_KEYS = ("qm", "K", "k_LDF", "D_ax", "L_min_m")
EXPECTED_TEST_NODE_COUNT = 10  # used by test.sh CTRF accounting documentation


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


def _require_unprivileged_agent(env: dict) -> dict:
    """When verifier is root, agent children MUST drop to nobody.

    Silent root fallback is disallowed: if nobody cannot import the scientific
    stack, fail closed rather than executing agent code as root.
    """
    if os.geteuid() != 0:
        return {}
    try:
        pw = pwd.getpwnam("nobody")
    except KeyError as exc:
        raise RuntimeError(
            "verifier is root but user 'nobody' is missing; refusing to run agent code as root"
        ) from exc
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
        raise RuntimeError(
            "nobody cannot import numpy/scipy/pandas; refusing silent root fallback for agent code. "
            f"stderr={probe.stderr!r}"
        )
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


def _killpg_quiet(pid: int) -> None:
    try:
        os.killpg(pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError, OSError):
        pass


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def max_outlet_under_truth(case: dict, params: dict, L: float, n_nodes: int = 64) -> float:
    """Independent MOL check used for L_min feasibility (BDF, N=64)."""
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


def answer_passes_gates(payload: dict, entry: dict, case: dict) -> tuple[bool, dict]:
    """Same accuracy + BDF/N=64 feasibility gates as test_hidden_case_accuracy."""
    tol = entry["tolerances"]
    params = entry["parameters"]
    detail = {"relative_errors": {}, "feasible": False, "not_oversized": False}
    for key in REQUIRED_KEYS:
        if key not in payload:
            return False, detail
        try:
            val = float(payload[key])
        except (TypeError, ValueError):
            return False, detail
        if not math.isfinite(val) or val <= 0.0:
            return False, detail
    ok_params = True
    for key, rel in (
        ("qm", tol["qm_rel"]),
        ("K", tol["K_rel"]),
        ("k_LDF", tol["k_LDF_rel"]),
        ("D_ax", tol["D_ax_rel"]),
    ):
        err = abs(float(payload[key]) - float(params[key])) / abs(float(params[key]))
        detail["relative_errors"][key] = err
        if err > rel:
            ok_params = False
    L_est = float(payload["L_min_m"])
    L_ref = float(entry["L_min_m"])
    not_huge = L_est <= float(tol["L_min_oversize_ceil"]) * L_ref
    detail["not_oversized"] = not_huge
    c_in = float(case["design_duty"]["c_in_mol_m3"])
    c_limit = float(case["design_duty"]["outlet_mole_fraction_limit"]) * c_in
    slack = float(tol["outlet_feasibility_slack"])
    try:
        max_out = max_outlet_under_truth(case, params, L_est, n_nodes=64)
        feasible = max_out <= c_limit * slack
        detail["max_outlet"] = max_out
        detail["c_limit"] = c_limit
    except Exception as exc:
        feasible = False
        detail["feasibility_error"] = str(exc)
    detail["feasible"] = feasible
    return bool(ok_params and feasible and not_huge), detail


def run_program_isolated(program_path: Path, case_src: Path, timeout_s: float) -> dict:
    """Copy only one case into a temp dir and execute the agent program there."""
    _harden_hidden_tree()
    with tempfile.TemporaryDirectory(prefix="adsorb_case_") as tmp:
        work = Path(tmp)
        case_dir = work / "case"
        tmp_dir = work / "tmp"
        out_path = work / "result.json"
        prog_path = work / "column_solver.py"
        stdout_path = work / "stdout.txt"
        stderr_path = work / "stderr.txt"
        shutil.copytree(case_src, case_dir)
        shutil.copy2(program_path, prog_path)
        tmp_dir.mkdir()
        for p in (work, case_dir, tmp_dir):
            os.chmod(p, 0o755)
        for f in case_dir.rglob("*"):
            if f.is_file():
                os.chmod(f, 0o644)
        os.chmod(prog_path, 0o755)

        cmd = [
            sys.executable,
            "-s",
            str(prog_path),
            "--case",
            str(case_dir / "case.json"),
            "--out",
            str(out_path),
        ]
        env = _minimal_env(work)
        user_kwargs = _require_unprivileged_agent(env)
        if not user_kwargs:
            # Non-root local runs: allow invoking user's site packages.
            user_site = site.getusersitepackages()
            if user_site:
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
            except OSError as exc:
                raise RuntimeError(
                    f"failed to chown sandbox for unprivileged agent: {exc}"
                ) from exc

        t0 = time.monotonic()
        proc = None
        pgid = None
        stdout_f = open(stdout_path, "w", encoding="utf-8")
        stderr_f = open(stderr_path, "w", encoding="utf-8")
        try:
            proc = subprocess.Popen(
                cmd,
                cwd=str(work),
                env=env,
                stdout=stdout_f,
                stderr=stderr_f,
                start_new_session=True,
                **user_kwargs,
            )
            pgid = proc.pid
            try:
                proc.communicate(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                if pgid is not None:
                    _killpg_quiet(pgid)
                try:
                    proc.wait(timeout=5)
                except Exception:
                    pass
                stdout = stdout_path.read_text(errors="ignore")
                stderr = stderr_path.read_text(errors="ignore")
                pytest.fail(
                    f"case {case_src.name}: program exceeded {timeout_s}s per-run limit "
                    f"(timeout fails the task); stdout={stdout!r} stderr={stderr!r}"
                )
        finally:
            stdout_f.close()
            stderr_f.close()
            # Always reap the process group, including on normal completion.
            if pgid is not None:
                _killpg_quiet(pgid)
            if proc is not None and proc.poll() is None:
                try:
                    proc.kill()
                except Exception:
                    pass
                try:
                    proc.wait(timeout=5)
                except Exception:
                    pass

        elapsed = time.monotonic() - t0
        stdout = stdout_path.read_text(errors="ignore")
        stderr = stderr_path.read_text(errors="ignore")
        if proc is None or proc.returncode != 0:
            rc = None if proc is None else proc.returncode
            pytest.fail(
                f"case {case_src.name}: program exited {rc} in {elapsed:.2f}s\n"
                f"stdout:\n{stdout}\nstderr:\n{stderr}"
            )
        assert out_path.is_file(), f"case {case_src.name}: no output written"
        assert not out_path.is_symlink(), f"case {case_src.name}: symlink output rejected"
        payload = json.loads(out_path.read_text())
        assert isinstance(payload, dict), f"case {case_src.name}: output must be a JSON object"
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
    constant = {
        "qm": 0.43,
        "K": 2.7,
        "k_LDF": 0.027,
        "D_ax": 2.3e-4,
        "L_min_m": 0.12,
    }
    failures = []
    for entry in manifest["hidden_cases"]:
        case = json.loads((HIDDEN_ROOT / entry["case_id"] / "case.json").read_text())
        ok, _ = answer_passes_gates(constant, entry, case)
        if ok:
            failures.append(entry["case_id"])
    assert not failures, f"constant baseline incorrectly passes cases: {failures}"


def test_visible_param_cheats_fail_all_hidden(manifest):
    """Copying any visible sample params (L×1.0/1.1/1.2) must fail every hidden case."""
    assert VISIBLE_TRUTH.is_file(), f"missing visible truth {VISIBLE_TRUTH}"
    visibles = json.loads(VISIBLE_TRUTH.read_text())
    leaks = []
    for entry in manifest["hidden_cases"]:
        case = json.loads((HIDDEN_ROOT / entry["case_id"] / "case.json").read_text())
        for vis in visibles:
            for scale in (1.0, 1.1, 1.2):
                payload = {
                    "qm": vis["parameters"]["qm"],
                    "K": vis["parameters"]["K"],
                    "k_LDF": vis["parameters"]["k_LDF"],
                    "D_ax": vis["parameters"]["D_ax"],
                    "L_min_m": float(vis["L_min_m"]) * scale,
                }
                ok, _ = answer_passes_gates(payload, entry, case)
                if ok:
                    leaks.append(
                        {
                            "hidden": entry["case_id"],
                            "visible": vis["case_id"],
                            "L_scale": scale,
                        }
                    )
    assert not leaks, f"visible-param cheats incorrectly pass: {leaks}"


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
    for key in REQUIRED_KEYS:
        assert key in payload, f"{case_id}: missing key {key}"
        _finite_positive(payload[key], f"{case_id}.{key}")
    assert float(blob["elapsed_s"]) <= float(manifest["per_run_timeout_s"]) + 1e-6
    ok, detail = answer_passes_gates(payload, truth, case)
    assert ok, f"{case_id}: failed verifier gates detail={detail}"
