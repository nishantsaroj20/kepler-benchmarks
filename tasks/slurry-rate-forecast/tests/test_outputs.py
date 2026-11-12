"""Run the forecast program on the held-out campaign and compare it with the latent truth."""

from __future__ import annotations

import json
import math
import os
import shutil
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

PROGRAM = Path("/app/forecast.py")
CAMPAIGN = Path("/tests/hidden_campaign")
TRUTH_PATH = Path("/tests/truth/expected.json")
RUN_LIMIT_SECONDS = 20
REL_RATE = 0.090
ABS_RATE = 0.00015
REL_CONSUMED = 0.060
ABS_CONSUMED = 0.20
FRAC_TOL = 0.030
SUM_TOL = 0.0001
MAX_BYTES = 2_000_000
NUMERIC_KEYS = (
    "rate_mol_m3_s",
    "hydrogen_consumed_mol_m3",
    "resistance_fraction_gas_liquid",
    "resistance_fraction_liquid_solid",
    "resistance_fraction_kinetic",
)
REQUIRED_KEYS = ("case_id",) + NUMERIC_KEYS
FRACTION_KEYS = NUMERIC_KEYS[2:]


def _reject_constant(token: str):
    raise ValueError(token)


def _pairs(pairs):
    obj = {}
    for key, value in pairs:
        if key in obj:
            raise ValueError(f"duplicate key {key}")
        obj[key] = value
    return obj


def load_document(path: Path):
    assert not path.is_symlink()
    assert path.is_file()
    resolved = path.resolve()
    assert resolved.parent == path.parent.resolve()
    raw = path.read_bytes()
    assert 0 < len(raw) <= MAX_BYTES
    return json.loads(raw.decode("utf-8"), object_pairs_hook=_pairs, parse_constant=_reject_constant)


def finite_number(value):
    assert not isinstance(value, bool)
    assert isinstance(value, (int, float))
    number = float(value)
    assert math.isfinite(number)
    return number


def _discard(path: Path):
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.exists():
        shutil.rmtree(path)


def _program_copy() -> Path:
    """Copy the artifact onto a normal file.

    The harness may mount /app/forecast.py read-only and private. Reading
    the bytes here and executing the copy does not depend on changing that
    mount or switching users.
    """
    payload = PROGRAM.read_bytes()
    directory = Path(tempfile.mkdtemp(prefix="slurry-program-"))
    os.chmod(directory, 0o755)
    destination = directory / "forecast.py"
    destination.write_bytes(payload)
    os.chmod(destination, 0o644)
    return destination


def _seal(path: Path):
    """Remove grading data before the submitted program starts."""
    try:
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)
    except OSError:
        return


def _run_program(campaign: Path, cases: Path, output: Path):
    program = _program_copy()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    env["HOME"] = "/tmp"
    env.pop("PYTHONPATH", None)
    command = [
        sys.executable,
        str(program),
        "--campaign",
        str(campaign),
        "--cases",
        str(cases),
        "--output",
        str(output),
    ]

    def _pin():
        try:
            os.sched_setaffinity(0, {0})
        except (AttributeError, OSError):
            return

    proc = None
    grouped = False
    launch_error = None
    for options in (
        {"start_new_session": True, "preexec_fn": _pin},
        {"start_new_session": True},
        {},
    ):
        try:
            proc = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
                **options,
            )
            grouped = bool(options.get("start_new_session"))
            break
        except OSError as exc:
            launch_error = exc
    if proc is None:
        pytest.fail(f"python3 /app/forecast.py could not start: {launch_error}")
    timed_out = False
    try:
        _stdout, stderr = proc.communicate(timeout=RUN_LIMIT_SECONDS)
    except subprocess.TimeoutExpired:
        timed_out = True
        if grouped:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            proc.kill()
        _stdout, stderr = proc.communicate()
    if timed_out:
        _discard(output)
        pytest.fail(f"python3 /app/forecast.py exceeded {RUN_LIMIT_SECONDS} seconds")
    if proc.returncode != 0:
        _discard(output)
        detail = stderr.decode("utf-8", "replace").strip()[-400:]
        pytest.fail(f"python3 /app/forecast.py exited {proc.returncode}: {detail}")


@pytest.fixture(scope="session")
def forecast_document():
    assert PROGRAM.is_file()
    assert not PROGRAM.is_symlink()
    assert CAMPAIGN.is_dir()
    assert TRUTH_PATH.is_file()
    truth = json.loads(TRUTH_PATH.read_text(encoding="utf-8"))
    _seal(TRUTH_PATH)
    work = Path(tempfile.mkdtemp(prefix="slurry-campaign-"))
    shutil.copytree(CAMPAIGN, work, dirs_exist_ok=True)
    out_dir = Path(tempfile.mkdtemp(prefix="slurry-out-"))
    output = out_dir / "predictions.json"
    _run_program(work, work / "predict_cases.csv", output)
    return load_document(output), truth


def test_program_is_a_file():
    assert PROGRAM.is_file()
    assert not PROGRAM.is_symlink()
    text = PROGRAM.read_text(encoding="utf-8")
    assert "--campaign" in text and "--cases" in text and "--output" in text


def test_hidden_forecast(forecast_document):
    forecast_document, truth = forecast_document
    assert set(forecast_document.keys()) == {"cases"}
    rows = forecast_document["cases"]
    assert isinstance(rows, list) and rows
    truth_rows = truth["cases"]
    truth_ids = [row["case_id"] for row in truth_rows]
    assert len(truth_ids) == len(set(truth_ids))
    seen = []
    mapped = {}
    for row in rows:
        assert isinstance(row, dict)
        assert set(row.keys()) == set(REQUIRED_KEYS)
        case_id = row["case_id"]
        assert isinstance(case_id, str) and case_id
        assert case_id not in mapped
        seen.append(case_id)
        mapped[case_id] = row
    assert seen
    assert set(seen) == set(truth_ids)
    assert len(seen) == len(truth_ids)
    missed = []
    for truth_row in truth_rows:
        got = mapped[truth_row["case_id"]]
        numbers = {key: finite_number(got[key]) for key in NUMERIC_KEYS}
        rate = numbers["rate_mol_m3_s"]
        consumed = numbers["hydrogen_consumed_mol_m3"]
        true_rate = float(truth_row["rate_mol_m3_s"])
        true_consumed = float(truth_row["hydrogen_consumed_mol_m3"])
        fractions = [numbers[key] for key in FRACTION_KEYS]
        if abs(sum(fractions) - 1.0) > SUM_TOL:
            missed.append(truth_row["case_id"])
            continue
        if abs(rate - true_rate) > max(REL_RATE * abs(true_rate), ABS_RATE):
            missed.append(truth_row["case_id"])
            continue
        if abs(consumed - true_consumed) > max(REL_CONSUMED * abs(true_consumed), ABS_CONSUMED):
            missed.append(truth_row["case_id"])
            continue
        for key in FRACTION_KEYS:
            if abs(numbers[key] - float(truth_row[key])) > FRAC_TOL:
                missed.append(truth_row["case_id"])
                break
    assert not missed, f"{len(missed)} charges outside the stated limits, including {missed[:8]}"
