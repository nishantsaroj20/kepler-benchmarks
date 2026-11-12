#!/usr/bin/env python3
"""Run the submitted analyzer under an unprivileged, sealed process.

Outputs are written inside a runner-owned work directory (not the pytest temp
tree), then copied back by the root verifier process. This avoids the common
failure where setuid children cannot traverse root-only TemporaryDirectory paths.
"""

from __future__ import annotations

import hashlib
import os
import pwd
import resource
import shutil
import signal
import subprocess
import time
from pathlib import Path

PROGRAM = Path("/app/solver/analyze.py")
RUN_USER = "runner"
PER_RUN_TIMEOUT_SEC = 8
# Virtual address space must be generous: NumPy/OpenBLAS reserve large VAS.
ADDRESS_LIMIT = 8 * 1024 * 1024 * 1024
PROCESS_LIMIT = 64
FILE_LIMIT = 256 * 1024 * 1024
RUN_ID_VAR = "CONTACT_RUN_ID"


def _limits(uid: int, gid: int):
    def apply():
        os.setsid()
        resource.setrlimit(resource.RLIMIT_CPU, (PER_RUN_TIMEOUT_SEC + 2, PER_RUN_TIMEOUT_SEC + 2))
        resource.setrlimit(resource.RLIMIT_AS, (ADDRESS_LIMIT, ADDRESS_LIMIT))
        resource.setrlimit(resource.RLIMIT_NPROC, (PROCESS_LIMIT, PROCESS_LIMIT))
        resource.setrlimit(resource.RLIMIT_FSIZE, (FILE_LIMIT, FILE_LIMIT))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        resource.setrlimit(resource.RLIMIT_NOFILE, (256, 256))
        os.setgroups([])
        os.setgid(gid)
        os.setuid(uid)

    return apply


def _kill_leftovers(pgid: int, run_id: str) -> None:
    try:
        os.killpg(pgid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass
    tag = f"{RUN_ID_VAR}={run_id}".encode()
    for pid in filter(str.isdigit, os.listdir("/proc")):
        try:
            with open(f"/proc/{pid}/environ", "rb") as fh:
                if tag in fh.read():
                    os.kill(int(pid), signal.SIGKILL)
        except (OSError, ValueError):
            continue


def run_analyzer(pack_dir: Path, out_json: Path, program: Path = PROGRAM) -> dict:
    """Execute analyzer as runner; copy JSON result to out_json as root."""
    if not program.is_file():
        raise FileNotFoundError(f"missing program: {program}")
    if out_json.is_symlink():
        raise RuntimeError("output path must not be a symlink")
    if out_json.exists():
        out_json.unlink()

    user = pwd.getpwnam(RUN_USER)
    run_id = hashlib.sha256(f"{pack_dir}:{time.time_ns()}".encode()).hexdigest()
    work = Path("/tmp") / f"run_{run_id[:16]}"
    if work.exists():
        shutil.rmtree(work, ignore_errors=True)
    work.mkdir(mode=0o700)
    os.chown(work, user.pw_uid, user.pw_gid)

    # Runner-owned output and stderr live inside work/
    inner_out = work / "out.json"
    err_path = work / "err.txt"
    for path in (inner_out, err_path):
        fd = os.open(str(path), os.O_CREAT | os.O_WRONLY | os.O_TRUNC | os.O_NOFOLLOW, 0o600)
        os.fchown(fd, user.pw_uid, user.pw_gid)
        os.close(fd)

    os.chmod(program, 0o644)
    # Pack must be readable/traversable by runner without being writable.
    for dirpath, _, filenames in os.walk(pack_dir):
        try:
            os.chmod(dirpath, 0o555)
        except OSError:
            pass
        for name in filenames:
            try:
                os.chmod(os.path.join(dirpath, name), 0o444)
            except OSError:
                pass

    env = {
        "PATH": "/usr/local/bin:/usr/bin:/bin",
        "HOME": str(work),
        "LANG": "C.UTF-8",
        "PYTHONDONTWRITEBYTECODE": "1",
        RUN_ID_VAR: run_id,
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    }
    cmd = [
        "python3",
        str(program),
        "--pack",
        str(pack_dir.resolve()),
        "--out",
        str(inner_out),
    ]
    start = time.time()
    with open(err_path, "wb") as err, open(os.devnull, "rb") as nul:
        proc = subprocess.Popen(
            cmd,
            stdin=nul,
            stdout=subprocess.DEVNULL,
            stderr=err,
            cwd=str(work),
            env=env,
            preexec_fn=_limits(user.pw_uid, user.pw_gid),
            close_fds=True,
        )
        timed_out = False
        try:
            proc.wait(timeout=PER_RUN_TIMEOUT_SEC)
        except subprocess.TimeoutExpired:
            timed_out = True
            proc.kill()
            proc.wait()
    _kill_leftovers(proc.pid, run_id)
    secs = time.time() - start

    err_tail = ""
    try:
        err_tail = err_path.read_text(encoding="utf-8", errors="replace")[-2000:]
    except OSError:
        pass

    # Copy result back for the root verifier (out_json may sit in a root-only temp dir).
    if proc.returncode == 0 and inner_out.is_file() and inner_out.stat().st_size > 0:
        if inner_out.is_symlink():
            shutil.rmtree(work, ignore_errors=True)
            raise RuntimeError("refusing symlink inner output")
        out_json.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(inner_out, out_json)

    shutil.rmtree(work, ignore_errors=True)

    return {
        "status": proc.returncode,
        "timeout": timed_out,
        "seconds": secs,
        "stderr_tail": err_tail,
    }


def load_output_json(path: Path) -> dict:
    if path.is_symlink():
        raise RuntimeError("refusing to read symlink output")
    fd = os.open(str(path), os.O_RDONLY | os.O_NOFOLLOW)
    try:
        data = os.read(fd, 2 * 1024 * 1024)
    finally:
        os.close(fd)
    import json

    obj = json.loads(data.decode("utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("output must be a JSON object")
    return obj
