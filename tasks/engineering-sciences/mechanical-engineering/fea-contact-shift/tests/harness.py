#!/usr/bin/env python3
"""Run the submitted analyzer under an unprivileged, sealed process."""

from __future__ import annotations

import hashlib
import os
import pwd
import resource
import signal
import subprocess
import time
from pathlib import Path

PROGRAM = Path("/app/solver/analyze.py")
RUN_USER = "runner"
PER_RUN_TIMEOUT_SEC = 8
ADDRESS_LIMIT = 2048 * 1024 * 1024
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
    """Execute analyzer as runner. Writes to out_json via O_NOFOLLOW-safe path."""
    if not program.is_file():
        raise FileNotFoundError(f"missing program: {program}")
    if out_json.exists() or out_json.is_symlink():
        # refuse symlink / pre-existing redirect targets
        if out_json.is_symlink():
            raise RuntimeError("output path must not be a symlink")
        out_json.unlink()

    user = pwd.getpwnam(RUN_USER)
    run_id = hashlib.sha256(str(pack_dir).encode()).hexdigest()
    work = Path("/tmp") / f"run_{run_id[:16]}"
    if work.exists():
        # best-effort clean
        import shutil

        shutil.rmtree(work, ignore_errors=True)
    work.mkdir(mode=0o700)
    os.chown(work, user.pw_uid, user.pw_gid)

    # Ensure pack is readable by runner; program is world-readable
    os.chmod(program, 0o644)
    for dirpath, _, filenames in os.walk(pack_dir):
        os.chmod(dirpath, 0o555)
        for name in filenames:
            os.chmod(os.path.join(dirpath, name), 0o444)

    # Create empty output file owned by runner (no symlink)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(str(out_json), os.O_CREAT | os.O_WRONLY | os.O_TRUNC | os.O_NOFOLLOW, 0o600)
    os.fchown(fd, user.pw_uid, user.pw_gid)
    os.close(fd)

    err_path = work / ".err"
    fd_err = os.open(str(err_path), os.O_CREAT | os.O_WRONLY | os.O_TRUNC, 0o600)
    os.fchown(fd_err, user.pw_uid, user.pw_gid)
    os.close(fd_err)

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
        str(pack_dir),
        "--out",
        str(out_json),
    ]
    start = time.time()
    with open(err_path, "wb") as err, open(os.devnull, "rb") as nul:
        proc = subprocess.Popen(
            cmd,
            stdin=nul,
            stdout=subprocess.DEVNULL,  # no capture pipes for agent stdout
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

    # Clean workdir
    import shutil

    shutil.rmtree(work, ignore_errors=True)

    return {
        "status": proc.returncode,
        "timeout": timed_out,
        "seconds": secs,
    }


def load_output_json(path: Path) -> dict:
    if path.is_symlink():
        raise RuntimeError("refusing to read symlink output")
    # O_NOFOLLOW open
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
