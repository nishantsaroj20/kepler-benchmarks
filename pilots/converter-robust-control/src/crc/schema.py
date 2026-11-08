"""Declarative controller.json contract. Research checker, not the future verifier."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from .params import D_SIGMA_MAX, PLANT_SET_ID, SCHEMA_VERSION

MAX_BYTES = 65536
BLOCK_IDS = ("converter_1", "converter_2")
TOP_KEYS = {"schema_version", "plant_set_id", "blocks"}
BLOCK_KEYS = {"id", "A", "B", "C", "D"}


class SchemaError(ValueError):
    pass


def _reject_duplicates(pairs):
    out = {}
    for key, value in pairs:
        if key in out:
            raise SchemaError(f"duplicate key {key}")
        out[key] = value
    return out


def _as_matrix(value, rows: int, cols: int, where: str) -> np.ndarray:
    if isinstance(value, bool) or not isinstance(value, list) or len(value) != rows:
        raise SchemaError(f"{where} must be a {rows}-by-{cols} array")
    out = np.zeros((rows, cols))
    for i, row in enumerate(value):
        if isinstance(row, bool) or not isinstance(row, list) or len(row) != cols:
            raise SchemaError(f"{where} row {i} has the wrong shape")
        for j, entry in enumerate(row):
            if isinstance(entry, bool) or not isinstance(entry, (int, float)):
                raise SchemaError(f"{where}[{i}][{j}] is not a finite number")
            number = float(entry)
            if not np.isfinite(number):
                raise SchemaError(f"{where}[{i}][{j}] is not finite")
            out[i, j] = number
    return out


def load_controller(path: Path) -> list[dict]:
    path = Path(path)
    if path.is_symlink() or os.path.islink(path):
        raise SchemaError("controller path is a symlink")
    if not path.is_file():
        raise SchemaError("controller path is not a regular file")
    size = path.stat().st_size
    if size > MAX_BYTES or size <= 0:
        raise SchemaError(f"file size {size} is outside 1..{MAX_BYTES}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicates)
    except SchemaError:
        raise
    except json.JSONDecodeError as exc:
        raise SchemaError(f"invalid JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise SchemaError("top level must be an object")
    extra = set(payload) - TOP_KEYS
    if extra or set(payload) != TOP_KEYS:
        raise SchemaError(f"top-level keys must be {sorted(TOP_KEYS)}")
    if payload["schema_version"] != SCHEMA_VERSION:
        raise SchemaError("schema_version mismatch")
    if payload["plant_set_id"] != PLANT_SET_ID:
        raise SchemaError("plant_set_id mismatch")
    blocks = payload["blocks"]
    if not isinstance(blocks, list) or len(blocks) != 2:
        raise SchemaError("blocks must be a list of two objects")
    parsed = []
    seen = set()
    for block in blocks:
        if not isinstance(block, dict) or set(block) != BLOCK_KEYS:
            raise SchemaError("each block must contain exactly id, A, B, C, D")
        ident = block["id"]
        if ident not in BLOCK_IDS or ident in seen:
            raise SchemaError(f"illegal or duplicate block id {ident}")
        seen.add(ident)
        parsed.append(
            {
                "id": ident,
                "A": _as_matrix(block["A"], 2, 2, f"{ident}.A"),
                "B": _as_matrix(block["B"], 2, 2, f"{ident}.B"),
                "C": _as_matrix(block["C"], 2, 2, f"{ident}.C"),
                "D": _as_matrix(block["D"], 2, 2, f"{ident}.D"),
            }
        )
    parsed.sort(key=lambda b: b["id"])
    return parsed


def dump_controller(path: Path, blocks: list[dict]) -> None:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "plant_set_id": PLANT_SET_ID,
        "blocks": [
            {
                "id": b["id"],
                "A": np.asarray(b["A"], float).tolist(),
                "B": np.asarray(b["B"], float).tolist(),
                "C": np.asarray(b["C"], float).tolist(),
                "D": np.asarray(b["D"], float).tolist(),
            }
            for b in blocks
        ],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def d_sigma_max(blocks: list[dict]) -> float:
    sig = 0.0
    for block in blocks:
        values = np.linalg.svd(np.asarray(block["D"], float), compute_uv=False)
        sig = max(sig, float(values[0]) if len(values) else 0.0)
    return sig


def gain_ok(blocks: list[dict]) -> bool:
    return d_sigma_max(blocks) <= D_SIGMA_MAX + 1e-9
