"""Turn a local clone of carmonalab/scGate_models into normalized corpus records.

Nothing is downloaded and nothing is shipped: the caller points to a clone that sits
at a named tag. Each gating model becomes one label record. The positive genes are the
signature of the model's last level; a gene written ``GENE-`` is a negative gene. The
gating path and the named populations a model excludes are kept as extra fields on the
record, so a reviewer sees them. They do not yet act as a gate.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from .contracts import ContractError

MODEL_SUFFIX = "_scGate_Model.tsv"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _genes(signature: str) -> tuple[list[str], list[str]]:
    positive, negative = [], []
    for item in signature.split(";"):
        item = item.strip()
        if item.endswith("-"):
            negative.append(item[:-1])
        elif item:
            positive.append(item)
    return positive, negative


def _check_tag(models_dir: Path, tag: str) -> None:
    done = subprocess.run(
        ["git", "-C", str(models_dir), "describe", "--tags", "--exact-match"],
        capture_output=True, text=True, check=False,
    )
    if done.returncode != 0 or done.stdout.strip() != tag:
        found = done.stdout.strip() or "no exact tag"
        raise ContractError(
            f"scGate_models clone is at {found}, not {tag}; check out {tag} "
            "or pass the tag you mean to import"
        )


def scgate_payload(
    models_dir: str | Path, collections: list[str], *, species: str, tissue: str, tag: str,
) -> bytes:
    """Build the records JSON for the chosen collections of one species."""
    models_dir = Path(models_dir)
    if not collections:
        raise ContractError("name at least one scGate_models collection")
    if not tissue.strip():
        raise ContractError("name the tissue scope of the imported labels; reviewed labels need one")
    _check_tag(models_dir, tag)
    dictionary = {}
    dictionary_path = models_dir / "CellOntology_dictionary.tsv"
    if dictionary_path.is_file():
        dictionary = {row["scGate_multi"]: row["CellOntology_ID"] for row in _rows(dictionary_path)}
    labels: list[dict[str, Any]] = []
    for collection in collections:
        folder = models_dir / species / collection
        master_path = folder / "master_table.tsv"
        if not master_path.is_file():
            raise ContractError(f"scGate_models collection has no master_table.tsv: {folder}")
        master = {row["name"]: row["signature"] for row in _rows(master_path)}
        for path in sorted(folder.glob(f"*{MODEL_SUFFIX}")):
            name = path.name[: -len(MODEL_SUFFIX)]
            rows = _rows(path)
            levels = sorted({row["levels"] for row in rows}, key=lambda item: int(item.removeprefix("level")))
            gate_path: dict[str, dict[str, list[str]]] = {}
            for level in levels:
                gate_path[level] = {"positive": [], "negative": []}
                for row in rows:
                    if row["levels"] == level:
                        gate_path[level][row["use_as"]].append(row["name"])
            positive: list[str] = []
            negative: list[str] = []
            for row in rows:
                if row["levels"] != levels[-1] or row["use_as"] != "positive":
                    continue
                signature = row["signature"].strip() or master.get(row["name"], "")
                if not signature:
                    raise ContractError(f"{path.name}: no signature for {row['name']}")
                plus, minus = _genes(signature)
                positive += plus
                negative += minus
            excludes = sorted({row["name"] for row in rows if row["use_as"] == "negative"})
            # Decision row 30: each excluded population is tested through its own signature.
            exclusions = []
            for excluded in excludes:
                plus, minus = _genes(master.get(excluded, ""))
                if plus:
                    exclusions.append({"name": excluded, "positive": plus, "negative": minus})
            record_id = f"{species}.{collection}.{name}"
            labels.append({
                "source_record_id": record_id, "name": name,
                "description": (
                    f"scGate model {name} ({collection}, {tag}). Gating path: "
                    + "; ".join(
                        f"{level} needs {', '.join(item['positive']) or 'nothing'}"
                        + (f" and excludes {', '.join(item['negative'])}" if item["negative"] else "")
                        for level, item in gate_path.items()
                    )
                ),
                "positive": list(dict.fromkeys(positive)),
                "negative": [gene for gene in dict.fromkeys(negative) if gene not in positive],
                "source_level": levels[-1],
                "cl_id": dictionary.get(name) or dictionary.get(name.removesuffix("s")),
                "exclusions": exclusions,
                "scgate_gate_path": gate_path, "scgate_excludes": excludes,
                "scgate_excludes_without_signature": sorted(set(excludes) - {item["name"] for item in exclusions}),
                "scgate_file": str(path.relative_to(models_dir)),
                "scgate_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            })
    return json.dumps(
        {"context": {"species": species, "tissue": tissue}, "labels": labels},
        sort_keys=True, indent=2,
    ).encode()
