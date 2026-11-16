#!/usr/bin/env python3
"""annotation-consensus: a deterministic consensus of standalone annotation tools.

Purpose: combine the per-cell calls of standalone annotation tools (CellTypist, popV,
         SingleR, scANVI) with cell-curator marker-program evidence, and record why each
         cluster and each flagged cell got its label.
Run:     python annotation_consensus.py --demo --output OUT
         python annotation_consensus.py run --input CONFIG.yaml --output OUT [--score] [--force]
         python annotation_consensus.py verify --input CONFIG.yaml --output OUT
         steps: validate votes smooth map-tools consensus markers score report
Inputs (all named in ONE strict YAML config, relative paths resolve against the config):
  cells        h5ad with obs[cluster_key] and obsm[embedding_key] (and X when votes are computed)
  vocabulary   TSV: label, cl_id, positive_markers, negative_markers (corpus vocabulary)
  ontology     cl.obo (referenced by path, never copied)
  tools        one TSV per tool: barcode, predicted_label, predicted_cl_id
  votes        precomputed per-cell marker votes, or computed from the expression matrix
  markers      optional levels/L1/markers.tsv of a cell-curator run
  truth        optional, read ONLY by the score step
Outputs (OUT/): calls.tsv, consensus.tsv, cell_flags.tsv, tool_term_map.tsv, marker_table.tsv,
  report.md, result.json, work/*, reproducibility/{commands.sh,checksums.sha256,
  environment.txt,input_checksums.sha256}, run_manifest.json (not hashed: timestamp, host,
  paths, thread count, input hashes, version drift).

The logic is a port of the measured offline scripts (cov95-2026-11-11 candidates c004, c009,
c010, c012, c014). Every impactful choice is a required key of the config; nothing has a
default. Single process, threads pinned, stable sorts, barcodes aligned by name.
Pinned environment: requirements.lock lists python and package versions; a run refuses any
other version unless --allow-version-drift. environment.txt is hashed, so `verify` passes only
on the same pinned environment. Exact kNN distance ties are broken by barcode rank.
"""

import os
import sys

# Threads are pinned before numpy loads. The config key determinism.threads decides. ------
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def _threads_from_argv(argv):
    for i, arg in enumerate(argv):
        path = None
        if arg in ("--input", "--config") and i + 1 < len(argv):
            path = argv[i + 1]
        elif arg.startswith("--input="):
            path = arg.split("=", 1)[1]
        elif arg == "--demo":
            path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "examples/demo/config.yaml")
        if path:
            try:
                import yaml

                return int(yaml.safe_load(open(path))["determinism"]["threads"])
            except Exception:  # noqa: BLE001 - validate reports a bad config properly
                return 1
    return 1


THREADS = _threads_from_argv(sys.argv[1:])
for _var in THREAD_VARS:
    os.environ[_var] = str(THREADS)

import argparse  # noqa: E402
import copy  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import platform  # noqa: E402
import shutil  # noqa: E402
import socket  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
from datetime import UTC, datetime  # noqa: E402
from importlib import metadata  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from scipy import sparse  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402

SKILL_DIR = Path(__file__).resolve().parent
VERSION = "0.1.0"
LOCK_FILE = SKILL_DIR / "requirements.lock"  # python and package versions a run must match
EXPECTED_FILE = SKILL_DIR / "expected_hashes.json"  # hashed outputs of the demo run
MARKER = ".annotation_consensus_output"  # written into a directory this pipeline created
SCORE_EPS = 9  # decimals: scores are rounded to 1e-9 before an argmax or argsort
EPS = 1e-9  # float tolerance on gate comparisons, as in the measured scripts
UNKNOWN = "Unknown"
CHUNK = 20_000  # as the package
VOTE_CHECK_TOL = 1e-9  # built-in check of the vote loop against the run's scores.tsv
STEPS = ["validate", "votes", "smooth", "map-tools", "consensus", "markers", "score", "report"]
RULES = ["all_agree", "majority", "tie_break_specificity", "tie_break_marker", "tie_break_tool", "abstain"]
RULE_CHOICES = {"majority", "h4", "h5"}
DEMO_CONFIG = SKILL_DIR / "examples/demo/config.yaml"
MARKERS_NOT_BUILT = (
    "Enrichment of marker genes against known markers is NOT built yet; the marker table "
    "lists the genes only."
)


class PipelineError(Exception):
    """A loud, expected failure: the message is the whole story."""


# Strict config schema ---------------------------------------------------------------------
# A leaf is (type, required). Unknown keys fail. An impactful choice is marked so the error
# says that it has no default.
IMPACTFUL = {
    "smoothing.k",
    "smoothing.embedding_key",
    "consensus.rule",
    "consensus.marker_weight",
    "consensus.gate",
    "consensus.tie_break_delta",
    "consensus.snap_to_vocabulary",
    "consensus.marker_scope",
    "consensus.tool_order",
    "consensus.values_status",
    "cluster_gate.share",
    "cluster_gate.lead",
    "cluster_gate.fallback_margin",
    "markers_table.top_n",
    "markers_table.comparator_role",
}
NUM = (int, float)
SCHEMA = {
    "schema_version": (int, True),
    "run": {"name": (str, True), "output_dir": (str, False)},
    "inputs": {
        "cells": (str, True),
        "cluster_key": (str, True),
        "vocabulary": (str, True),
        "ontology": (str, True),
        "tools": (dict, True),
        "markers": (str, False),
    },
    "votes": {
        "mode": (str, True),
        "path": (str, False),
        "compute": {
            "matrix_source": (str, False),
            "matrix_key": (str, False),
            "feature_id_column": (str, False),
            "gene_map": (str, False),
            "level": (str, False),
            "parent_scope": (str, False),
            "check_scores": (str, False),
        },
    },
    "smoothing": {"k": (int, True), "embedding_key": (str, True)},
    "consensus": {
        "rule": (str, True),
        "marker_weight": (NUM, True),
        "gate": (NUM, True),
        "tie_break_delta": (NUM, True),
        "snap_to_vocabulary": (bool, True),
        "marker_scope": (str, True),
        "tool_order": (list, True),
        "values_status": (str, True),
    },
    "cluster_gate": {"share": (NUM, True), "lead": (NUM, True), "fallback_margin": (NUM, True)},
    "markers_table": {"top_n": (int, True), "comparator_role": (str, True)},
    "scoring": {"script": (str, True), "truth": (str, True), "dataset": (str, True)},
    "determinism": {"threads": (int, True), "float_digits": (int, True)},
}
OPTIONAL_BLOCKS = {"votes.compute", "markers_table", "scoring"}


def _check_schema(node, schema, prefix, errors):
    for key in node:
        if key not in schema:
            errors.append(f"unknown key: {prefix}{key}")
    for key, spec in schema.items():
        path = f"{prefix}{key}"
        if key not in node:
            if isinstance(spec, dict):
                if path not in OPTIONAL_BLOCKS:
                    errors.append(f"missing required block: {path}")
            elif spec[1]:
                note = " (impactful choice: no default, set it explicitly)" if path in IMPACTFUL else ""
                errors.append(f"missing required key: {path}{note}")
            continue
        value = node[key]
        if isinstance(spec, dict):
            if not isinstance(value, dict):
                errors.append(f"{path} must be a mapping")
            else:
                _check_schema(value, spec, f"{path}.", errors)
        else:
            types = spec[0] if isinstance(spec[0], tuple) else (spec[0],)
            bad_bool = isinstance(value, bool) and bool not in types
            if bad_bool or not isinstance(value, types):
                errors.append(f"{path} has the wrong type: {type(value).__name__}")


def load_config(path):
    path = Path(path).resolve()
    if not path.is_file():
        raise PipelineError(f"config not found: {path}")
    cfg = yaml.safe_load(path.read_text())
    if not isinstance(cfg, dict):
        raise PipelineError("config must be one YAML mapping")
    errors = []
    _check_schema(cfg, SCHEMA, "", errors)
    if not errors:
        errors += _check_values(cfg)
    if errors:
        raise PipelineError("config is not valid:\n  " + "\n  ".join(errors))
    cfg["_path"] = path
    cfg["_base"] = path.parent
    return cfg


def _check_values(cfg):
    e = []
    c = cfg["consensus"]
    if cfg["schema_version"] != 1:
        e.append("schema_version must be 1")
    if c["rule"] not in RULE_CHOICES:
        e.append(f"consensus.rule must be one of {sorted(RULE_CHOICES)}")
    if c["marker_scope"] not in {"cell", "cluster"}:
        e.append("consensus.marker_scope must be cell or cluster")
    if not 0 <= c["gate"] <= 1:
        e.append("consensus.gate must be in [0, 1]")
    if c["marker_weight"] < 0 or c["tie_break_delta"] < 0:
        e.append("consensus.marker_weight and consensus.tie_break_delta must be >= 0")
    if set(c["tool_order"]) != set(cfg["inputs"]["tools"]) or len(c["tool_order"]) != len(
        cfg["inputs"]["tools"]
    ):
        e.append("consensus.tool_order must list each tool of inputs.tools exactly once")
    if not str(c["values_status"]).strip():
        e.append("consensus.values_status must say whether the values are approved")
    g = cfg["cluster_gate"]
    if not all(0 <= g[k] <= 1 for k in ("share", "lead", "fallback_margin")):
        e.append("cluster_gate values must be in [0, 1]")
    if cfg["smoothing"]["k"] < 1:
        e.append("smoothing.k must be >= 1")
    if cfg["votes"]["mode"] not in {"precomputed", "compute"}:
        e.append("votes.mode must be precomputed or compute")
    if cfg["votes"]["mode"] == "precomputed" and not cfg["votes"].get("path"):
        e.append("votes.path is required when votes.mode is precomputed")
    if cfg["votes"]["mode"] == "compute":
        need = ("matrix_source", "matrix_key", "feature_id_column", "gene_map", "level", "parent_scope")
        comp = cfg["votes"].get("compute", {})
        e += [f"votes.compute.{k} is required when votes.mode is compute" for k in need if k not in comp]
        if comp.get("matrix_source") not in (None, "X", "layer"):
            e.append("votes.compute.matrix_source must be X or layer")
    if cfg["inputs"].get("markers") and "markers_table" not in cfg:
        e.append("markers_table block is required when inputs.markers is set")
    if cfg["determinism"]["threads"] < 1 or not 3 <= cfg["determinism"]["float_digits"] <= 15:
        e.append("determinism.threads >= 1 and 3 <= determinism.float_digits <= 15")
    return e


def resolve(cfg, value):
    p = Path(value)
    return p if p.is_absolute() else (cfg["_base"] / p).resolve()


# Small helpers ----------------------------------------------------------------------------
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fmt(x, digits):
    """A float as a fixed number of significant digits; NaN as NA."""
    return "NA" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{float(x):.{digits}g}"


def fmt_col(values, digits):
    return np.array([fmt(v, digits) for v in values], dtype=object)


def jdump(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def write_tsv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def read_table(path, **kw):
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str, **kw)


def natural_key(value):
    return (0, int(value), "") if str(value).lstrip("-").isdigit() else (1, 0, str(value))


def rkey(values):
    """Scores rounded to 1e-9: the key of every argmax and argsort, so a near-tie is a tie."""
    return np.round(np.asarray(values, dtype=np.float64), SCORE_EPS)


# Pinned environment -----------------------------------------------------------------------
def read_lock(path=None):
    """requirements.lock -> {"python": version, package: version}."""
    path = Path(LOCK_FILE if path is None else path)
    if not path.is_file():
        raise PipelineError(f"lock file not found: {path}")
    pinned = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("# python "):
            pinned["python"] = line.split()[-1]
        elif line and not line.startswith("#"):
            name, _, version = line.partition("==")
            pinned[name.strip()] = version.strip()
    if "python" not in pinned:
        raise PipelineError(f"lock file {path} has no '# python X.Y.Z' line")
    return pinned


def dist_version(name):
    """Installed version of a package. A source checkout on sys.path (pytest adds src/) can carry a
    stale *.egg-info of another version, so an installed *.dist-info wins over it."""
    found = list(metadata.distributions(name=name))
    if not found:
        return "not installed"
    installed = [d for d in found if str(getattr(d, "_path", "")).endswith(".dist-info")]
    return (installed or found)[0].version


def installed_versions(names):
    out = {"python": platform.python_version()}
    out.update({name: dist_version(name) for name in names if name != "python"})
    return out


def version_differences(lock=None):
    """Lines 'name: lock X, installed Y' for every pinned version that differs."""
    pinned = read_lock() if lock is None else lock
    have = installed_versions(pinned)
    return [f"{k}: lock {pinned[k]}, installed {have[k]}" for k in sorted(pinned) if have[k] != pinned[k]]


LOCK_HASH_KEY = "requirements.lock-sha256"  # the line of environment.txt that holds the lock hash


def check_versions(allow):
    """Refuse a run on a version that differs from the lock file, unless allowed (and recorded)."""
    diffs = version_differences()
    if diffs and not allow:
        raise PipelineError(
            "the environment differs from requirements.lock (use the pinned environment, or "
            "--allow-version-drift, which is recorded in run_manifest.json and result.json):\n  "
            + "\n  ".join(diffs)
        )
    return diffs


# Inputs and the output directory ----------------------------------------------------------
PATH_KEYS = (
    "inputs.cells", "inputs.vocabulary", "inputs.ontology", "inputs.markers", "votes.path",
    "votes.compute.gene_map", "votes.compute.check_scores", "scoring.script", "scoring.truth",
)  # fmt: skip


def canonical_config(cfg):
    """The config as hashed: no thread count, no output dir, file names instead of paths."""
    c = copy.deepcopy({k: v for k, v in cfg.items() if not k.startswith("_")})
    c["determinism"].pop("threads")
    c["run"].pop("output_dir", None)
    c["inputs"]["tools"] = {t: Path(p).name for t, p in c["inputs"]["tools"].items()}
    for dotted in PATH_KEYS:
        *parents, leaf = dotted.split(".")
        node = c
        for part in parents:
            node = node.get(part, {}) if isinstance(node, dict) else {}
        if isinstance(node, dict) and node.get(leaf):
            node[leaf] = Path(node[leaf]).name
    return c


def config_sha256(cfg):
    return hashlib.sha256(jdump(canonical_config(cfg)).encode()).hexdigest()


def input_hashes(cfg):
    """key -> (file name, sha256) of every input; the config is hashed in canonical form."""
    out = {}
    for key, p in input_files(cfg).items():
        if key == "config":  # the config file's own name and location are arbitrary: hash its canonical form
            out[key] = ("canonical-config", config_sha256(cfg))
        else:
            out[key] = (p.name, sha256_file(p))
    return out


def unsafe_output(out):
    """The cell-curator run directory that holds `out`, or None."""
    out = Path(os.path.abspath(out))
    for p in (out, *out.parents):
        if (p / "run_state.json").exists() or (p / "annotation" / "levels").is_dir():
            return p
    return None


class Ctx:
    """Config, output directory and the overwrite rule."""

    def __init__(self, cfg, out, force, argv=None, drift=(), allow_drift=False):
        self.cfg, self.out, self.force = cfg, Path(out), force
        self.argv, self.drift, self.allow_drift = list(argv or []), list(drift), allow_drift
        self.digits = cfg["determinism"]["float_digits"]

    def owned(self):
        """True when the pipeline created the output directory (it wrote the marker)."""
        return (self.out / MARKER).is_file()

    def mark(self):
        """Write the marker into a directory that is new or empty; leave other directories alone."""
        self.out.mkdir(parents=True, exist_ok=True)
        if not self.owned() and not any(self.out.iterdir()):
            (self.out / MARKER).write_text("created by annotation_consensus.py; --force may remove its outputs\n")

    def claim(self, *relative):
        """Refuse to overwrite an existing output unless --force."""
        paths = [self.out / r for r in relative]
        clash = [str(p.relative_to(self.out)) for p in paths if p.exists()]
        if clash and not self.force:
            raise PipelineError(f"refusing to overwrite {clash}; use --force or a new --output")
        for p in paths:
            p.parent.mkdir(parents=True, exist_ok=True)
        return paths


# Loading ----------------------------------------------------------------------------------
def load_ontology(cfg):
    from cell_curator.ontology import CellOntology

    return CellOntology.load(resolve(cfg, cfg["inputs"]["ontology"]))


def load_cells(cfg, expression=False):
    """Barcodes in the h5ad's order, cluster ids, embedding."""
    path = resolve(cfg, cfg["inputs"]["cells"])
    adata = ad.read_h5ad(path, backed=None if expression else "r")
    key = cfg["inputs"]["cluster_key"]
    if key not in adata.obs:
        raise PipelineError(f"cells h5ad has no obs column {key!r}")
    names = pd.Index(adata.obs_names.astype(str))
    if names.has_duplicates:
        raise PipelineError("cells h5ad has duplicate barcodes")
    return adata, names, adata.obs[key].astype(str).to_numpy()


def load_vocabulary(cfg):
    vocab = pd.read_csv(
        resolve(cfg, cfg["inputs"]["vocabulary"]), sep="\t", keep_default_na=False, dtype=str
    )
    for col in ("label", "cl_id"):
        if col not in vocab:
            raise PipelineError(f"vocabulary has no column {col!r}")
    small = vocab[["label", "cl_id"]].drop_duplicates().sort_values("label", kind="stable")
    if small["label"].duplicated().any():
        raise PipelineError("a vocabulary label has more than one cl_id")
    return vocab, small.reset_index(drop=True)


def load_votes(ctx, names):
    cfg = ctx.cfg
    if cfg["votes"]["mode"] == "precomputed":
        path = resolve(cfg, cfg["votes"]["path"])
    else:
        path = ctx.out / "work/cell_votes.tsv"
    votes = read_table(path)
    need = {"barcode", "top1", "voted"}
    if not need <= set(votes.columns):
        raise PipelineError(f"votes file lacks columns {sorted(need - set(votes.columns))}")
    votes["voted"] = votes["voted"].astype(str).str.lower().isin({"true", "1"})
    votes = votes.set_index("barcode")
    missing = names.difference(votes.index)
    if len(missing):
        raise PipelineError(f"votes miss {len(missing)} barcodes, first: {list(missing[:3])}")
    return votes.loc[names]  # aligned by name


def load_tool(path, names, tool):
    pred = pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str).set_index("barcode")
    missing = names.difference(pred.index)
    if len(missing):
        raise PipelineError(f"tool {tool} misses {len(missing)} barcodes, first: {list(missing[:3])}")
    pred = pred.loc[names]
    pred.loc[pred["predicted_label"].eq(UNKNOWN), "predicted_cl_id"] = ""  # no id = abstention
    return pred


def input_files(cfg):
    """name -> absolute path of every input the run reads, truth excluded."""
    files = {"config": cfg["_path"]}
    inp = cfg["inputs"]
    files["cells"] = resolve(cfg, inp["cells"])
    files["vocabulary"] = resolve(cfg, inp["vocabulary"])
    files["ontology"] = resolve(cfg, inp["ontology"])
    for tool, p in inp["tools"].items():
        files[f"tool:{tool}"] = resolve(cfg, p)
    if cfg["votes"]["mode"] == "precomputed":
        files["votes"] = resolve(cfg, cfg["votes"]["path"])
    else:
        comp = cfg["votes"]["compute"]
        if comp.get("gene_map"):
            files["gene_map"] = resolve(cfg, comp["gene_map"])
        if comp.get("check_scores"):
            files["votes_check_scores"] = resolve(cfg, comp["check_scores"])
    if inp.get("markers"):
        files["markers"] = resolve(cfg, inp["markers"])
    return files


# Step: validate ---------------------------------------------------------------------------
def step_validate(ctx):
    cfg = ctx.cfg
    errors, info = [], {}
    files = input_files(cfg)
    for name, p in files.items():
        if not p.is_file():
            errors.append(f"input file not found: {name} -> {p}")
    if errors:
        raise PipelineError("validate failed:\n  " + "\n  ".join(errors))
    ont = load_ontology(cfg)
    adata, names, clusters = load_cells(cfg)
    emb = cfg["smoothing"]["embedding_key"]
    if emb not in adata.obsm:
        errors.append(f"cells h5ad has no obsm[{emb!r}]")
    if cfg["smoothing"]["k"] + 1 > len(names):
        errors.append("smoothing.k + 1 exceeds the number of cells")
    vocab, small = load_vocabulary(cfg)
    bad = [x for x in small["cl_id"] if x and not _resolves(ont, x)]
    if bad:
        errors.append(f"vocabulary cl_id not in the ontology: {sorted(set(bad))}")
    ids = [ont.resolve(x).accession for x in small["cl_id"] if x and _resolves(ont, x)]
    if len(ids) != len(set(ids)):
        errors.append("two vocabulary labels resolve to the same ontology term")
    if cfg["votes"]["mode"] == "precomputed":
        votes = read_table(files["votes"])
        only_votes, only_cells = set(votes["barcode"]) - set(names), set(names) - set(votes["barcode"])
        if only_cells:
            errors.append(f"{len(only_cells)} cell barcodes have no vote row (matched by name)")
        if only_votes:
            info["votes_extra_barcodes_ignored"] = len(only_votes)
        outside = sorted(set(votes["top1"]) - set(small["label"]))
        if outside:
            errors.append(f"vote labels outside the vocabulary: {outside[:5]}")
    info["tools"] = {}
    for tool in cfg["consensus"]["tool_order"]:
        pred = pd.read_csv(files[f"tool:{tool}"], sep="\t", keep_default_na=False, dtype=str)
        if not {"barcode", "predicted_label", "predicted_cl_id"} <= set(pred.columns):
            errors.append(f"tool {tool}: needs columns barcode, predicted_label, predicted_cl_id")
            continue
        gap = set(names) - set(pred["barcode"])
        if gap:
            errors.append(f"tool {tool}: {len(gap)} cell barcodes missing (matched by name)")
        called = pred[pred["predicted_label"].ne(UNKNOWN) & pred["predicted_cl_id"].ne("")]
        unmapped = sorted({x for x in called["predicted_cl_id"].unique() if not _resolves(ont, x)})
        if unmapped:
            errors.append(f"tool {tool}: labels with no ontology term: {unmapped[:5]}")
        info["tools"][tool] = {
            "n_rows": int(len(pred)),
            "n_extra_barcodes_ignored": int(len(set(pred["barcode"]) - set(names))),
            "n_abstain": int((pred["predicted_label"].eq(UNKNOWN) | pred["predicted_cl_id"].eq("")).sum()),
        }
    if cfg["inputs"].get("markers"):
        head = pd.read_csv(files["markers"], sep="\t", nrows=1, dtype=str)
        need = {"cluster_id", "method", "feature", "rank", "comparator_role"}
        if not need <= set(head.columns):
            errors.append(f"markers table lacks columns {sorted(need - set(head.columns))}")
    if errors:
        raise PipelineError("validate failed:\n  " + "\n  ".join(errors))
    (path,) = ctx.claim("validation.json")
    info.update(
        {
            "ok": True,
            "n_cells": int(len(names)),
            "n_clusters": int(len(set(clusters))),
            "n_vocabulary_labels": int(len(small)),
            "n_vocabulary_labels_without_cl_id": int((small["cl_id"] == "").sum()),
            "votes_mode": cfg["votes"]["mode"],
        }
    )
    path.write_text(json.dumps(info, indent=1, sort_keys=True) + "\n")
    print(f"validate: ok, {len(names)} cells, {len(set(clusters))} clusters, {len(small)} labels")


def _resolves(ont, value):
    try:
        ont.resolve(value)
    except Exception:  # noqa: BLE001 - the ontology raises its own ContractError
        return False
    return True


# Step: votes (port of c004/c009 01_cell_votes.py and evidence.score_marker_programs) -------
def step_votes(ctx):
    cfg = ctx.cfg
    (status_path,) = ctx.claim("work/votes_status.json")
    if cfg["votes"]["mode"] == "precomputed":
        p = resolve(cfg, cfg["votes"]["path"])
        status = {"mode": "precomputed", "file": p.name, "sha256": sha256_file(p)}
        status_path.write_text(json.dumps(status, indent=1, sort_keys=True) + "\n")
        print(f"votes: skipped, per-cell votes supplied precomputed ({p.name})")
        return
    comp = cfg["votes"]["compute"]
    (votes_path,) = ctx.claim("work/cell_votes.tsv")
    vocab, _ = load_vocabulary(cfg)
    for col in ("axis", "level", "parent_scope", "positive_markers", "negative_markers"):
        if col not in vocab:
            raise PipelineError(f"vocabulary has no column {col!r}, needed to compute votes")
    sel = vocab["axis"].eq("cell_type") & vocab["level"].eq(comp["level"])
    labels = set(vocab.loc[sel, "label"])
    programs = vocab[sel & vocab["parent_scope"].eq(comp["parent_scope"])]
    adata, names, clusters = load_cells(cfg, expression=True)
    matrix = adata.X if comp["matrix_source"] == "X" else adata.layers[comp["matrix_key"]]
    if comp["gene_map"]:
        gm = pd.read_csv(resolve(cfg, comp["gene_map"]), sep="\t", dtype=str)
        lookup = dict(zip(gm["feature_id"], gm["gene_symbol"], strict=True))
        features = pd.Index([lookup[v] for v in adata.var_names.astype(str)])
    elif comp["feature_id_column"]:
        features = pd.Index(adata.var[comp["feature_id_column"]].astype(str))
    else:
        features = pd.Index(adata.var_names.astype(str))
    if features.has_duplicates:
        raise PipelineError("feature ids are not unique")
    from cell_curator.contracts import feature_positions

    positions = feature_positions(features)
    scored = []
    for _, row in programs.iterrows():
        if row["label"] not in labels:
            continue
        pos = [positions[g] for g in row["positive_markers"].split("|") if g and g in positions]
        neg = [positions[g] for g in row["negative_markers"].split("|") if g and g in positions]
        if pos:
            scored.append((row["label"], pos, neg))
    scored.sort(key=lambda x: x[0])  # label name order: a tie goes to the first label name
    cols = sorted({i for _, a, b in scored for i in (*a, *b)})
    local = {f: k for k, f in enumerate(cols)}

    def block(start, stop):
        rows = matrix[start:stop][:, cols]
        return rows.toarray() if sparse.issparse(rows) else np.asarray(rows, dtype=float)

    n = len(names)
    total, total_sq = np.zeros(len(cols)), np.zeros(len(cols))
    for s in range(0, n, CHUNK):
        v = block(s, min(s + CHUNK, n)).astype(float)
        total += v.sum(0)
        total_sq += (v**2).sum(0)
    mean = total / n
    spread = np.sqrt(np.maximum(total_sq / n - mean**2, 0.0))
    spread = np.where(spread <= 1e-12, 1.0, spread)
    blocks = []
    for s in range(0, n, CHUNK):
        z = (block(s, min(s + CHUNK, n)) - mean) / spread
        blocks.append(
            np.column_stack(
                [
                    z[:, [local[i] for i in pos]].mean(1)
                    - (np.maximum(z[:, [local[i] for i in neg]], 0.0).mean(1) if neg else 0.0)
                    for _, pos, neg in scored
                ]
            )
        )
    scores = np.vstack(blocks)
    prog_names = np.asarray([p for p, _, _ in scored])
    order = np.argsort(-rkey(scores), axis=1, kind="stable")
    best, second, rows = order[:, 0], order[:, 1], np.arange(n)
    cells = pd.DataFrame(
        {
            "barcode": names,
            "cluster_id": clusters,
            "top1": prog_names[best],
            "top2": prog_names[second],
            "score1": fmt_col(scores[rows, best], ctx.digits),
            "score2": fmt_col(scores[rows, second], ctx.digits),
            "voted": scores[rows, best] > 0,  # a cell with max score <= 0 does not vote
        }
    )
    status = {"mode": "compute", "n_programs": len(scored)}
    if comp.get("check_scores"):
        saved = pd.read_csv(resolve(cfg, comp["check_scores"]), sep="\t", dtype={"cluster_id": str})
        saved = saved.set_index(["cluster_id", "label"])["raw_vote_fraction"]
        mine = (
            cells[cells["voted"]].groupby(["cluster_id", "top1"]).size()
            / cells.groupby("cluster_id").size().rename("n")
        ).rename("mine")
        mine.index.names = ["cluster_id", "label"]
        joined = pd.concat([saved, mine], axis=1).fillna({"mine": 0.0})
        diff = float((joined["raw_vote_fraction"] - joined["mine"]).abs().max())
        status["check_max_abs_diff"] = fmt(diff, 6)
        if not diff < VOTE_CHECK_TOL:
            raise PipelineError(
                f"vote check failed: max abs diff {diff:.3g} vs raw_vote_fraction (limit {VOTE_CHECK_TOL})"
            )
        print(f"votes: check against scores.tsv passed, max abs diff {diff:.3g}")
    write_tsv(cells, votes_path)
    status_path.write_text(json.dumps(status, indent=1, sort_keys=True) + "\n")
    print(f"votes: computed for {n} cells with {len(scored)} programs")


# Step: smooth (port of c009/c010/c012 01_evidence.py) --------------------------------------
def exact_sq_dist(emb, idx):
    """Squared distance from row i of emb to each row idx[i], summed over features in float64.

    The sum depends only on the two points, never on where either row sits, so the same pair
    gives the same bits in any row order and at any thread count.
    """
    out = np.empty(idx.shape)
    step = max(1, 2_000_000 // max(1, idx.shape[1] * emb.shape[1]))
    for s in range(0, len(idx), step):
        diff = emb[s : s + step, None, :] - emb[idx[s : s + step]]
        out[s : s + step] = (diff * diff).sum(axis=-1)
    return out


def knn_indices(emb, m):
    """The m nearest rows of emb to each row, ordered by exact distance, then by row index.

    Callers pass rows sorted by barcode, so the row index is the barcode rank and a tie does not
    depend on the input row order. emb is cast to float64. sklearn finds candidates on one thread
    (chunking cannot depend on determinism.threads) but its distances carry rounding noise and
    break exact ties in an undefined order. So ask for a window of candidates, recompute their
    distances with exact_sq_dist, sort by (distance, index) and keep m. The window widens until
    its farthest candidate is more than a rounding slack (1e-10 of the largest squared norm)
    beyond the m-th, or it covers every row. Two points tie only if their float64 squared
    distances are equal bit for bit.
    """
    emb = np.ascontiguousarray(emb, dtype=np.float64)
    n = len(emb)
    window = min(n, m + 8)
    slack = 1e-10 * max(1.0, float((emb * emb).sum(axis=1).max()))
    nn = NearestNeighbors(n_neighbors=window, algorithm="brute", n_jobs=1).fit(emb)
    while True:
        idx = nn.kneighbors(emb, n_neighbors=window, return_distance=False)
        d2 = exact_sq_dist(emb, idx)
        order = np.lexsort((idx, d2), axis=1)  # the last key sorts first: distance, then index
        idx, d2 = np.take_along_axis(idx, order, axis=1), np.take_along_axis(d2, order, axis=1)
        if window >= n or not (d2[:, -1] <= d2[:, m - 1] + slack).any():
            return idx[:, :m]
        window = min(n, window * 2)


def step_smooth(ctx):
    cfg = ctx.cfg
    paths = ctx.claim("work/smoothed_votes.npy", "work/smoothed_labels.tsv", "work/cells.tsv")
    adata, names, clusters = load_cells(cfg)
    votes = load_votes(ctx, names)
    _, small = load_vocabulary(cfg)
    labels = small["label"].to_numpy()  # sorted, every vocabulary label
    col = np.searchsorted(labels, votes["top1"].to_numpy())
    col = np.minimum(col, len(labels) - 1)
    if not (labels[col] == votes["top1"].to_numpy()).all():
        raise PipelineError("a top1 vote label is outside the vocabulary")
    n, k = len(names), cfg["smoothing"]["k"]
    onehot = sparse.csr_matrix(
        (votes["voted"].to_numpy(float), (np.arange(n), col)), shape=(n, len(labels))
    )
    emb = np.asarray(adata.obsm[cfg["smoothing"]["embedding_key"]])
    # Search on rows sorted by barcode: a distance tie goes to the lower barcode, so the neighbour
    # set does not depend on the order of the input rows. Tables stay in the input order.
    perm = np.argsort(names.to_numpy(), kind="stable")
    idx = perm[knn_indices(emb[perm], k + 1)]  # includes the cell itself; original row numbers
    rows = np.repeat(perm, k + 1)
    graph = sparse.csr_matrix((np.ones(rows.size) / (k + 1), (rows, idx.ravel())), shape=(n, n))
    smooth = (graph @ onehot).toarray()  # float64, as c009/c010
    np.save(paths[0], smooth)
    write_tsv(small, paths[1])
    write_tsv(pd.DataFrame({"barcode": names, "cluster_id": clusters}), paths[2])
    print(f"smooth: k={k} over {cfg['smoothing']['embedding_key']}, {n} cells x {len(labels)} labels")


# Step: map-tools (port of c012 snap map) ---------------------------------------------------
def up_dist(ont, term):
    """Shortest is_a steps from term up to each of its ancestors."""
    dist, frontier, d = {}, [term], 0
    while frontier:
        d += 1
        nxt = []
        for x in frontier:
            for p in ont.terms[x].parents:
                if p in ont.terms and p not in dist:
                    dist[p] = d
                    nxt.append(p)
        frontier = nxt
    return dist


def vocab_ids(ont, small):
    keep = small["cl_id"].ne("").to_numpy()
    return keep, [ont.resolve(x).accession for x in small["cl_id"][keep]]


def step_map_tools(ctx):
    cfg = ctx.cfg
    paths = ctx.claim("tool_term_map.tsv", "work/tools_mapped.tsv")
    ont = load_ontology(cfg)
    _, names, _ = load_cells(cfg)
    _, small = load_vocabulary(cfg)
    _, vids = vocab_ids(ont, small)
    vset = set(vids)
    snap = cfg["consensus"]["snap_to_vocabulary"]
    raw = {}
    for tool in cfg["consensus"]["tool_order"]:
        pred = load_tool(resolve(cfg, cfg["inputs"]["tools"][tool]), names, tool)
        raw[tool] = pred["predicted_cl_id"].map(lambda x: ont.resolve(x).accession if x else "")
    universe = sorted({t for s in raw.values() for t in s if t} | vset)
    up = {t: up_dist(ont, t) for t in universe}
    mapped_to, rows = {}, []
    for t in universe:
        status, target, steps = "in_vocabulary", t, ""
        if t not in vset:
            best = None
            for v in sorted(vset):  # nearest vocabulary term on the lineage, first sorted wins ties
                dist = up[v].get(t, up[t].get(v))
                if dist is not None and (best is None or dist < best[0]):
                    best = (dist, v)
            if not snap:
                status, target = "not_snapped_option_off", t
            elif best:
                status, target, steps = "snapped", best[1], str(best[0])
            else:
                status, target = "no_vocabulary_term_on_lineage", t
        mapped_to[t] = target
        rows.append((t, status, target, steps))
    table = []
    for tool, series in raw.items():
        counts = series[series.ne("")].value_counts()
        for t, nn in sorted(counts.items()):
            _, status, target, steps = next(r for r in rows if r[0] == t)
            table.append(
                {
                    "tool": tool,
                    "term": t,
                    "name": ont.resolve(t).name,
                    "status": status,
                    "mapped_to": target,
                    "mapped_name": ont.resolve(target).name,
                    "is_a_steps": steps,
                    "n_cells": int(nn),
                }
            )
        table.append(
            {
                "tool": tool, "term": "", "name": UNKNOWN, "status": "abstain", "mapped_to": "",
                "mapped_name": "", "is_a_steps": "", "n_cells": int((series == "").sum()),
            }
        )
    write_tsv(pd.DataFrame(table), paths[0])
    wide = {"barcode": names}
    for tool, series in raw.items():
        wide[f"{tool}__raw"] = series.to_numpy()
        wide[f"{tool}__term"] = series.map(lambda t: mapped_to[t] if t else "").to_numpy()
    write_tsv(pd.DataFrame(wide), paths[1])
    n_snapped = sum(r["n_cells"] for r in table if r["status"] == "snapped")
    print(f"map-tools: {len(raw)} tools, snap={snap}, {n_snapped} tool calls snapped")


# Step: consensus (port of c012/c014 voting; c009 marker-only cluster rule) ------------------
def order_terms(ont, ids):
    """Term ids in label-name order (then id): the lowest index is the first label name."""
    return sorted(set(ids), key=lambda t: (ont.resolve(t).name, t))


class Voting:
    """Term universe, ancestry matrices and the vote functions of c014.

    `terms` must come from order_terms: every argmax and argsort below takes the lowest term
    index on a tie, which is then the first label name. Scores are rounded to 1e-9 first.
    """

    def __init__(self, ont, terms, w_mark):
        self.ont, self.terms, self.w = ont, terms, w_mark
        self.tix = {t: i for i, t in enumerate(terms)}
        self.nT = len(terms)
        self.DESC = np.zeros((self.nT, self.nT), bool)  # DESC[a, b]: b is a proper descendant of a
        for b, t in enumerate(terms):
            for a in ont.ancestors(t):
                if a in self.tix:
                    self.DESC[self.tix[a], b] = True
        self.names = np.array([ont.resolve(t).name for t in terms])

    def tool_counts(self, C):
        n, M = C.shape
        cnt = np.zeros((n, self.nT), np.float32)
        bonus = np.zeros((n, self.nT), np.float32)  # tiny rank bonus: the better tool wins exact ties
        for m in range(M):
            ok = np.where(C[:, m] >= 0)[0]
            cnt[ok, C[ok, m]] += 1
            bonus[ok, C[ok, m]] = np.maximum(bonus[ok, C[ok, m]], (M - m) * 1e-3)
        return cnt, bonus

    def plain(self, C, S=None):
        """majority (S None) and h4 (S given): weighted vote, argmax, no tie logic."""
        cnt, bonus = self.tool_counts(C)
        w = self.w if S is not None else 0.0
        tot = cnt + (w * S if S is not None else 0)
        win = rkey(tot + bonus).argmax(1)
        share = tot[np.arange(len(C)), win] / (C.shape[1] + w)
        return np.where(tot.max(1) > 0, win, -1), share, np.zeros(len(C), int)

    def h5(self, C, S, delta):
        """h4 vote; where top total - second total <= delta, break the tie by
        (1) specificity, (2) marker share, (3) rank of the backing tool."""
        n, M = C.shape
        cnt, bonus = self.tool_counts(C)
        tot = cnt + self.w * S
        top = tot.max(1)
        mask = (tot >= top[:, None] - delta - EPS) & (tot > EPS)
        ntied = mask.sum(1)
        win = np.where(top > EPS, rkey(tot + bonus).argmax(1), -1)
        rule = np.full(n, "majority", object)
        for i in np.where(ntied >= 2)[0]:
            tied = np.where(mask[i])[0]
            sub = self.DESC[np.ix_(tied, tied)]  # sub[a, b]: tied[b] below tied[a]
            below = tied[sub.any(0)]
            if len(below):  # step 1: a tied term that descends from another tied term
                win[i] = below[rkey(tot[i, below]).argmax()]
                rule[i] = "tie_break_specificity"
                continue
            s = rkey(S[i, tied])
            if s.max() > EPS and (np.sort(s)[-1] - np.sort(s)[-2]) > EPS:  # step 2: unique top marker share
                win[i] = tied[s.argmax()]
                rule[i] = "tie_break_marker"
                continue
            rank = [min([m for m in range(M) if C[i, m] == t], default=99) for t in tied]
            if min(rank) < 99:  # step 3: best tool
                win[i] = tied[int(np.argmin(rank))]
                rule[i] = "tie_break_tool"
            else:
                win[i] = -1
                rule[i] = "abstain"
        share = np.where(win >= 0, tot[np.arange(n), np.where(win >= 0, win, 0)] / (M + self.w), 0.0)
        return win, share, ntied, rule


def cluster_mean(S, cl):
    out = np.zeros((cl.max() + 1, S.shape[1]), S.dtype)
    np.add.at(out, cl, S)
    out /= np.bincount(cl, minlength=cl.max() + 1)[:, None]
    return out


def cluster_plurality(term, cl, nC, nT):
    """Per cluster: most common term (first label name on ties), its count; -1 / 0 if none."""
    called = term >= 0
    cnt = np.bincount(cl[called] * nT + term[called], minlength=nC * nT).reshape(nC, nT)
    return np.where(cnt.sum(1) > 0, cnt.argmax(1), -1), cnt


def step_consensus(ctx):
    cfg = ctx.cfg
    dg = ctx.digits
    c, cg = cfg["consensus"], cfg["cluster_gate"]
    paths = ctx.claim("calls.tsv", "consensus.tsv", "cell_flags.tsv", "work/consensus_summary.json")
    ont = load_ontology(cfg)
    _, names, _ = load_cells(cfg)
    cells = pd.read_csv(ctx.out / "work/cells.tsv", sep="\t", dtype=str)
    if not cells["barcode"].equals(pd.Series(names, name="barcode")):
        raise PipelineError("work/cells.tsv does not match the h5ad barcodes; rerun smooth")
    _, small = load_vocabulary(cfg)
    keep, vids = vocab_ids(ont, small)
    smooth = np.load(ctx.out / "work/smoothed_votes.npy")
    tools = c["tool_order"]
    mapped = pd.read_csv(ctx.out / "work/tools_mapped.tsv", sep="\t", keep_default_na=False, dtype=str)
    if not mapped["barcode"].equals(cells["barcode"]):
        raise PipelineError("work/tools_mapped.tsv does not match the cells; rerun map-tools")
    terms = order_terms(ont, {t for tool in tools for t in mapped[f"{tool}__term"] if t} | set(vids))
    V = Voting(ont, terms, float(c["marker_weight"]))
    tix, nT = V.tix, V.nT
    n = len(cells)
    cl = pd.factorize(cells["cluster_id"])[0]
    nC = cl.max() + 1
    S = np.zeros((n, nT), np.float32)
    S[:, [tix[t] for t in vids]] = smooth[:, keep].astype(np.float32)
    C = np.stack(
        [np.array([tix[t] if t else -1 for t in mapped[f"{tool}__term"]]) for tool in tools], 1
    )

    # Votes ------------------------------------------------------------------------------
    S_used = cluster_mean(S, cl)[cl] if c["marker_scope"] == "cluster" else S
    if c["rule"] == "majority":
        win, share, ntied = V.plain(C)
        rule = np.full(n, "majority", object)
    elif c["rule"] == "h4":
        win, share, ntied = V.plain(C, S_used)
        rule = np.full(n, "majority", object)
    else:
        win, share, ntied, rule = V.h5(C, S_used, float(c["tie_break_delta"]))
    all_same = (C >= 0).all(1) & (C == C[:, [0]]).all(1)
    rule[all_same & (win == C[:, 0])] = "all_agree"
    rule[win < 0] = "abstain"
    called = (share >= c["gate"] - EPS) & (win >= 0)
    rule_final = np.where(called, rule, "abstain")
    strict = np.where(called, win, -1)

    # Marker-only call (c009 cluster rule, c010 frozen): share, lead, fallback margin -----
    full = smooth  # float64, every vocabulary label, as c009
    flabels = small["label"].to_numpy()
    fcl = dict(zip(small["label"], small["cl_id"], strict=True))
    order = np.argsort(-rkey(full), axis=1, kind="stable")  # flabels are sorted: a tie goes to the first name
    cell_best = full[np.arange(n), order[:, 0]]
    cell_margin = np.clip(np.round(cell_best - full[np.arange(n), order[:, 1]], 10), 0.0, 1.0)
    cell_call = np.where(cell_best > 0, flabels[order[:, 0]], UNKNOWN)
    mean = pd.DataFrame(full).groupby(cells["cluster_id"].to_numpy()).mean()
    mord = np.argsort(-rkey(mean.to_numpy()), axis=1, kind="stable")
    top = mean.to_numpy()[np.arange(len(mean)), mord[:, 0]]
    sec = mean.to_numpy()[np.arange(len(mean)), mord[:, 1]]
    clus_marker = pd.DataFrame(
        {"call": flabels[mord[:, 0]], "share": top, "lead": top - sec}, index=mean.index
    )
    per_cell = clus_marker.loc[cells["cluster_id"]].reset_index(drop=True)
    m_ok = ((per_cell["share"] >= cg["share"] - EPS) & (per_cell["lead"] >= cg["lead"] - EPS)).to_numpy()
    take = ~m_ok & (cell_margin >= cg["fallback_margin"] - EPS) & (cell_best > 0)
    marker_label = np.where(m_ok, per_cell["call"].to_numpy(), np.where(take, cell_call, UNKNOWN))
    marker_cl = np.array([fcl.get(x, "") for x in marker_label])

    # Cluster table -------------------------------------------------------------------------
    cluster_ids = list(pd.unique(cells["cluster_id"]))  # order of factorize
    cl_strict, cl_cnt = cluster_plurality(strict, cl, nC, nT)
    cl_guess, cl_gcnt = cluster_plurality(win, cl, nC, nT)
    sizes = np.bincount(cl, minlength=nC)
    cl_term = np.full(nC, -1)
    cluster_rows = []
    tool_major = {tool: cluster_plurality(C[:, m], cl, nC, nT) for m, tool in enumerate(tools)}
    for j in range(nC):
        share_j = cl_cnt[j] / sizes[j]
        ranked = np.argsort(-cl_cnt[j], kind="stable")
        top1, top2 = ranked[0], ranked[1] if nT > 1 else ranked[0]
        s1, s2 = share_j[top1], share_j[top2]
        lead = s1 - s2
        fails = []
        if s1 < cg["share"] - EPS:
            fails.append("cluster_share")
        if lead < cg["lead"] - EPS:
            fails.append("cluster_lead")
        passed = s1 > 0 and not fails
        if cl_cnt[j].sum() == 0:
            fails = ["no_cell_called"]
            passed = False
        if passed:
            cl_term[j] = top1
        tool_calls, tool_terms = {}, []
        for tool in tools:
            tt, tc = tool_major[tool][0][j], tool_major[tool][1][j]
            tool_terms.append(int(tt))
            tool_calls[tool] = {
                "label": V.names[tt] if tt >= 0 else UNKNOWN,
                "cl_id": terms[tt] if tt >= 0 else "",
                "cluster_share": float(fmt(tc[tt] / sizes[j], dg)) if tt >= 0 else 0.0,
            }
        gmean = mean.loc[cluster_ids[j]].to_numpy()
        gtop = np.argsort(-rkey(gmean), kind="stable")[:3]
        marker_votes = [
            {"label": flabels[g], "cl_id": fcl[flabels[g]], "vote_share": float(fmt(gmean[g], dg))}
            for g in gtop
        ]
        cluster_rows.append(
            {
                "j": j, "top1": top1, "top2": top2, "s1": s1, "s2": s2, "lead": lead, "passed": passed,
                "fails": fails, "tool_calls": tool_calls, "tool_terms": tool_terms,
                "marker_votes": marker_votes,
            }
        )
    cell_cluster_term = cl_term[cl]

    # Cluster decision rule and text ----------------------------------------------------
    out_rows = []
    for r in cluster_rows:
        j = r["j"]
        cid = cluster_ids[j]
        label_term = cl_term[j]
        guess_term = label_term if r["passed"] else cl_guess[j]
        members = np.where(cl == j)[0]
        agree_all = len(set(r["tool_terms"])) == 1 and r["tool_terms"][0] >= 0
        if not r["passed"]:
            decision = "abstain"
        elif agree_all and r["tool_terms"][0] == label_term:
            decision = "all_agree"
        else:
            carriers = members[strict[members] == label_term]
            counts = pd.Series(rule_final[carriers]).value_counts()
            decision = sorted(counts.index, key=lambda x: (-counts[x], RULES.index(x)))[0]
        label = V.names[label_term] if label_term >= 0 else UNKNOWN
        guess = V.names[guess_term] if guess_term >= 0 else UNKNOWN
        mk = clus_marker.loc[cid]
        marker_gate_pass = bool(mk["share"] >= cg["share"] - EPS and mk["lead"] >= cg["lead"] - EPS)
        flags = []
        if not agree_all:
            flags.append("tools_disagree")
        if not marker_gate_pass:
            flags.append("marker_gate_fail")
        if label_term >= 0 and mk["call"] != V.names[label_term] and mk["call"] != label:
            flags.append("marker_call_differs")
        tool_txt = ", ".join(f"{t}={v['label']}" for t, v in r["tool_calls"].items())
        top_txt = f"{V.names[r['top1']]} {fmt(r['s1'], dg)}" if r["s1"] > 0 else "no called cell"
        if r["passed"]:
            why = (
                f"{label}: {int(round(r['s1'] * sizes[j]))} of {sizes[j]} cells call it "
                f"(share {fmt(r['s1'], dg)} >= {fmt(cg['share'], dg)}, lead {fmt(r['lead'], dg)} >= "
                f"{fmt(cg['lead'], dg)}). Tool cluster majorities: {tool_txt}. Marker-only call: "
                f"{mk['call']} ({fmt(mk['share'], dg)})."
            )
        else:
            why = (
                f"No label passes the cluster gate (share >= {fmt(cg['share'], dg)}, lead >= "
                f"{fmt(cg['lead'], dg)}). Top cell call: {top_txt}, lead {fmt(r['lead'], dg)}. "
                f"Tool cluster majorities: {tool_txt}. Marker-only call: {mk['call']} "
                f"({fmt(mk['share'], dg)})."
            )
        rejected = [
            {"label": V.names[t], "cl_id": terms[t], "cell_share": float(fmt(cl_cnt[j][t] / sizes[j], dg))}
            for t in np.argsort(-cl_cnt[j], kind="stable")[:4]
            if cl_cnt[j][t] > 0 and t != label_term
        ][:3]
        flagged = 0  # filled below
        out_rows.append(
            {
                "cluster_id": cid,
                "label": label,
                "cl_id": terms[label_term] if label_term >= 0 else "",
                "label_best_guess": guess,
                "cl_id_best_guess": terms[guess_term] if guess_term >= 0 else "",
                "best_guess_low_confidence": bool(not r["passed"] and guess_term >= 0),
                "decision_rule": decision,
                "tool_calls": jdump(r["tool_calls"]),
                "marker_votes": jdump(r["marker_votes"]),
                "marker_only_call": mk["call"],
                "marker_gate_pass": marker_gate_pass,
                "gate_failures": "|".join(r["fails"]),
                "gate_flags": "|".join(flags),
                "rationale": why,
                "rejected_alternatives": jdump(rejected),
                "n_cells": int(sizes[j]),
                "n_cells_flagged": flagged,
                "cluster_cell_share": fmt(r["s1"], dg),
                "cluster_cell_lead": fmt(r["lead"], dg),
            }
        )

    # Per-cell flags -------------------------------------------------------------------------
    cluster_known = cell_cluster_term >= 0
    departs = cluster_known & (strict != cell_cluster_term)
    fallback = ~cluster_known & (strict >= 0)
    tie = np.isin(rule_final, ["tie_break_specificity", "tie_break_marker", "tie_break_tool"])
    reason = np.full(n, "", object)
    for mask, text in (
        (departs & (strict >= 0), "departs_from_cluster_label"),
        (departs & (strict < 0), "abstain_in_labelled_cluster"),
        (fallback, "fallback_call_in_unresolved_cluster"),
        (tie, "tie_break"),
    ):
        reason = np.where(mask, np.where(reason == "", text, reason + "|" + text), reason)
    flagged_mask = reason != ""
    per_cluster_flagged = np.bincount(cl[flagged_mask], minlength=nC)
    for r in out_rows:
        r["n_cells_flagged"] = int(per_cluster_flagged[cluster_ids.index(r["cluster_id"])])
    name_of = lambda idx: np.where(idx >= 0, V.names[np.where(idx >= 0, idx, 0)], UNKNOWN)  # noqa: E731
    id_of = lambda idx: np.where(idx >= 0, np.array(terms)[np.where(idx >= 0, idx, 0)], "")  # noqa: E731
    tool_label = {t: name_of(C[:, m]) for m, t in enumerate(tools)}

    # Write -------------------------------------------------------------------------------
    calls = pd.DataFrame(
        {
            "barcode": cells["barcode"],
            "cluster_id": cells["cluster_id"],
            "label_strict": name_of(strict),
            "cl_id_strict": id_of(strict),
            "label_best_guess": name_of(win),
            "cl_id_best_guess": id_of(win),
            "best_guess_low_confidence": (~called) & (win >= 0),
            "decision_rule": rule_final,
            "winner_share": fmt_col(share, dg),
            "n_tied_terms": ntied if c["rule"] == "h5" else "NA",
            "marker_only_label": marker_label,
            "marker_only_cl_id": marker_cl,
            **{f"tool_{t}": tool_label[t] for t in tools},
        }
    )
    write_tsv(calls, paths[0])
    cons = pd.DataFrame(out_rows)
    cons["_k"] = cons["cluster_id"].map(natural_key)
    cons = cons.sort_values("_k", kind="stable").drop(columns="_k")
    write_tsv(cons, paths[1])
    flags_tbl = pd.DataFrame(
        {
            "barcode": cells["barcode"][flagged_mask].to_numpy(),
            "cluster_id": cells["cluster_id"][flagged_mask].to_numpy(),
            "cell_label": name_of(strict)[flagged_mask],
            "cell_best_guess": name_of(win)[flagged_mask],
            "cluster_label": name_of(cell_cluster_term)[flagged_mask],
            "reason": reason[flagged_mask],
            "decision_rule": rule_final[flagged_mask],
            "winner_share": fmt_col(share[flagged_mask], dg),
            "tool_labels": [
                ";".join(f"{t}={tool_label[t][i]}" for t in tools) for i in np.where(flagged_mask)[0]
            ],
        }
    )
    write_tsv(flags_tbl, paths[2])
    summary = {
        "n_cells": int(n),
        "n_clusters": int(nC),
        "n_called": int(called.sum()),
        "coverage": float(fmt(called.mean(), dg)),
        "n_flagged_cells": int(flagged_mask.sum()),
        "decision_rule_counts": {r_: int((rule_final == r_).sum()) for r_ in RULES},
        "n_clusters_labelled": int((cl_term >= 0).sum()),
        "all_tools_agree_cells": int(all_same.sum()),
    }
    paths[3].write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(
        f"consensus: rule={c['rule']}, coverage {summary['coverage']}, {summary['n_clusters_labelled']}/"
        f"{nC} clusters labelled, {summary['n_flagged_cells']} cells flagged"
    )


# Step: markers ---------------------------------------------------------------------------
def step_markers(ctx):
    cfg = ctx.cfg
    (status_path,) = ctx.claim("work/markers_status.json")
    if not cfg["inputs"].get("markers"):
        status_path.write_text(json.dumps({"status": "skipped, inputs.markers is not set"}) + "\n")
        print("markers: skipped, inputs.markers is not set")
        return
    (path,) = ctx.claim("marker_table.tsv")
    mt = cfg["markers_table"]
    raw = pd.read_csv(resolve(cfg, cfg["inputs"]["markers"]), sep="\t", keep_default_na=False, dtype=str)
    raw = raw[raw["comparator_role"].eq(mt["comparator_role"])].copy()
    raw["rank"] = raw["rank"].astype(int)
    raw = raw[raw["rank"] <= mt["top_n"]]
    methods = ["wilcoxon", "logreg"]
    parts = {}
    for m in methods:
        sub = raw[raw["method"].eq(m)]
        if sub.empty:
            raise PipelineError(f"markers table has no {m} rows for role {mt['comparator_role']!r}")
        parts[m] = sub.set_index(["cluster_id", "feature"])[["rank", "effect", "adjusted_pvalue"]].rename(
            columns=lambda x, m=m: f"{x}_{m}"
        )
    both = parts["wilcoxon"].join(parts["logreg"], how="inner")  # in the top N under BOTH methods
    ref = raw[raw["method"].eq("wilcoxon")].set_index(["cluster_id", "feature"])
    both = both.join(
        ref[["target_fraction_detected", "reference_fraction_detected"]], how="left", validate="1:1"
    ).reset_index()
    both["rank_wilcoxon"] = both["rank_wilcoxon"].astype(int)
    both["rank_logreg"] = both["rank_logreg"].astype(int)
    both["_k"] = both["cluster_id"].map(natural_key)
    both = both.sort_values(["_k", "rank_wilcoxon", "feature"], kind="stable").drop(columns="_k")
    both.insert(2, "comparator_role", mt["comparator_role"])
    for col in both.columns:
        if col.startswith(("effect_", "adjusted_pvalue_")) or col.endswith("_fraction_detected"):
            both[col] = [fmt(float(v), ctx.digits) if v != "" else "NA" for v in both[col]]
    write_tsv(both, path)
    status_path.write_text(
        json.dumps({"status": "written", "n_rows": int(len(both)), "enrichment": "not built"}) + "\n"
    )
    print(f"markers: {len(both)} cluster-gene rows in the top {mt['top_n']} under both methods")
    print(f"markers: {MARKERS_NOT_BUILT}")


# Step: score (the only step that reads truth) ---------------------------------------------
def step_score(ctx):
    cfg = ctx.cfg
    if "scoring" not in cfg:
        raise PipelineError("the score step needs a scoring block in the config")
    sc = cfg["scoring"]
    script, truth = resolve(cfg, sc["script"]), resolve(cfg, sc["truth"])
    for p in (script, truth):
        if not p.is_file():
            raise PipelineError(f"scoring input not found: {p}")
    outdir = ctx.out / "score"
    (table_path,) = ctx.claim("scores.tsv")  # scores.tsv is written last: no table, no finished score
    if outdir.exists():
        if not ctx.owned():
            raise PipelineError(f"{outdir} exists and {ctx.out} was not created by this pipeline; remove it yourself")
        shutil.rmtree(outdir)  # left by an interrupted score step
    (outdir / "preds").mkdir(parents=True)
    calls = pd.read_csv(ctx.out / "calls.tsv", sep="\t", keep_default_na=False, dtype=str)
    tools = cfg["consensus"]["tool_order"]
    mapped = pd.read_csv(ctx.out / "work/tools_mapped.tsv", sep="\t", keep_default_na=False, dtype=str)
    ont = load_ontology(cfg)

    def pred_frame(label, cl):
        return pd.DataFrame(
            {
                "barcode": calls["barcode"],
                "predicted_label": label,
                "predicted_cl_id": cl,
                "confidence": np.where(np.asarray(label) == UNKNOWN, 0.0, 1.0),
            }
        )

    methods = {
        "consensus_strict": pred_frame(calls["label_strict"], calls["cl_id_strict"]),
        "consensus_best_guess": pred_frame(calls["label_best_guess"], calls["cl_id_best_guess"]),
        "marker_only": pred_frame(calls["marker_only_label"], calls["marker_only_cl_id"]),
    }
    for t in tools:
        terms = mapped[f"{t}__raw"]
        label = terms.map(lambda x: ont.resolve(x).name if x else UNKNOWN)
        methods[f"tool_{t}"] = pred_frame(label, terms)
    rows = []
    for name, frame in methods.items():
        p = outdir / "preds" / f"predictions_{name}.tsv"
        frame.to_csv(p, sep="\t", index=False, lineterminator="\n")
        subprocess.run(
            [sys.executable, str(script), str(p), str(truth), str(outdir), "--method", name,
             "--dataset", sc["dataset"]],
            check=True, capture_output=True, cwd=script.parent,
        )
        d = json.loads((outdir / f"{sc['dataset']}.{name}.json").read_text())
        rows.append(
            {
                "method": name,
                "coverage": fmt(d["coverage"], ctx.digits),
                "descendant_accuracy": fmt(d["descendant_accuracy"], ctx.digits),
                "strict_selective_accuracy": fmt(d["selective_accuracy"], ctx.digits),
                "strict_accuracy_all_cells": fmt(d["accuracy"], ctx.digits),
            }
        )
    write_tsv(pd.DataFrame(rows), table_path)
    print("score: " + "; ".join(f"{r['method']} {r['descendant_accuracy']} at {r['coverage']}" for r in rows))


# Step: report -----------------------------------------------------------------------------
HASHED_EXCLUDE = {"run_manifest.json", MARKER}


def hashed_outputs(out):
    """Every output file that is part of the hash comparison, relative paths, sorted."""
    files = []
    for p in sorted(out.rglob("*")):
        rel = p.relative_to(out).as_posix()
        if not p.is_file() or rel in HASHED_EXCLUDE or rel.startswith("score/"):
            continue
        if rel == "reproducibility/checksums.sha256":
            continue
        files.append(rel)
    return files


def step_report(ctx):
    cfg = ctx.cfg
    dg = ctx.digits
    c, cg = cfg["consensus"], cfg["cluster_gate"]
    paths = ctx.claim(
        "report.md", "result.json", "reproducibility/commands.sh", "reproducibility/environment.txt",
        "reproducibility/input_checksums.sha256", "reproducibility/checksums.sha256",
    )
    cons = pd.read_csv(ctx.out / "consensus.tsv", sep="\t", keep_default_na=False, dtype=str)
    summary = json.loads((ctx.out / "work/consensus_summary.json").read_text())
    validation = json.loads((ctx.out / "validation.json").read_text())
    flags = pd.read_csv(ctx.out / "cell_flags.tsv", sep="\t", keep_default_na=False, dtype=str)
    tmap = pd.read_csv(ctx.out / "tool_term_map.tsv", sep="\t", keep_default_na=False, dtype=str)
    ins = input_hashes(cfg)  # key -> (file name, sha256): names only, paths live in run_manifest.json
    params = {
        "consensus": {k: c[k] for k in sorted(c)},
        "cluster_gate": dict(cg),
        "smoothing": dict(cfg["smoothing"]),
        "determinism": {"float_digits": cfg["determinism"]["float_digits"]},  # threads: manifest only
    }

    # report.md ---------------------------------------------------------------------------
    L = [f"# Annotation consensus: {cfg['run']['name']}", ""]
    L.append(f"**Status of the settings: {c['values_status']}.**")
    L.append("")
    L.append(
        f"{summary['n_cells']} cells in {summary['n_clusters']} clusters. Consensus rule `{c['rule']}` "
        f"over {', '.join(c['tool_order'])} plus cell-curator marker votes. "
        f"{summary['n_called']} cells ({fmt(summary['coverage'], dg)}) keep a strict label; "
        f"{summary['n_clusters_labelled']} of {summary['n_clusters']} clusters pass the cluster gate; "
        f"{summary['n_flagged_cells']} cells are flagged."
    )
    L += ["", "## Settings", "", "| key | value |", "| --- | --- |"]
    for k in sorted(c):
        L.append(f"| consensus.{k} | {c[k] if not isinstance(c[k], list) else ', '.join(c[k])} |")
    for blk in ("cluster_gate", "smoothing"):
        for k, v in cfg[blk].items():
            L.append(f"| {blk}.{k} | {v} |")
    L += ["", "## Cluster labels", "",
          "Strict label keeps `Unknown`. Best guess is the plurality of per-cell winners and is low "
          "confidence wherever the strict label is `Unknown`.", "",
          "| cluster | cells | strict label | best guess | low confidence | rule | flagged | flags |",
          "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for _, r in cons.iterrows():
        L.append(
            f"| {r['cluster_id']} | {r['n_cells']} | {r['label']} | {r['label_best_guess']} | "
            f"{r['best_guess_low_confidence']} | {r['decision_rule']} | {r['n_cells_flagged']} | "
            f"{r['gate_flags'] or '-'} |"
        )
    L += ["", "## Why each cluster got its label", ""]
    for _, r in cons.iterrows():
        L.append(f"- Cluster {r['cluster_id']}: {r['rationale']}")
    L += ["", "## Decision rules by cell", "", "| decision_rule | cells |", "| --- | --- |"]
    for k in RULES:
        L.append(f"| {k} | {summary['decision_rule_counts'][k]} |")
    L += ["", f"All tools agree on {summary['all_tools_agree_cells']} cells. Agreement is not proof: "
          "tools that share training data can agree and be wrong.", ""]
    L += ["## Flagged cells", ""]
    if len(flags):
        counts = flags["reason"].value_counts().sort_index()
        L += ["| reason | cells |", "| --- | --- |"] + [f"| {k} | {v} |" for k, v in counts.items()]
        L += ["", "Every flagged cell is in `cell_flags.tsv` with its tool labels."]
    else:
        L.append("No cell is flagged.")
    L += ["", "## Tool label mapping", ""]
    snapped = tmap[tmap["status"].eq("snapped")]
    if len(snapped):
        agg = snapped.groupby(["tool", "name", "mapped_name"], sort=True)["n_cells"].apply(
            lambda s: int(s.astype(int).sum())
        )
        L += ["| tool | tool label | mapped to | cells |", "| --- | --- | --- | --- |"]
        L += [f"| {t} | {a} | {b} | {v} |" for (t, a, b), v in agg.items()]
    else:
        L.append(
            f"No tool label was moved (snap_to_vocabulary is {c['snap_to_vocabulary']}, or every label "
            "is already in the vocabulary)."
        )
    marker_file = ctx.out / "marker_table.tsv"
    L += ["", "## Marker genes", ""]
    if marker_file.exists():
        mtab = pd.read_csv(marker_file, sep="\t", keep_default_na=False, dtype=str)
        per = mtab.groupby("cluster_id").size()
        L.append(
            f"`marker_table.tsv` lists genes in the top {cfg['markers_table']['top_n']} under both "
            f"wilcoxon and logreg against `{cfg['markers_table']['comparator_role']}`: {len(mtab)} "
            f"cluster-gene rows over {len(per)} clusters."
        )
    else:
        L.append("No marker table: `inputs.markers` is not set.")
    L.append("")
    L.append(f"**{MARKERS_NOT_BUILT}**")
    if (ctx.out / "scores.tsv").exists():
        sc = pd.read_csv(ctx.out / "scores.tsv", sep="\t", keep_default_na=False, dtype=str)
        L += ["", "## Scores against truth", "",
              "Descendant accuracy among calls and coverage, from `score_v8.py`. Truth was read only "
              "by this step.", "", "| method | coverage | descendant accuracy |", "| --- | --- | --- |"]
        L += [f"| {r['method']} | {r['coverage']} | {r['descendant_accuracy']} |" for _, r in sc.iterrows()]
    L += ["", "## Evidence and limits", "",
          "table fallback: this report is a provenance record of per-cluster and per-cell decisions, "
          "not an analysis result with a figure. No figure was made, so no conclusion here has an "
          "inspected figure: unsupported: no inspected figure.",
          "",
          "A consensus label is evidence, like any other. It never writes back by itself and it does "
          "not pass any cell-curator gate. No LLM is used and nothing leaves this machine.", ""]
    L += ["## Reproducibility", "",
          "Inputs and outputs are hashed in `reproducibility/`. `annotation_consensus.py verify` reruns "
          "the pipeline into a temporary directory and compares every hashed output byte by byte.", ""]
    paths[0].write_text("\n".join(L) + "\n")

    # reproducibility bundle -----------------------------------------------------------------
    cmds = ["#!/usr/bin/env bash", "# Rerun this analysis. Set CONFIG to the strict config YAML.",
            "# The script pins the thread variables from determinism.threads itself.",
            "set -euo pipefail",
            'CONFIG="${CONFIG:?set CONFIG to the strict config YAML}"',
            'OUT="${OUT:-annotation_consensus_out}"',
            'SKILL_DIR="${SKILL_DIR:-.}"']
    for s in STEPS:
        if s == "score" and "scoring" not in cfg:
            continue
        flag = "" if s != "score" else " --force"
        cmds.append(f'python "$SKILL_DIR/annotation_consensus.py" {s} --input "$CONFIG" --output "$OUT"{flag}')
    paths[2].write_text("\n".join(cmds) + "\n")
    have = installed_versions(read_lock())  # same names as the lock file
    env = [f"python {have.pop('python')}", f"platform {platform.system()} {platform.machine()}"]
    env += [f"{k} {v}" for k, v in sorted(have.items())]
    env.append(f"{LOCK_HASH_KEY} {sha256_file(LOCK_FILE)}")  # the lock file this run was checked against
    paths[3].write_text("\n".join(env) + "\n")  # hashed: verify passes only on the pinned environment
    paths[4].write_text("".join(f"{ins[k][1]}  {k}={ins[k][0]}\n" for k in sorted(ins)))

    # result.json ------------------------------------------------------------------------------
    data_files = [f for f in hashed_outputs(ctx.out) if f not in ("report.md", "result.json")
                  and not f.startswith("reproducibility/")]
    cluster_rows = [
        {
            **{k: r[k] for k in ("cluster_id", "label", "cl_id", "label_best_guess", "decision_rule",
                                 "gate_failures", "gate_flags")},
            "n_cells": int(r["n_cells"]),
            "n_cells_flagged": int(r["n_cells_flagged"]),
        }
        for _, r in cons.iterrows()
    ]
    result = {
        "skill": "annotation-consensus",
        "version": VERSION,
        "run_name": cfg["run"]["name"],
        "parameters": params,
        "inputs": {k: {"file": ins[k][0], "sha256": ins[k][1]} for k in sorted(ins)},
        "version_drift": list(ctx.drift),
        "validation": validation,
        "summary": summary,
        "clusters": cluster_rows,
        "marker_enrichment": "not built",
        "outputs": {f: sha256_file(ctx.out / f) for f in data_files},
    }
    if (ctx.out / "scores.tsv").exists():
        sc = pd.read_csv(ctx.out / "scores.tsv", sep="\t", keep_default_na=False, dtype=str)
        result["scores"] = sc.to_dict(orient="records")
    paths[1].write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")
    lines = [f"{sha256_file(ctx.out / f)}  {f}\n" for f in hashed_outputs(ctx.out)]
    paths[5].write_text("".join(lines))
    write_manifest(ctx)
    print(f"report: wrote report.md, result.json and the reproducibility bundle ({len(lines)} hashed files)")


def write_manifest(ctx):
    """run_manifest.json is not hashed: time, host, paths, threads, input hashes, version drift."""
    cfg = ctx.cfg
    inputs = {
        key: {"path": str(p), "sha256": sha256_file(p)} for key, p in input_files(cfg).items()
    }
    inputs["config"]["sha256_canonical"] = config_sha256(cfg)  # what input_checksums.sha256 holds
    if (ctx.out / "scores.tsv").exists():  # truth is read by the score step only
        truth = resolve(cfg, cfg["scoring"]["truth"])
        inputs["truth"] = {"path": str(truth), "sha256": sha256_file(truth)}
    manifest = {
        "written_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "output_dir": str(ctx.out.resolve()),
        "config": str(cfg["_path"]),
        "argv": ctx.argv,
        "threads": THREADS,
        "env": {v: os.environ.get(v) for v in THREAD_VARS},
        "inputs": inputs,
        "lock_file": str(LOCK_FILE),
        "lock_file_sha256": sha256_file(LOCK_FILE),
        "versions": installed_versions(read_lock()),
        "version_drift": {"allowed": ctx.allow_drift, "differences": list(ctx.drift)},
        "container_sif_sha256": os.environ.get("CONSENSUS_SIF_SHA256"),  # set by container/run_container.sh
        "note": "not hashed: holds time, host, paths, thread count and input hashes",
    }
    (ctx.out / "run_manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


# run and verify ---------------------------------------------------------------------------
STEP_FUNCS = {
    "validate": step_validate,
    "votes": step_votes,
    "smooth": step_smooth,
    "map-tools": step_map_tools,
    "consensus": step_consensus,
    "markers": step_markers,
    "score": step_score,
    "report": step_report,
}
RUN_OUTPUTS = [
    "validation.json", "work", "tool_term_map.tsv", "calls.tsv", "consensus.tsv", "cell_flags.tsv",
    "marker_table.tsv", "scores.tsv", "score", "report.md", "result.json", "reproducibility",
    "run_manifest.json",
]


def run_all(ctx, with_score):
    out = ctx.out
    present = [r for r in RUN_OUTPUTS if (out / r).exists()]
    if present and not ctx.force:
        raise PipelineError(f"refusing to overwrite {present} in {out}; use --force or a new --output")
    if present and not ctx.owned():
        raise PipelineError(
            f"--force refused: {out} holds {present} but no {MARKER} marker, so this pipeline did not "
            "create it; remove those paths yourself or use a new --output"
        )
    for r in present:  # only a directory the pipeline created reaches this line
        p = out / r
        shutil.rmtree(p) if p.is_dir() else p.unlink(missing_ok=True)
    ctx.mark()
    for s in STEPS:
        if s == "score" and not (with_score and "scoring" in ctx.cfg):
            continue
        STEP_FUNCS[s](ctx)


def scrubbed_env(cfg):
    """The only environment a verify rerun gets: PATH and the pinned thread variables."""
    env = {"PATH": os.environ.get("PATH", "")}
    env.update({v: str(cfg["determinism"]["threads"]) for v in THREAD_VARS})
    return env


def read_sums(path):
    return dict(line.split("  ", 1)[::-1] for line in Path(path).read_text().splitlines() if line)


def lock_problems(out, allow_drift):
    """Lines saying the lock file on disk is not the one the run recorded in environment.txt.

    Refused unless --allow-version-drift; allowed, the line is printed and the hash comparison
    of environment.txt still fails, because the rerun records the lock it sees.
    """
    env = Path(out) / "reproducibility/environment.txt"
    recorded = dict(row.split(" ", 1) for row in env.read_text().splitlines() if " " in row).get(LOCK_HASH_KEY)
    now = sha256_file(LOCK_FILE)
    if recorded == now:
        return []
    line = (
        f"requirements.lock changed since the run: recorded {str(recorded)[:12]}, "
        f"on disk {now[:12]} ({LOCK_FILE})"
    )
    if allow_drift:
        print(f"verify: {line} (allowed by --allow-version-drift)", file=sys.stderr)
        return []
    return [line + "; use the lock file of the run, or --allow-version-drift"]


def verify(cfg, out, allow_drift=False):
    """Hash the stored files, rerun in a fresh scrubbed process, and compare both to the record."""
    out = Path(out)
    sums = out / "reproducibility/checksums.sha256"
    if not sums.is_file():
        raise PipelineError(f"{sums} not found; run the pipeline first")
    recorded = read_sums(sums)
    problems = []
    stored = {f: sha256_file(out / f) for f in hashed_outputs(out)}  # the files as they lie on disk
    for f, digest in recorded.items():
        if f not in stored:
            problems.append(f"{f}: recorded but missing on disk")
        elif stored[f] != digest:
            problems.append(f"{f}: stored file edited, recorded {digest[:12]}, on disk {stored[f][:12]}")
    problems += [f"{f}: on disk but not recorded" for f in sorted(set(stored) - set(recorded))]
    problems += lock_problems(out, allow_drift)
    if problems and any(p.startswith("requirements.lock changed") for p in problems):
        raise PipelineError("verify FAILED:\n  " + "\n  ".join(problems))  # no rerun on another lock
    now = input_hashes(cfg)
    for line in (out / "reproducibility/input_checksums.sha256").read_text().splitlines():
        digest, label = line.split("  ", 1)
        key = label.split("=", 1)[0]
        if key not in now or now[key][1] != digest:
            problems.append(f"input changed since the run: {label}")
    with tempfile.TemporaryDirectory(prefix="annotation_consensus_verify_") as tmp:
        rerun = Path(tmp) / "rerun"
        cmd = [sys.executable, str(Path(__file__).resolve()), "run", "--input", str(cfg["_path"]),
               "--output", str(rerun)]  # fmt: skip
        cmd += ["--score"] * (out / "scores.tsv").exists() + ["--allow-version-drift"] * allow_drift
        done = subprocess.run(cmd, env=scrubbed_env(cfg), cwd=tmp, capture_output=True, text=True)
        if done.returncode != 0:
            raise PipelineError(f"verify: the rerun failed:\n{done.stderr.strip()[-2000:]}")
        fresh = {f: sha256_file(rerun / f) for f in hashed_outputs(rerun)}
    for f, digest in recorded.items():
        if fresh.get(f) != digest:
            problems.append(f"{f}: recorded {digest[:12]}, rerun {str(fresh.get(f))[:12]}")
    problems += [f"{f}: produced by the rerun but not recorded" for f in sorted(set(fresh) - set(recorded))]
    if problems:
        raise PipelineError("verify FAILED:\n  " + "\n  ".join(problems))
    print(f"verify: OK, {len(recorded)} stored files and the rerun are byte-identical to the record")


def write_expected(out, path=None):
    """expected_hashes.json from a finished demo run: the committed record of the demo outputs."""
    path = Path(EXPECTED_FILE if path is None else path)
    sums = read_sums(Path(out) / "reproducibility/checksums.sha256")
    doc = {
        "about": "sha256 of every hashed output of the demo run on the pinned environment",
        "config": "examples/demo/config.yaml",
        "lock_file_sha256": sha256_file(LOCK_FILE),
        "files": dict(sorted(sums.items())),
    }
    path.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    return path


def expected_diff(out, path=None):
    """Lines that differ between a demo run and expected_hashes.json; empty when they match."""
    want = json.loads(Path(EXPECTED_FILE if path is None else path).read_text())["files"]
    have = read_sums(Path(out) / "reproducibility/checksums.sha256")
    lines = [f"{f}: expected {want[f][:12]}, got {str(have.get(f))[:12]}" for f in want if have.get(f) != want[f]]
    return lines + [f"{f}: not in expected_hashes.json" for f in sorted(set(have) - set(want))]


# CLI --------------------------------------------------------------------------------------
HELP_EPILOG = (
    "Every impactful choice (rule, marker weight, gate, tie-break delta, smoothing k, "
    "cluster-gate share/lead/fallback margin, snap option, marker scope, tool order) is a "
    "required config key; there are no defaults. " + MARKERS_NOT_BUILT + " "
    "Pinned environment: requirements.lock lists python and package versions, and a run refuses "
    "any other version unless --allow-version-drift. reproducibility/environment.txt is hashed, "
    "so verify passes only on the same pinned environment. Exact kNN distance ties are broken by "
    "barcode rank, so neither thread count nor input row order can change a neighbour. "
    "--force removes outputs only from a "
    "directory this pipeline created (it holds a " + MARKER + " marker). An --output inside a "
    "cell-curator run directory is refused."
)
STEP_HELP = {
    "validate": "check inputs exist, barcodes match by name, labels map, required config keys",
    "votes": "port of the marker-program vote loop; skipped when votes are supplied precomputed",
    "smooth": "kNN smoothing of one-hot top-1 votes over the embedding",
    "map-tools": "snap tool labels to the corpus vocabulary along the ontology",
    "consensus": "tool + marker vote, gate, tie-break; writes calls, consensus, cell_flags",
    "markers": "tidy table of genes in the top N under wilcoxon and logreg (enrichment NOT built yet)",
    "score": "optional: score against truth with score_v8.py (the only step that reads truth)",
    "report": "report.md, result.json and the reproducibility bundle",
    "run": "run every step into an empty output directory",
    "verify": "hash the stored outputs, rerun in a fresh scrubbed subprocess, compare all hashed outputs",
    "expected-hashes": "write expected_hashes.json from a finished demo run in --output",
}


def build_parser():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--input", "--config", dest="input", help="the one strict config YAML")
    common.add_argument("--output", help="output directory (default: run.output_dir of the config)")
    common.add_argument("--demo", action="store_true", help="use the bundled demo config")
    common.add_argument("--force", action="store_true", help="overwrite existing outputs")
    common.add_argument("--score", action="store_true", help="also run the score step (needs truth)")
    common.add_argument(
        "--allow-version-drift", action="store_true",
        help="run on versions that differ from requirements.lock; the drift is recorded in the outputs",
    )
    parser = argparse.ArgumentParser(
        prog="annotation_consensus.py", description=__doc__.split("\n\n")[0], epilog=HELP_EPILOG,
        parents=[common],
    )
    sub = parser.add_subparsers(dest="command")
    for name, text in STEP_HELP.items():
        sp = sub.add_parser(name, help=text, description=text, epilog=HELP_EPILOG, parents=[common])
        if name == "expected-hashes":
            sp.add_argument("--expected-file", help="where to write (default: expected_hashes.json of the skill)")
    return parser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    known = set(STEP_HELP)
    if not any(a in known for a in argv) and not any(a in ("-h", "--help") for a in argv):
        argv = ["run", *argv]
    args = build_parser().parse_args(argv)
    command = args.command or "run"
    try:
        if args.demo:
            args.input = str(DEMO_CONFIG)
        if not args.input:
            raise PipelineError("give --input CONFIG.yaml, or --demo")
        cfg = load_config(args.input)
        drift = check_versions(args.allow_version_drift)
        out = args.output or cfg["run"].get("output_dir")
        if not out:
            raise PipelineError("give --output, or set run.output_dir in the config")
        out = resolve(cfg, out) if not args.output else Path(args.output)
        if command == "verify":
            verify(cfg, out, args.allow_version_drift)
            return 0
        if command == "expected-hashes":
            print(f"expected-hashes: wrote {write_expected(out, args.expected_file)}")
            return 0
        if (run_dir := unsafe_output(out)) is not None:
            raise PipelineError(f"--output {out} is inside the cell-curator run directory {run_dir}; use another path")
        ctx = Ctx(cfg, out, args.force, argv, drift, args.allow_version_drift)
        if command == "run":
            run_all(ctx, args.score)
        else:
            if command != "validate" and not (ctx.out / "validation.json").exists():
                raise PipelineError("run validate first: validation.json is missing")
            ctx.mark()
            STEP_FUNCS[command](ctx)
        return 0
    except PipelineError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
