"""Tests for the annotation-consensus pipeline.

Purpose: check that the demo runs, reruns are byte-identical, bad inputs fail loudly, truth is
         read only by the score step, the provenance tables keep their schema, and (slow) the
         full AIDA run reproduces the frozen c012/c014 numbers.
Run:     .venv/bin/python -m pytest skills/annotation-consensus/tests -q
         add -m "not slow" to skip the full-data regression
"""

import importlib.util
import json
import os
import shutil
import socket
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

SKILL = Path(__file__).resolve().parents[1]
DEMO = SKILL / "examples/demo"
W = Path("/exports/para-lipg-hpc/mdmanurung/cc-v8-work")
CAND = W / "tusoskill/cov95-2026-11-11/candidates"
REPO = SKILL.parents[1]

spec = importlib.util.spec_from_file_location("annotation_consensus", SKILL / "annotation_consensus.py")
ac = importlib.util.module_from_spec(spec)
sys.modules["annotation_consensus"] = ac
spec.loader.exec_module(ac)

# Output columns promised by plans/v8-consensus-mode.md ("Provenance schema").
PLAN_CONSENSUS_COLS = [
    "cluster_id", "label", "cl_id", "decision_rule", "tool_calls", "marker_votes",
    "gate_failures", "gate_flags", "rationale", "rejected_alternatives", "n_cells",
    "n_cells_flagged",
]  # fmt: skip
PLAN_FLAG_COLS = ["barcode", "cluster_id", "reason", "tool_labels"]
DECISION_RULES = {
    "all_agree", "majority", "tie_break_specificity", "tie_break_marker", "tie_break_tool", "abstain",
}  # fmt: skip
RESULT_KEYS = {
    "clusters", "inputs", "marker_enrichment", "outputs", "parameters", "run_name",
    "skill", "summary", "validation", "version",
}  # fmt: skip
DETERMINISTIC_FILES = [
    "calls.tsv", "consensus.tsv", "cell_flags.tsv", "tool_term_map.tsv",
    "work/cells.tsv", "work/tools_mapped.tsv", "work/smoothed_labels.tsv",
    "work/smoothed_votes.npy", "work/consensus_summary.json",
]  # fmt: skip


# Helpers ----------------------------------------------------------------------------------
def run(*argv):
    """Call main() in-process; returns the exit code."""
    return ac.main([str(a) for a in argv])


def read(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str)


def sha(path):
    return ac.sha256_file(Path(path))


def make_cfg(tmp, edit=None, name="config.yaml"):
    """Copy of the demo config with absolute paths; `edit(cfg)` may change it. Returns its path."""
    cfg = yaml.safe_load((DEMO / "config.yaml").read_text())

    def absolute(value):
        return str((DEMO / value).resolve())

    cfg["run"].pop("output_dir", None)
    inp = cfg["inputs"]
    for key in ("cells", "vocabulary", "ontology", "markers"):
        inp[key] = absolute(inp[key])
    inp["tools"] = {t: absolute(p) for t, p in inp["tools"].items()}
    cfg["votes"]["path"] = absolute(cfg["votes"]["path"])
    cfg["scoring"]["script"] = absolute(cfg["scoring"]["script"])
    cfg["scoring"]["truth"] = absolute(cfg["scoring"]["truth"])
    if edit:
        edit(cfg)
    path = Path(tmp) / name
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def copy_input(tmp, src):
    dest = Path(tmp) / Path(src).name
    shutil.copy(src, dest)
    return dest


# Shared runs ------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def demo_a(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo_a") / "out"
    assert run("--demo", "--output", out) == 0
    return out


@pytest.fixture(scope="session")
def demo_b(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo_b") / "out"
    assert run("--demo", "--output", out) == 0
    return out


@pytest.fixture(scope="session")
def demo_scored(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo_scored") / "out"
    assert run("--demo", "--output", out, "--score") == 0
    return out


# 1. The demo runs and writes its documented outputs -----------------------------------------
def test_demo_writes_report_result_and_reproducibility(demo_a):
    assert (demo_a / "report.md").read_text().startswith("#")
    result = json.loads((demo_a / "result.json").read_text())
    assert RESULT_KEYS <= set(result)
    assert result["skill"] == "annotation-consensus"
    assert result["summary"]["n_cells"] == 654 and result["summary"]["n_clusters"] == 23
    assert result["marker_enrichment"] == "not built"
    assert set(result["parameters"]["consensus"]) >= {
        "rule", "marker_weight", "gate", "tie_break_delta", "snap_to_vocabulary", "values_status",
    }  # fmt: skip
    assert "demo values, not approved defaults" in result["parameters"]["consensus"]["values_status"]
    repro = demo_a / "reproducibility"
    for name in ("commands.sh", "checksums.sha256", "environment.txt", "input_checksums.sha256"):
        assert (repro / name).stat().st_size > 0, name
    assert (demo_a / "run_manifest.json").is_file()


def test_demo_cli_subprocess(tmp_path):
    """The real entry point, with threads pinned by the script itself."""
    import subprocess

    out = tmp_path / "cli"
    done = subprocess.run(
        [sys.executable, str(SKILL / "annotation_consensus.py"), "--demo", "--output", str(out)],
        capture_output=True, text=True, env={k: v for k, v in os.environ.items() if k != "OMP_NUM_THREADS"},
    )  # fmt: skip
    assert done.returncode == 0, done.stderr
    assert (out / "report.md").is_file() and (out / "result.json").is_file()
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["env"]["OMP_NUM_THREADS"] == "1"  # the thread count lives in the manifest only
    assert "OMP_NUM_THREADS" not in (out / "reproducibility/environment.txt").read_text()


def test_demo_score_step(demo_scored):
    scores = read(demo_scored / "scores.tsv").set_index("method")
    assert {"consensus_strict", "consensus_best_guess", "marker_only", "tool_popv"} <= set(scores.index)
    strict = scores.loc["consensus_strict"]
    # builder-reported demo numbers on the 654-cell fixture: a smoke check, not a frozen benchmark
    assert float(strict["descendant_accuracy"]) == pytest.approx(0.8248, abs=5e-4)
    assert float(strict["coverage"]) == pytest.approx(0.8991, abs=5e-4)
    assert (demo_scored / "score").is_dir()
    assert "scores" in json.loads((demo_scored / "result.json").read_text())  # only after --score


# 2. Determinism ---------------------------------------------------------------------------
def test_two_runs_have_identical_hashes(demo_a, demo_b):
    sums_a = (demo_a / "reproducibility/checksums.sha256").read_text()
    sums_b = (demo_b / "reproducibility/checksums.sha256").read_text()
    assert sums_a == sums_b
    files = ac.hashed_outputs(demo_a)
    assert files == ac.hashed_outputs(demo_b) and len(files) >= 15
    for rel in files:
        assert sha(demo_a / rel) == sha(demo_b / rel), rel


def test_checksum_file_matches_disk(demo_a):
    recorded = dict(
        line.split("  ", 1)[::-1] for line in (demo_a / "reproducibility/checksums.sha256").read_text().splitlines()
    )
    assert set(recorded) == set(ac.hashed_outputs(demo_a))
    assert all(sha(demo_a / rel) == digest for rel, digest in recorded.items())
    inputs = (demo_a / "reproducibility/input_checksums.sha256").read_text()
    for name in ("cells.h5ad", "predictions_popv.tsv", "cell_votes.tsv", "config="):
        assert name in inputs
    assert "truth.tsv" not in inputs


def test_hashed_outputs_hold_no_host_time_or_temp_path(demo_a):
    host = socket.gethostname()
    for rel in ac.hashed_outputs(demo_a):
        if rel.endswith(".npy"):
            continue
        text = (demo_a / rel).read_text()
        assert str(demo_a) not in text, rel
        assert str(demo_a.parent) not in text, rel
        assert host not in text or len(host) < 4, rel
    manifest = (demo_a / "run_manifest.json").read_text()
    assert host in manifest or "host" in manifest  # host details live only in the manifest


def test_threads_recorded_and_pinned(demo_a):
    manifest = json.loads((demo_a / "run_manifest.json").read_text())
    assert manifest["threads"] == 1
    assert all(manifest["env"][v] == "1" for v in ac.THREAD_VARS)
    assert ac.THREADS == 1
    assert all(os.environ[v] == "1" for v in ac.THREAD_VARS)


def test_verify_passes(demo_a, capsys):
    assert run("verify", "--demo", "--output", demo_a) == 0
    assert "verify: OK" in capsys.readouterr().out


def test_verify_fails_on_tampered_checksum(demo_a, tmp_path, capsys):
    copy = tmp_path / "copy"
    shutil.copytree(demo_a, copy)
    sums = copy / "reproducibility/checksums.sha256"
    lines = sums.read_text().splitlines()
    lines[0] = "0" * 64 + lines[0][64:]
    sums.write_text("\n".join(lines) + "\n")
    assert run("verify", "--demo", "--output", copy) == 2
    assert "verify FAILED" in capsys.readouterr().err


def test_verify_fails_when_an_input_changes(tmp_path, capsys):
    tool = copy_input(tmp_path, DEMO / "predictions_celltypist.tsv")
    cfg = make_cfg(tmp_path, lambda c: c["inputs"]["tools"].update(celltypist=str(tool)))
    out = tmp_path / "out"
    assert run("run", "--input", cfg, "--output", out) == 0
    assert run("verify", "--input", cfg, "--output", out) == 0
    frame = read(tool)
    frame.loc[0, "predicted_label"] = "Unknown"
    frame.loc[0, "predicted_cl_id"] = ""
    frame.to_csv(tool, sep="\t", index=False)
    capsys.readouterr()
    assert run("verify", "--input", cfg, "--output", out) == 2
    assert "input changed since the run" in capsys.readouterr().err


def test_row_order_of_inputs_does_not_matter(demo_a, tmp_path):
    """Barcodes align by name: shuffled tool and vote files give the same calls."""
    votes = copy_input(tmp_path, DEMO / "cell_votes.tsv")
    tool = copy_input(tmp_path, DEMO / "predictions_popv.tsv")
    for path in (votes, tool):
        read(path).sample(frac=1, random_state=3).to_csv(path, sep="\t", index=False)

    def edit(c):
        c["votes"]["path"] = str(votes)
        c["inputs"]["tools"]["popv"] = str(tool)

    cfg = make_cfg(tmp_path, edit)
    out = tmp_path / "out"
    assert run("run", "--input", cfg, "--output", out) == 0
    for rel in DETERMINISTIC_FILES:
        assert sha(out / rel) == sha(demo_a / rel), rel


# 3. Failure modes fail loudly -----------------------------------------------------------------
def _required_paths(schema=ac.SCHEMA, prefix=""):
    for key, spec_ in schema.items():
        path = f"{prefix}{key}"
        if isinstance(spec_, dict):
            if path not in ac.OPTIONAL_BLOCKS:
                yield from _required_paths(spec_, f"{path}.")
        elif spec_[1]:
            yield path


# markers_table is an optional block, but its keys are required whenever the block is present.
BLOCK_KEYS = ["markers_table.top_n", "markers_table.comparator_role"]
REQUIRED = sorted(_required_paths()) + BLOCK_KEYS


def test_required_key_list_covers_every_impactful_choice():
    assert set(ac.IMPACTFUL) <= set(REQUIRED)
    assert not set(BLOCK_KEYS) & set(_required_paths())
    assert len(REQUIRED) >= 25


@pytest.mark.parametrize("dotted", REQUIRED)
def test_missing_required_key_fails_loudly(dotted, tmp_path, capsys):
    def edit(cfg):
        *parents, leaf = dotted.split(".")
        node = cfg
        for p in parents:
            node = node[p]
        del node[leaf]

    cfg = make_cfg(tmp_path, edit)
    out = tmp_path / "out"
    assert run("validate", "--input", cfg, "--output", out) == 2
    err = capsys.readouterr().err
    assert f"missing required key: {dotted}" in err
    if dotted in ac.IMPACTFUL:
        assert "no default" in err
    assert not out.exists()


def test_unknown_key_fails_loudly(tmp_path, capsys):
    cfg = make_cfg(tmp_path, lambda c: c["consensus"].update(typo_gate=0.5))
    assert run("validate", "--input", cfg, "--output", tmp_path / "out") == 2
    assert "unknown key: consensus.typo_gate" in capsys.readouterr().err


def test_a_bad_choice_value_fails_loudly(tmp_path, capsys):
    cfg = make_cfg(tmp_path, lambda c: c["consensus"].update(rule="referee"))
    assert run("validate", "--input", cfg, "--output", tmp_path / "out") == 2
    assert "consensus.rule must be one of" in capsys.readouterr().err


def test_barcode_mismatch_in_a_tool_file(tmp_path, capsys):
    tool = copy_input(tmp_path, DEMO / "predictions_popv.tsv")
    frame = read(tool)
    frame.loc[:4, "barcode"] = frame.loc[:4, "barcode"] + "_renamed"
    frame.to_csv(tool, sep="\t", index=False)
    cfg = make_cfg(tmp_path, lambda c: c["inputs"]["tools"].update(popv=str(tool)))
    out = tmp_path / "out"
    assert run("run", "--input", cfg, "--output", out) == 2
    err = capsys.readouterr().err
    assert "tool popv: 5 cell barcodes missing" in err
    assert not (out / "consensus.tsv").exists()


def test_barcode_mismatch_in_the_vote_file(tmp_path, capsys):
    votes = copy_input(tmp_path, DEMO / "cell_votes.tsv")
    read(votes).iloc[3:].to_csv(votes, sep="\t", index=False)
    cfg = make_cfg(tmp_path, lambda c: c["votes"].update(path=str(votes)))
    out = tmp_path / "out"
    assert run("run", "--input", cfg, "--output", out) == 2
    assert "3 cell barcodes have no vote row" in capsys.readouterr().err
    assert not (out / "consensus.tsv").exists()


@pytest.mark.parametrize("bad_id", ["CL:9999999", "NOTANID:1"])
def test_unmappable_tool_label(bad_id, tmp_path, capsys):
    tool = copy_input(tmp_path, DEMO / "predictions_singler.tsv")
    frame = read(tool)
    frame.loc[frame["predicted_label"].ne("Unknown").idxmax(), "predicted_cl_id"] = bad_id
    frame.to_csv(tool, sep="\t", index=False)
    cfg = make_cfg(tmp_path, lambda c: c["inputs"]["tools"].update(singler=str(tool)))
    out = tmp_path / "out"
    assert run("run", "--input", cfg, "--output", out) == 2
    err = capsys.readouterr().err
    assert "tool singler: labels with no ontology term" in err and bad_id in err
    assert not (out / "consensus.tsv").exists()


def test_existing_output_directory_needs_force(demo_a, tmp_path, capsys):
    out = tmp_path / "out"
    shutil.copytree(demo_a, out)
    before = {rel: sha(out / rel) for rel in ac.hashed_outputs(out)}
    assert run("--demo", "--output", out) == 2
    assert "refusing to overwrite" in capsys.readouterr().err
    assert before == {rel: sha(out / rel) for rel in ac.hashed_outputs(out)}
    assert run("validate", "--demo", "--output", out) == 2  # single steps refuse too
    assert "refusing to overwrite" in capsys.readouterr().err
    assert run("--demo", "--output", out, "--force") == 0
    assert before == {rel: sha(out / rel) for rel in ac.hashed_outputs(out)}


# 4. Truth is read only by the score step ---------------------------------------------------
def test_truth_is_not_read_before_scoring(demo_a, tmp_path):
    scrambled = tmp_path / "truth_scrambled.tsv"
    truth = read(DEMO / "truth.tsv")
    truth["cell_type"] = truth["cell_type"].sample(frac=1, random_state=1).to_numpy()
    truth.to_csv(scrambled, sep="\t", index=False)
    runs = {"missing": tmp_path / "no_such_truth.tsv", "scrambled": scrambled}
    for name, truth_path in runs.items():
        cfg = make_cfg(tmp_path, lambda c, p=truth_path: c["scoring"].update(truth=str(p)), f"{name}.yaml")
        out = tmp_path / name
        assert run("run", "--input", cfg, "--output", out) == 0, name
        for rel in DETERMINISTIC_FILES:
            assert sha(out / rel) == sha(demo_a / rel), f"{name}: {rel}"
        assert not (out / "scores.tsv").exists()


def test_scoring_does_not_change_the_consensus_outputs(demo_a, demo_scored):
    for rel in DETERMINISTIC_FILES:
        assert sha(demo_a / rel) == sha(demo_scored / rel), rel
    assert not any("truth" in line for line in (demo_scored / "reproducibility/input_checksums.sha256").read_text().splitlines())


def test_score_step_without_truth_fails_loudly(tmp_path, capsys):
    cfg = make_cfg(tmp_path, lambda c: c["scoring"].update(truth=str(tmp_path / "gone.tsv")))
    out = tmp_path / "out"
    assert run("run", "--input", cfg, "--output", out, "--score") == 2
    assert "scoring input not found" in capsys.readouterr().err


# 5. Provenance schema ---------------------------------------------------------------------
@pytest.fixture(scope="module")
def tables(demo_a):
    return {
        "cons": read(demo_a / "consensus.tsv"),
        "flags": read(demo_a / "cell_flags.tsv"),
        "calls": read(demo_a / "calls.tsv"),
    }


def test_consensus_columns_and_one_row_per_cluster(tables):
    cons = tables["cons"]
    assert cons.columns[: len(PLAN_CONSENSUS_COLS)].tolist() == PLAN_CONSENSUS_COLS[:3] + [
        c for c in cons.columns[3 : len(PLAN_CONSENSUS_COLS)]
    ]
    assert set(PLAN_CONSENSUS_COLS) <= set(cons.columns)
    import anndata as ad

    clusters = set(ad.read_h5ad(DEMO / "cells.h5ad").obs["bench_cluster"].astype(str))
    assert len(clusters) == 23
    assert cons["cluster_id"].is_unique and set(cons["cluster_id"]) == clusters


def test_every_cluster_has_a_decision_rule(tables):
    rules = tables["cons"]["decision_rule"]
    assert (rules != "").all()
    assert set(rules) <= DECISION_RULES
    unknown = tables["cons"]["label"] == "Unknown"
    assert unknown.any() and set(rules[unknown]) == {"abstain"}
    for col in ("tool_calls", "marker_votes", "rejected_alternatives"):
        for value in tables["cons"][col]:
            json.loads(value)  # every JSON column parses
    assert (tables["cons"]["rationale"] != "").all()


def test_cell_decision_rules_are_known(tables):
    assert (tables["calls"]["decision_rule"] != "").all()
    assert set(tables["calls"]["decision_rule"]) <= DECISION_RULES


def test_cell_flags_hold_only_flagged_cells(tables):
    flags, calls, cons = tables["flags"], tables["calls"], tables["cons"]
    assert PLAN_FLAG_COLS == [c for c in flags.columns if c in PLAN_FLAG_COLS]
    assert flags["barcode"].is_unique and set(flags["barcode"]) <= set(calls["barcode"])
    assert (flags["reason"] != "").all()
    assert len(flags) == cons["n_cells_flagged"].astype(int).sum()
    per_cluster = flags.groupby("cluster_id").size()
    expected = cons.set_index("cluster_id")["n_cells_flagged"].astype(int)
    assert (per_cluster.reindex(expected.index, fill_value=0) == expected).all()
    # a cell outside the flag file carries its cluster's label and took no tie-break
    merged = calls.merge(cons[["cluster_id", "label"]], on="cluster_id", validate="m:1")
    quiet = merged[~merged["barcode"].isin(flags["barcode"])]
    assert (quiet["label_strict"] == quiet["label"]).all()
    assert not quiet["decision_rule"].str.startswith("tie_break").any()
    # the reverse: a call that departs from its cluster, or a tie-break, is flagged
    departs = merged[merged["label_strict"] != merged["label"]]
    assert set(departs["barcode"]) <= set(flags["barcode"])
    ties = merged[merged["decision_rule"].str.startswith("tie_break")]
    assert set(ties["barcode"]) <= set(flags["barcode"])


def test_strict_keeps_unknown_and_best_guess_has_no_empty_value(tables):
    calls, cons = tables["calls"], tables["cons"]
    assert (calls["label_strict"] == "Unknown").sum() > 0
    assert (calls["label_strict"] != "").all() and (cons["label"] != "").all()
    assert (calls["label_best_guess"] != "").all() and (cons["label_best_guess"] != "").all()
    assert (calls["label_best_guess"] != "Unknown").all()
    assert (cons["label_best_guess"] != "Unknown").all()
    low = calls["best_guess_low_confidence"] == "True"
    assert set(calls["best_guess_low_confidence"]) <= {"True", "False"}
    assert (low == (calls["label_strict"] == "Unknown")).all()
    low_c = cons["best_guess_low_confidence"] == "True"
    assert (low_c == (cons["label"] == "Unknown")).all()
    called = ~low
    assert (calls.loc[called, "label_strict"] == calls.loc[called, "label_best_guess"]).all()


def test_calls_cover_every_demo_cell_by_name(tables):
    import anndata as ad

    names = list(ad.read_h5ad(DEMO / "cells.h5ad").obs_names)
    assert sorted(tables["calls"]["barcode"]) == sorted(names)


# 7. Determinism contract (plans/v8-consensus-mode.md, rows 1-11) -------------------------------
import importlib.util as _ilu  # noqa: E402
import random  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
import types  # noqa: E402

import numpy as np  # noqa: E402

CLEAN_ENV = {k: v for k, v in os.environ.items() if k not in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
EXPECTED = SKILL / "expected_hashes.json"


def cli(cfg, out, *extra, cwd=None, env=None):
    """The real entry point in a fresh process (threads are pinned at import, so this is the only way
    to see a different thread count)."""
    return subprocess.run(
        [sys.executable, str(SKILL / "annotation_consensus.py"), "run", "--input", str(cfg), "--output", str(out), *extra],
        capture_output=True, text=True, cwd=cwd, env=CLEAN_ENV if env is None else env,
    )  # fmt: skip


def hashed_bytes(out):
    return {rel: (Path(out) / rel).read_bytes() for rel in ac.hashed_outputs(Path(out))}


def lock_with(tmp, **edit):
    """Copy of requirements.lock with one version changed; returns its path."""
    lines = []
    for line in (SKILL / "requirements.lock").read_text().splitlines():
        if line.startswith("# python ") and "python" in edit:
            line = f"# python {edit['python']}"
        elif "==" in line and line.split("==")[0] in edit:
            line = f"{line.split('==')[0]}=={edit[line.split('==')[0]]}"
        lines.append(line)
    path = Path(tmp) / "drift.lock"
    path.write_text("\n".join(lines) + "\n")
    return path


# Row 1: pinned environment ------------------------------------------------------------------
def test_lock_file_matches_this_environment():
    assert ac.version_differences() == []
    pinned = ac.read_lock()
    assert {"python", "numpy", "pandas", "scipy", "scikit-learn", "anndata", "pyyaml", "cell-curator"} <= set(pinned)


@pytest.mark.parametrize("pkg", ["numpy", "python"])
def test_version_drift_is_refused(pkg, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(ac, "LOCK_FILE", lock_with(tmp_path, **{pkg: "0.0.1"}))
    out = tmp_path / "out"
    assert run("--demo", "--output", out) == 2
    err = capsys.readouterr().err
    assert f"{pkg}: lock 0.0.1" in err and "--allow-version-drift" in err
    assert not (out / "validation.json").exists()  # refused before any step ran


def test_allow_version_drift_is_recorded(tmp_path, monkeypatch):
    monkeypatch.setattr(ac, "LOCK_FILE", lock_with(tmp_path, numpy="0.0.1"))
    out = tmp_path / "out"
    assert run("--demo", "--output", out, "--allow-version-drift") == 0
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["version_drift"]["allowed"] is True
    assert any(d.startswith("numpy: lock 0.0.1") for d in manifest["version_drift"]["differences"])
    result = json.loads((out / "result.json").read_text())
    assert any(d.startswith("numpy: lock 0.0.1") for d in result["version_drift"])


def test_verify_fails_when_the_lock_file_changes_after_the_run(demo_a, tmp_path, monkeypatch, capsys):
    lock_sha = sha(SKILL / "requirements.lock")
    env = (demo_a / "reproducibility/environment.txt").read_text()
    assert f"{ac.LOCK_HASH_KEY} {lock_sha}" in env  # recorded in a hashed output
    assert json.loads((demo_a / "run_manifest.json").read_text())["lock_file_sha256"] == lock_sha
    changed = tmp_path / "changed.lock"
    changed.write_text((SKILL / "requirements.lock").read_text() + "# one more line\n")
    monkeypatch.setattr(ac, "LOCK_FILE", changed)
    assert run("verify", "--demo", "--output", demo_a) == 2
    err = capsys.readouterr().err
    assert "requirements.lock changed since the run" in err and lock_sha[:12] in err
    # allowed: the warning is printed; a real edited lock would also change the rerun's environment.txt
    assert run("verify", "--demo", "--output", demo_a, "--allow-version-drift") == 0
    assert "allowed by --allow-version-drift" in capsys.readouterr().err


# Row 2: every input is hashed in the manifest ------------------------------------------------
def test_manifest_holds_input_hashes(demo_a, demo_scored):
    manifest = json.loads((demo_a / "run_manifest.json").read_text())
    files = {
        "cells": "cells.h5ad", "vocabulary": "vocabulary.tsv", "votes": "cell_votes.tsv", "markers": "markers.tsv",
        "tool:celltypist": "predictions_celltypist.tsv", "tool:popv": "predictions_popv.tsv",
        "tool:singler": "predictions_singler.tsv",
    }  # fmt: skip
    got = manifest["inputs"]
    assert set(got) == {*files, "config", "ontology"}
    for key, name in files.items():
        assert got[key]["sha256"] == sha(DEMO / name), key
    assert got["config"]["sha256"] == sha(DEMO / "config.yaml")
    assert got["ontology"]["sha256"] == sha(Path(got["ontology"]["path"])) and got["ontology"]["path"].endswith("cl.obo")
    assert "truth" not in got  # truth is hashed only when the score step ran
    scored = json.loads((demo_scored / "run_manifest.json").read_text())["inputs"]
    assert scored["truth"]["sha256"] == sha(DEMO / "truth.tsv")
    assert manifest["lock_file_sha256"] == sha(SKILL / "requirements.lock")


def test_one_changed_byte_changes_the_manifest(demo_a, tmp_path):
    tool = copy_input(tmp_path, DEMO / "predictions_popv.tsv")
    tool.write_bytes(tool.read_bytes() + b"\n")
    cfg = make_cfg(tmp_path, lambda c: c["inputs"]["tools"].update(popv=str(tool)))
    assert run("run", "--input", cfg, "--output", tmp_path / "out") == 0
    before = json.loads((demo_a / "run_manifest.json").read_text())["inputs"]
    after = json.loads((tmp_path / "out/run_manifest.json").read_text())["inputs"]
    assert after["tool:popv"]["sha256"] != before["tool:popv"]["sha256"]
    assert after["tool:singler"]["sha256"] == before["tool:singler"]["sha256"]
    assert "popv=" in (tmp_path / "out/reproducibility/input_checksums.sha256").read_text()


# Row 3: no global RNG ------------------------------------------------------------------------
def test_no_global_random_state():
    code = (SKILL / "annotation_consensus.py").read_text()
    assert not re.search(r"np\.random\.|numpy\.random|^\s*import random|random\.(seed|random|shuffle)", code, re.M)


# Row 4: the thread count does not reach a hashed file ----------------------------------------
def test_hashed_outputs_identical_at_threads_1_and_8(tmp_path):
    results = {}
    for n in (1, 8):
        cfg = make_cfg(tmp_path, lambda c, n=n: c["determinism"].update(threads=n), name=f"t{n}.yaml")
        done = cli(cfg, tmp_path / f"out{n}")
        assert done.returncode == 0, done.stderr
        results[n] = tmp_path / f"out{n}"
    manifests = {n: json.loads((out / "run_manifest.json").read_text()) for n, out in results.items()}
    assert [manifests[n]["threads"] for n in (1, 8)] == [1, 8]
    assert manifests[8]["env"]["OMP_NUM_THREADS"] == "8"  # it was really pinned to 8
    one, eight = hashed_bytes(results[1]), hashed_bytes(results[8])
    assert one.keys() == eight.keys() and len(one) >= 15
    assert [rel for rel in one if one[rel] != eight[rel]] == []
    assert (results[1] / "reproducibility/checksums.sha256").read_text() == (
        results[8] / "reproducibility/checksums.sha256"
    ).read_text()


def test_knn_ties_break_by_index_and_do_not_depend_on_row_order_of_equal_points():
    grid = np.array([[x, y] for x in range(6) for y in range(6)], dtype=float)  # many exact distance ties
    m = 7
    d2 = ((grid[:, None, :] - grid[None, :, :]) ** 2).sum(-1)
    want = np.argsort(d2, axis=1, kind="stable")[:, :m]  # distance, then index
    assert (ac.knn_indices(grid, m) == want).all()
    dup = np.vstack([grid, grid[:5]])  # exact duplicates of five points
    d2 = ((dup[:, None, :] - dup[None, :, :]) ** 2).sum(-1)
    assert (ac.knn_indices(dup, m) == np.argsort(d2, axis=1, kind="stable")[:, :m]).all()


def duplicate_cells(tmp, shuffle_seed=None):
    """Demo cells with X_pca in groups of 5 identical rows (exact kNN ties); optionally row-shuffled."""
    import anndata as ad

    adata = ad.read_h5ad(DEMO / "cells.h5ad")
    base = np.random.default_rng(0).normal(size=(len(adata) // 5 + 1, 50)).astype(np.float32)
    adata.obsm["X_pca"] = base[np.arange(len(adata)) // 5]
    if shuffle_seed is not None:
        adata = adata[np.random.default_rng(shuffle_seed).permutation(len(adata))].copy()
    path = Path(tmp) / "cells.h5ad"
    adata.write_h5ad(path)
    return path


def smooth_with_threads(tmp, cells, threads):
    """validate + smooth in a fresh process; returns (smoothed votes, barcodes) in the file's row order."""
    tmp = Path(tmp)
    tmp.mkdir()
    cfg = make_cfg(tmp, lambda c: (c["inputs"].update(cells=str(cells)), c["determinism"].update(threads=threads)))
    out = tmp / "out"
    for step in ("validate", "smooth"):
        done = subprocess.run(
            [sys.executable, str(SKILL / "annotation_consensus.py"), step, "--input", str(cfg), "--output", str(out)],
            capture_output=True, text=True, env=CLEAN_ENV,
        )  # fmt: skip
        assert done.returncode == 0, done.stderr
    return np.load(out / "work/smoothed_votes.npy"), read(out / "work/cells.tsv")["barcode"].to_numpy()


def test_smoothed_votes_with_exact_duplicate_embeddings_ignore_row_order_and_threads(tmp_path):
    for sub in ("a", "b"):
        (tmp_path / sub).mkdir()
    plain = duplicate_cells(tmp_path / "a")
    mixed = duplicate_cells(tmp_path / "b", shuffle_seed=5)
    v1, n1 = smooth_with_threads(tmp_path / "t1", plain, 1)
    v8, n8 = smooth_with_threads(tmp_path / "t8", plain, 8)
    vs, ns = smooth_with_threads(tmp_path / "ts", mixed, 8)
    assert (n1 == n8).all() and v1.tobytes() == v8.tobytes()  # threads 1 and 8
    assert not (n1 == ns).all()  # the shuffle really moved rows
    aligned = pd.DataFrame(vs, index=ns).loc[n1].to_numpy()
    assert aligned.tobytes() == v1.tobytes()  # same votes per barcode


# Row 5: ties break by label name, after rounding to 1e-9 -------------------------------------
class StubOntology:
    names = {"CL:0000001": "Zeta cell", "CL:0000009": "Alpha cell", "CL:0000005": "Mid cell"}

    def resolve(self, term):
        return types.SimpleNamespace(name=self.names[term], accession=term)

    def ancestors(self, term):
        return []


def test_vote_tie_breaks_by_label_name_not_by_cl_id_and_survives_shuffling():
    ids = list(StubOntology.names)  # the lowest CL id is the LAST label name
    rng = random.Random(7)
    for _ in range(6):
        shuffled = ids[:]
        rng.shuffle(shuffled)
        terms = ac.order_terms(StubOntology(), shuffled)
        assert [StubOntology.names[t] for t in terms] == ["Alpha cell", "Mid cell", "Zeta cell"]
        voting = ac.Voting(StubOntology(), terms, 1.0)
        n = 5
        order = np.array(rng.sample(range(n), n))  # shuffled cell rows
        C = np.full((n, 2), -1)  # both tools abstain: only the marker vote is left
        S = np.zeros((n, 3))
        S[:, [0, 2]] = 0.5  # Alpha and Zeta tie exactly
        S[1, 2] += 4e-10  # a near-tie below 1e-9 is still a tie
        S[2] = [0.5, 0.0, 0.5 + 2e-9]  # a real lead above 1e-9 wins
        win, _, _ = voting.plain(C[order], S[order])
        assert [terms[w] for w in win] == [
            ["CL:0000009", "CL:0000009", "CL:0000001", "CL:0000009", "CL:0000009"][i] for i in order
        ]


def test_cluster_plurality_tie_goes_to_the_first_label_name():
    term = np.array([2, 0, 0, 2])  # Zeta (index 2) and Alpha (index 0) each have two cells
    best, counts = ac.cluster_plurality(term, np.zeros(4, dtype=int), 1, 3)
    assert best[0] == 0 and counts[0].tolist() == [2, 0, 2]


# Row 6: fixed float formatting in every table --------------------------------------------------
FLOAT = re.compile(r"^-?\d+\.\d+(?:[eE][+-]?\d+)?$|^-?\d+[eE][+-]?\d+$")


def float_tokens(text):
    for token in re.split(r"[\s\"',:;()\[\]{}|=]+", text):
        token = token.rstrip(".")
        if FLOAT.match(token):
            yield token


def sig_digits(token):
    mantissa = re.split("[eE]", token.lstrip("-"))[0].replace(".", "").strip("0")
    return len(mantissa)


def test_every_float_in_every_hashed_table_has_at_most_the_configured_digits(demo_a, demo_scored):
    digits = yaml.safe_load((DEMO / "config.yaml").read_text())["determinism"]["float_digits"]
    seen = {}
    for out in (demo_a, demo_scored):
        for rel in ac.hashed_outputs(out):
            if rel.endswith(".npy") or rel == "reproducibility/environment.txt":
                continue
            for token in float_tokens((out / rel).read_text()):
                seen[rel] = seen.get(rel, 0) + 1
                assert sig_digits(token) <= digits, f"{rel}: {token} has more than {digits} significant digits"
    assert seen["marker_table.tsv"] > 1000  # the check is not vacuous: marker values are floats
    assert "scores.tsv" in seen and "calls.tsv" in seen


# Row 7: no checkout path, host, user, time or temp path in a hashed file ------------------------
def test_hashed_files_are_the_same_from_another_cwd_dir_name_user_and_config_location(demo_a, tmp_path):
    elsewhere = tmp_path / "somewhere else"
    elsewhere.mkdir()
    cfg = make_cfg(tmp_path, name="a different config name.yaml")  # absolute paths, another place
    env = {**CLEAN_ENV, "USER": "nobody-else", "LOGNAME": "nobody-else", "HOSTNAME": "another-host", "TMPDIR": str(tmp_path)}
    out = tmp_path / "output named differently"
    done = cli(cfg, out, cwd=elsewhere, env=env)
    assert done.returncode == 0, done.stderr
    a, b = hashed_bytes(demo_a), hashed_bytes(out)
    assert a.keys() == b.keys()
    assert [rel for rel in a if a[rel] != b[rel]] == []


def test_hashed_files_hold_no_path_host_user_or_time(demo_a):
    forbidden = [str(REPO), str(demo_a.parent), "..", socket.gethostname(), os.environ.get("USER", "\0"), "/exports", "/tmp"]
    stamp = re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}")
    for rel in ac.hashed_outputs(demo_a):
        if rel.endswith(".npy"):
            continue
        text = (demo_a / rel).read_text()
        for bad in forbidden:
            assert bad not in text, f"{rel} holds {bad!r}"
        assert not stamp.search(text), rel
    manifest = (demo_a / "run_manifest.json").read_text()
    assert socket.gethostname() in manifest and str(demo_a) in manifest  # they live here instead


# Row 9: Snakemake and the CLI give the same hashed files -----------------------------------------
@pytest.mark.skipif(_ilu.find_spec("snakemake") is None, reason="snakemake is not installed")
def test_snakemake_and_cli_give_identical_hashed_outputs(tmp_path):
    cfg = make_cfg(tmp_path, lambda c: c["run"].update(output_dir="snakemake_out"))
    done = subprocess.run(
        [sys.executable, "-m", "snakemake", "-s", str(SKILL / "Snakefile"), "--configfile", str(cfg), "-c1",
         "--directory", str(tmp_path), "--quiet"],
        capture_output=True, text=True, env=CLEAN_ENV,
    )  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    cli_run = cli(cfg, tmp_path / "cli_out")
    assert cli_run.returncode == 0, cli_run.stderr
    a, b = hashed_bytes(tmp_path / "snakemake_out"), hashed_bytes(tmp_path / "cli_out")
    assert a.keys() == b.keys() and len(a) >= 15
    assert [rel for rel in a if a[rel] != b[rel]] == []
    assert (tmp_path / "snakemake_out/run_manifest.json").is_file()  # the DAG writes it too


# Row 10: verify reruns in a scrubbed subprocess and hashes the stored files ----------------------
def test_verify_reruns_in_a_subprocess_with_a_scrubbed_environment(demo_a, monkeypatch, capsys):
    calls = []
    real = subprocess.run

    def spy(cmd, **kw):
        calls.append((cmd, kw))
        return real(cmd, **kw)

    monkeypatch.setattr(ac.subprocess, "run", spy)
    monkeypatch.setenv("PYTHONHASHSEED", "12345")
    monkeypatch.setenv("SOME_USER_SETTING", "x")
    assert run("verify", "--demo", "--output", demo_a) == 0
    assert "verify: OK" in capsys.readouterr().out
    (cmd, kw), = calls
    assert cmd[0] == sys.executable and "run" in cmd and Path(cmd[cmd.index("--output") + 1]) != demo_a
    assert set(kw["env"]) == {"PATH", *ac.THREAD_VARS}
    assert {kw["env"][v] for v in ac.THREAD_VARS} == {"1"}
    assert ac.scrubbed_env({"determinism": {"threads": 4}})["OMP_NUM_THREADS"] == "4"


def test_verify_does_not_inherit_the_callers_environment(demo_a, tmp_path, monkeypatch):
    poison = tmp_path / "poison"
    poison.mkdir()
    (poison / "anndata.py").write_text("raise RuntimeError('the caller PYTHONPATH leaked into the rerun')\n")
    monkeypatch.setenv("PYTHONPATH", str(poison))
    assert run("verify", "--demo", "--output", demo_a) == 0


def test_verify_fails_when_a_stored_output_file_is_edited(demo_a, tmp_path, capsys):
    copy = tmp_path / "copy"
    shutil.copytree(demo_a, copy)
    with open(copy / "calls.tsv", "a") as fh:
        fh.write("tampered\n")  # the recorded checksums are left alone
    assert run("verify", "--demo", "--output", copy) == 2
    err = capsys.readouterr().err
    assert "calls.tsv: stored file edited" in err
    assert "rerun" not in err.split("calls.tsv")[0]  # the rerun itself still matches the record


# Row 11: the committed expected hashes of the demo -----------------------------------------------
def test_expected_hashes_match_the_demo(demo_a):
    assert ac.expected_diff(demo_a, EXPECTED) == []
    doc = json.loads(EXPECTED.read_text())
    assert doc["lock_file_sha256"] == sha(SKILL / "requirements.lock")
    assert set(doc["files"]) == set(ac.hashed_outputs(demo_a))


def test_expected_hashes_file_is_what_the_pipeline_writes(demo_a, tmp_path):
    again = ac.write_expected(demo_a, tmp_path / "expected.json")
    assert again.read_bytes() == EXPECTED.read_bytes()


def test_a_changed_rule_fails_the_expected_hashes(tmp_path):
    cfg = make_cfg(tmp_path, lambda c: c["consensus"].update(gate=0.9))
    out = tmp_path / "changed"
    assert run("run", "--input", cfg, "--output", out) == 0
    diff = ac.expected_diff(out, EXPECTED)
    assert any(line.startswith("calls.tsv") for line in diff) and any(line.startswith("consensus.tsv") for line in diff)


# Output directory safety ----------------------------------------------------------------------------
def test_force_refuses_a_directory_the_pipeline_did_not_create(tmp_path, capsys):
    mine = tmp_path / "mine"
    (mine / "work").mkdir(parents=True)
    (mine / "work/precious.txt").write_text("keep me\n")
    assert run("--demo", "--output", mine, "--force") == 2
    assert "--force refused" in capsys.readouterr().err
    assert (mine / "work/precious.txt").read_text() == "keep me\n"


def test_force_replaces_the_outputs_of_a_directory_the_pipeline_created(demo_a, tmp_path):
    copy = tmp_path / "copy"
    shutil.copytree(demo_a, copy)
    (copy / "work/stale.txt").write_text("left over\n")
    assert run("--demo", "--output", copy, "--force") == 0
    assert not (copy / "work/stale.txt").exists() and (copy / ac.MARKER).is_file()
    assert (copy / "reproducibility/checksums.sha256").read_text() == (demo_a / "reproducibility/checksums.sha256").read_text()
    assert ac.MARKER not in (copy / "reproducibility/checksums.sha256").read_text()


@pytest.mark.parametrize("make", [lambda r: (r / "run_state.json").write_text("{}"), lambda r: (r / "annotation/levels").mkdir(parents=True)])
def test_output_inside_a_cell_curator_run_directory_is_refused(make, tmp_path, capsys):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    make(run_dir)
    assert run("--demo", "--output", run_dir / "consensus_out") == 2
    assert "cell-curator run directory" in capsys.readouterr().err
    assert not (run_dir / "consensus_out").exists()


def test_help_states_the_environment_and_tie_rules(capsys):
    with pytest.raises(SystemExit):
        ac.main(["--help"])
    text = " ".join(capsys.readouterr().out.split())
    assert "environment.txt is hashed" in text and "only on the same pinned environment" in text
    assert "kNN distance ties are broken by barcode rank" in text


# 6. Full-data regression against the frozen numbers (slow) ----------------------------------
needs_data = pytest.mark.skipif(not W.is_dir(), reason=f"{W} is absent")


def frozen_row(table, **where):
    frame = pd.read_csv(table, sep="\t")
    for key, value in where.items():
        frame = frame[frame[key] == value]
    assert len(frame) == 1, where
    return frame.iloc[0]


def aida_config(tmp, rule, snap):
    dataset = W / "prep-aida_split_a-2026-11-09"

    def edit(cfg):
        cfg["inputs"].update(
            cells=str(dataset / "input.h5ad"),
            vocabulary=str(W / "runs/v8-dev-aida-split-a-003/plan/vocabulary.tsv"),
            tools={
                "celltypist": str(W / "mvp-aida-2026-11-10/predictions_celltypist.tsv"),
                "popv": str(W / "mvp-aida-2026-11-10/predictions_popv.tsv"),
                "singler": str(CAND / "c012_harmonizer/preds/aida_predictions_singler.tsv"),
            },
        )
        cfg["inputs"].pop("markers")
        cfg.pop("markers_table")
        cfg["votes"]["path"] = str(CAND / "c004_celllevel/cell_votes.tsv")
        cfg["consensus"].update(
            rule=rule, snap_to_vocabulary=snap, values_status="frozen c012/c014 values, not approved defaults"
        )
        cfg["scoring"].update(truth=str(dataset / "scoring/truth.tsv"), dataset="aida_split_a")

    return make_cfg(tmp, edit)


@pytest.mark.slow
@needs_data
@pytest.mark.parametrize(
    ("rule", "snap", "frozen", "where"),
    [
        ("h4", False, CAND / "c012_harmonizer/frozen.tsv",
         {"run": "aida3", "variant": "H4cell", "w_mode": "equal", "rule": "R1"}),
        ("h5", True, CAND / "c014_tiebreak/frozen.tsv",
         {"run": "aida3", "variant": "H5cell"}),
    ],
    ids=["c012_H4cell_snap_off", "c014_H5cell_snap_on"],
)  # fmt: skip
def test_aida_split_a_reproduces_frozen_rows(rule, snap, frozen, where, tmp_path):
    cfg = aida_config(tmp_path, rule, snap)
    out = tmp_path / "out"
    assert run("run", "--input", cfg, "--output", out, "--score") == 0
    got = json.loads((out / "score/aida_split_a.consensus_strict.json").read_text())
    row = frozen_row(frozen, **where)
    assert got["coverage"] == pytest.approx(row["coverage"], abs=1e-6)
    assert got["descendant_accuracy"] == pytest.approx(row["desc_accuracy"], abs=1e-6)


@pytest.mark.slow
@needs_data
@pytest.mark.skipif(
    os.environ.get("ANNOTATION_CONSENSUS_AIFI") != "1",
    reason="AIFI (339,530 cells, 9 GB h5ad, full brute kNN) is opt-in: set ANNOTATION_CONSENSUS_AIFI=1",
)
def test_aifi_reproduces_frozen_h4_row(tmp_path):
    dataset = W / "prep-aifi_dev-2026-11-09"

    def edit(cfg):
        cfg["inputs"].update(
            cells=str(dataset / "input.h5ad"),
            vocabulary=str(W / "runs/v8-dev-aifi-dev-002/plan/vocabulary.tsv"),
            tools={m: str(CAND / f"c005_aifi_baseline/preds/predictions_{m}.tsv") for m in ("popv", "singler", "celltypist")},
        )
        cfg["inputs"].pop("markers")
        cfg.pop("markers_table")
        cfg["votes"]["path"] = str(CAND / "c009_aifi_votes/cell_votes.parquet")
        cfg["consensus"].update(
            rule="h4", snap_to_vocabulary=False, gate=0.40, tool_order=["popv", "singler", "celltypist"],
            values_status="frozen c012 values, not approved defaults",
        )  # fmt: skip
        cfg["scoring"].update(truth=str(dataset / "scoring/truth.tsv"), dataset="aifi_dev")

    cfg = make_cfg(tmp_path, edit)
    out = tmp_path / "out"
    assert run("run", "--input", cfg, "--output", out, "--score") == 0
    got = json.loads((out / "score/aifi_dev.consensus_strict.json").read_text())
    row = frozen_row(
        CAND / "c012_harmonizer/frozen.tsv", run="aifi3", variant="H4cell", w_mode="equal", rule="R1"
    )
    assert got["coverage"] == pytest.approx(row["coverage"], abs=1e-6)
    assert got["descendant_accuracy"] == pytest.approx(row["desc_accuracy"], abs=1e-6)
