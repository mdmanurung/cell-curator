# Purpose: turn one tusoskill candidate JSON into a candidate-only corpus proposal and, once the
#          user has reviewed it, into per-dev-set configs and the commands to run and score them.
#          The script never approves anything: it stops before `corpus review` (AGENTS.md
#          gate 3) and before `assumptions approve` (Gate A). A candidate corpus is a
#          development artifact and never replaces benchmarks/data/v8-corpus/reviewed_v3.json.
# Run:     .venv/bin/python benchmarks/scripts/tuso_candidate_v8.py CAND.json
#              engine-only candidate (only "name" and "engine"): no corpus step; it goes straight to
#              stage 2 on the frozen corpus.
#              corpus candidate, stage 1: learn, import (candidate-only), propose; prints the
#              review command and stops.
#          .venv/bin/python benchmarks/scripts/tuso_candidate_v8.py CAND.json --reviewed R.json
#              corpus candidate, stage 2: check isolation, write configs with fresh run_ids,
#              print the commands in a pre-gate and a post-gate file per dev set.
#          --learn-dir DIR uses an existing learn_corpus_v8.py output instead of submitting one;
#          --root DIR moves CANDIDATES_ROOT (for tests).
# Inputs:  CAND.json   {"name": "...", and any key of learn_corpus_v8.SELECTION_KEYS,
#                       and optional "engine": {min_confidence, min_margin, ...}}
#          LEARN_PREVIOUS/ranked.tsv.gz  the frozen Wilcoxon table, reused (selection-only)
#          BASE_CORPUS, FROZEN_CORPUS, the two dev-set config templates
#          R.json      stage 2 only: the candidate proposal after the user's `corpus review`
# Outputs: CANDIDATES_ROOT/NAME/candidate.json        the candidate as run
#          CANDIDATES_ROOT/NAME/learn/                learn_corpus_v8.py outputs
#          CANDIDATES_ROOT/NAME/corpus/imported.json  imported, not reviewed, not callable
#          CANDIDATES_ROOT/NAME/corpus/proposal.json  for the user's review
#          CANDIDATES_ROOT/NAME/CANDIDATE_ONLY.txt    the marker
#          CANDIDATES_ROOT/NAME/configs/{dataset}.yaml  stage 2
#          CANDIDATES_ROOT/NAME/commands/{dataset}_pre_gate.sh   scene, plan, assumptions list
#          CANDIDATES_ROOT/NAME/commands/{dataset}_post_gate.sh  execute, validate, export, score

import argparse
import datetime
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

from cell_curator.contracts import sha256
from cell_curator.corpus import check_test_isolation, load_corpus, propose_corpus

# Parameters -----------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
WORK = Path("/exports/para-lipg-hpc/mdmanurung/cc-v8-work")
CANDIDATES_ROOT = WORK / "candidates"
LEARN_PREVIOUS = WORK / "learn-2026-11-09c"  # ranked.tsv.gz of the frozen sample
BASE_CORPUS = WORK / "corpus/base.json"  # the empty snapshot the frozen corpus grew from
FROZEN_CORPUS = ROOT / "benchmarks/data/v8-corpus/reviewed_v3.json"  # which 34 labels are in
SOURCE_URL = "https://cellxgene.cziscience.com/collections/77f9d7e9-5675-49c3-abed-ce02f39eef1b"
# dataset -> template config, prep dir pattern (cluster count, truth), method name for scoring
DEV_SETS = {
    "aifi_dev": ("benchmarks/data/v8-dev/aifi_dev.yaml", "prep-aifi_dev-*"),
    "aida_split_a": ("benchmarks/data/v8-dev/aida_split_a.yaml", "prep-aida_split_a-2026-*"),
}
# the 3'-only note the user chose for the frozen records (amendment 6, revision v8-r3)
ASSAY_NOTE = (
    "Learned on 10x 3' v3 only; applied to any RNA assay by the user's choice (2026-11-09)."
)
# engine settings a candidate may change, all under `evidence:` in the config
ENGINE_KEYS = (
    "adjusted_pvalue_max",
    "max_negative_violation",
    "min_confidence",
    "min_margin",
    "max_normalized_entropy",
    "min_positive_coverage",
)
DATASET_ID = "e522d2cd-7927-4e59-a4ed-064009569279"  # the atlas file the donors come from
ENGINE_ONLY = {"name", "note", "engine"}  # a candidate with only these keys builds no corpus
CLI = ".venv/bin/cell-curator"
PYTHON = ".venv/bin/python"


def cli(*args):
    subprocess.run(
        [sys.executable, "-m", "cell_curator", *map(str, args)], check=True, cwd=ROOT
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--reviewed", type=Path)
    parser.add_argument("--learn-dir", type=Path)
    parser.add_argument("--root", type=Path, default=CANDIDATES_ROOT)
    args = parser.parse_args()
    candidate = json.loads(args.candidate.read_text())
    name = candidate["name"]
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]{1,40}", name):
        raise ValueError("candidate name is 2 to 41 characters of a-z, 0-9 and -")
    out = args.root / name
    if set(candidate) <= ENGINE_ONLY:
        finish(candidate, out, FROZEN_CORPUS, frozen=True, candidate_path=args.candidate)
    elif args.reviewed:
        finish(candidate, out, args.reviewed, frozen=False)
    else:
        start(candidate, args.candidate, out, args.learn_dir)


# Stage 1 ---------------------------------------------------------------------------
def start(candidate, candidate_path, out, learn_dir):
    name = candidate["name"]
    out.mkdir(parents=True, exist_ok=False)  # a candidate dir is never overwritten
    (out / "corpus").mkdir()
    (out / "candidate.json").write_text(candidate_path.read_text())
    (out / "CANDIDATE_ONLY.txt").write_text(
        f"Candidate {name}: development artifact. Never reviewed into the frozen corpus.\n"
    )

    # Learn: selection only, the Wilcoxon table is reused -----------------------------------
    if learn_dir is None:
        learn_dir = out / "learn"
        subprocess.run(
            [
                "sbatch", "--wait", "--mem=96G", "--time=03:00:00",
                "benchmarks/scripts/learn_corpus_v8.sbatch", learn_dir,
                "--candidate", out / "candidate.json", "--reuse-ranked", LEARN_PREVIOUS,
            ],
            check=True,
            cwd=ROOT,
        )  # fmt: skip

    # Import: status imported, so the package cannot call any record ------------------------
    imported = out / "corpus/imported.json"
    cli(
        "corpus", "import", "--input", BASE_CORPUS, "--records", learn_dir / "corpus_records.json",
        "--output", imported, "--source-kind", "generic", "--source-id", f"hiha-v8-cand-{name}",
        "--name", f"CANDIDATE {name}: Human Immune Health Atlas markers, development only",
        "--url", SOURCE_URL, "--source-version", f"cand-{name}",
        "--retrieved-at", datetime.date.today().isoformat(),
    )  # fmt: skip

    # Propose: the frozen corpus's 34 accepted labels, same scope edits as revision v8-r3 ---
    frozen = {
        label.cl_id for label in load_corpus(FROZEN_CORPUS).labels if label.status == "reviewed"
    }
    donors = (learn_dir / "corpus_donors.txt").read_text().split()  # derived_from, so the
    snapshot = load_corpus(imported)  # isolation check has donors to compare
    changes = [
        {
            **label.model_dump(mode="json"),
            "assays": [],
            "limitations": [*label.limitations, ASSAY_NOTE],
            "derived_from": [{"dataset_id": DATASET_ID, "donors": donors}],
        }
        for label in snapshot.labels
        if label.cl_id in frozen
    ]
    (out / "corpus/changes.json").write_text(json.dumps({"labels": changes}, indent=2))
    propose_corpus(
        imported, out / "corpus/proposal.json", labels=changes,
        reason=f"Development candidate {name}; not for the frozen corpus", proposer="tuso_candidate_v8",
    )  # fmt: skip

    callable_now = len(snapshot.select({"species": "human", "tissue": "blood"}))
    print(
        f"\ncandidate {name}: {len(snapshot.labels)} labels imported, {callable_now} callable."
    )
    print("STOP. The package only calls records with status `reviewed` (corpus.py select());")
    print("there is no flag for development use. A run on this corpus would label nothing.")
    print("To go on, the user reviews the proposal in chat, then runs:\n")
    print(
        f"  {CLI} corpus review --input {imported} --proposal {out}/corpus/proposal.json \\\n"
        f"    --output {out}/corpus/reviewed.json --reviewer 'YOUR NAME' --decision accept \\\n"
        f"    --revision cand-{name} --rationale 'development candidate {name}, not frozen'\n"
    )
    print(f"Then: {PYTHON} benchmarks/scripts/tuso_candidate_v8.py {candidate_path} \\")
    print(f"        --reviewed {out}/corpus/reviewed.json")


# Stage 2 ---------------------------------------------------------------------------
def finish(candidate, out, corpus, frozen, candidate_path=None):
    name = candidate["name"]
    # a finished prep (cluster counts) must exist for every dev set before anything is written
    preps = {}
    for dataset, (_, prep_pattern) in DEV_SETS.items():
        found = sorted(WORK.glob(prep_pattern))
        if not found or not (found[-1] / "cluster_sizes.tsv").is_file():
            raise FileNotFoundError(f"{dataset}: no finished prep dir {WORK}/{prep_pattern}")
        preps[dataset] = found[-1]
    snapshot = load_corpus(corpus)
    version = "v8-r3" if frozen else f"cand-{name}"
    if snapshot.revision != version or any(
        s.version != f"cand-{name}" for s in snapshot.sources if not frozen
    ):
        raise ValueError(f"{corpus} is not this candidate's corpus (revision {version})")
    digest = sha256(corpus)
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    out.mkdir(parents=True, exist_ok=True)
    if candidate_path:
        (out / "candidate.json").write_text(candidate_path.read_text())
    (out / "configs").mkdir(exist_ok=False)
    (out / "commands").mkdir(exist_ok=False)
    note = (
        "the frozen reviewed corpus"
        if frozen
        else f"CANDIDATE-ONLY {name}: development, never the frozen corpus"
    )
    for dataset, (template, _) in DEV_SETS.items():
        text = (ROOT / template).read_text()
        isolation = yaml.safe_load(text)["annotation"]["isolation"]
        verdict = check_test_isolation(
            snapshot, isolation["donor_split"], test_split=isolation["test_split"],
            dataset_id=isolation["dataset_id"],
        )  # fmt: skip
        prep = preps[dataset]
        clusters = len((prep / "cluster_sizes.tsv").read_text().splitlines()) - 1
        run_id = f"v8-dev-{dataset.replace('_', '-')}-{name}-{stamp}"

        # fresh run_id, the corpus, the cluster count of the prepared input
        text = re.sub(r"(?m)^  run_id: .*$", f"  run_id: {run_id}", text)
        text = re.sub(
            r"(?s)(    - kind: corpus\n).*?(\n    - kind: ontology)",
            rf"\1      # {note}\n"
            rf"      path: {corpus}\n      version: {version}\n"
            rf'      sha256: "{digest}"\2',
            text,
        )
        text = re.sub(r"(?m)^(    expected_parent_count: )\d+.*$", rf"\g<1>{clusters}", text)
        for key, value in candidate.get("engine", {}).items():
            if key not in ENGINE_KEYS:
                raise ValueError(f"engine key {key} is not one of {ENGINE_KEYS}")
            text, hits = re.subn(rf"(?m)^  {key}: .*$", f"  {key}: {value}", text)
            if hits != 1:
                raise ValueError(f"engine key {key} found {hits} times in {template}")
        config = out / f"configs/{dataset}.yaml"
        config.write_text(text)

        truth = prep / "scoring/truth.tsv"
        export_truth = out / f"scoring_{dataset}.truth_for_export.tsv"
        predictions = out / f"predictions_{dataset}.tsv"
        pre_gate = [
            f"# {dataset}: isolation {verdict['status']}, run_id {run_id}",
            f"{CLI} run scene --config {config}",
            f"{CLI} run plan --config {config}",
            f"{CLI} assumptions list --config {config}",
        ]
        post_gate = [
            "# Gate A first: the user reads the assumptions and runs `assumptions approve`",
            f"# for run_id {run_id}. An agent must not.",
            f"{CLI} run execute --config {config}",
            f"{CLI} run validate --config {config}",
            "# scoring (dev set truth, read at scoring time only):",
            f"{PYTHON} -c \"import pandas as pd,sys; pd.read_csv(sys.argv[1],sep='\\t',dtype=str,"
            "keep_default_na=False).rename(columns={'cell_type':'true_label',"
            "'cell_type_ontology_term_id':'true_cl_id'}).to_csv(sys.argv[2],sep='\\t',index=False)\" "
            f"{truth} {export_truth}",
            f"{CLI} run export-predictions --config {config} --truth {export_truth} "
            f"--output {predictions}",
            f"{PYTHON} benchmarks/scripts/score_v8.py {predictions} {truth} {out}/scores "
            f"--method cell-curator-{name} --dataset {dataset}",
        ]
        for part, lines in (("pre_gate", pre_gate), ("post_gate", post_gate)):
            (out / f"commands/{dataset}_{part}.sh").write_text("\n".join(lines) + "\n")
        print("\n".join([*pre_gate, "# --- STOP: Gate A (user) ---", *post_gate, ""]))


if __name__ == "__main__":
    main()
