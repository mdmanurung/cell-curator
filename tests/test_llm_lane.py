from __future__ import annotations

import json
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from test_end_to_end import _config, _confirm, _record_critics

from cell_curator.cli import build_parser, dispatch
from cell_curator.config import load_config
from cell_curator.contracts import ContractError
from cell_curator.llm import (
    PREFIX,
    Proposal,
    ProposalSet,
    apply_proposals,
    build_packet,
    lane_bottom_up,
    packet_from_run,
    packet_sha256,
    propose,
)
from cell_curator.pipeline import run_annotation, run_root
from cell_curator.writeback import write_back


def _tied_run(tmp_path: Path) -> Path:
    """Cluster 0 is half Type A and half Type B cells, so it abstains on a tie."""

    rng = np.random.default_rng(5)
    counts = np.zeros((60, 6), dtype=np.float32)
    counts[:15, :2] = rng.poisson(8, size=(15, 2)) + 1
    counts[15:30, 2:4] = rng.poisson(8, size=(15, 2)) + 1
    counts[30:, 4:] = rng.poisson(8, size=(30, 2)) + 1
    genes = ["A1", "A2", "B1", "B2", "N1", "N2"]
    adata = ad.AnnData(
        X=np.log1p(counts),
        obs=pd.DataFrame(
            {
                "canonical_cluster": ["0"] * 30 + ["1"] * 30,
                "lineage": ["immune"] * 60,
                "donor": ["donor-secret-1"] * 60,
            },
            index=[f"barcode-{index:03d}" for index in range(60)],
        ),
        var=pd.DataFrame({"gene_symbol": genes}, index=genes),
    )
    adata.layers["counts"] = counts
    adata.layers["logcounts"] = np.log1p(counts)
    adata.obsm["X_pca_harmony"] = np.vstack(
        [rng.normal(-3, 0.2, size=(30, 3)), rng.normal(3, 0.2, size=(30, 3))]
    ).astype(np.float32)
    source = tmp_path / "source.h5ad"
    adata.write_h5ad(source)
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            {
                "Type A": {"pos": ["A1", "A2"], "neg": ["N1", "N2"]},
                "Type B": {"pos": ["B1", "B2"], "neg": ["N1", "N2"]},
                "Type N": {"pos": ["N1", "N2"], "neg": ["A1", "B1"]},
            }
        )
    )
    return _config(tmp_path, source, markers)


def _answer(packet: dict, label: str, genes: tuple[str, ...], model: str = "m") -> Proposal:
    return Proposal(cluster_id="0", label=label, cited_genes=genes, model=model)


@pytest.fixture
def first_run(tmp_path: Path):
    config_path = _tied_run(tmp_path)
    run_annotation(config_path)
    config = load_config(config_path)
    root = run_root(config)
    proposals = pd.read_csv(
        root / "evidence/summary/annotation_proposals.tsv", sep="\t", keep_default_na=False
    ).set_index("cluster_id")
    assert proposals.loc[0, "cell_type"] == "Unknown"
    return config_path, config, root, packet_from_run(config, root)


def test_the_packet_carries_cluster_statistics_and_nothing_per_cell(first_run) -> None:
    _, config, root, packet = first_run
    raw = json.dumps(packet)
    assert [entry["cluster_id"] for entry in packet["clusters"]] == ["0"]
    entry = packet["clusters"][0]
    assert {item["label"] for item in entry["candidates"]} == {"Type A", "Type B"}
    assert "pseudobulk" in entry
    assert "barcode-" not in raw and "donor-secret" not in raw

    # Below the floor a cluster sends no pseudobulk means.
    small = config.model_copy(
        update={
            "evidence": config.evidence.model_copy(
                update={"llm": config.evidence.llm.model_copy(update={"min_cells": 31})}
            )
        }
    )
    assert "pseudobulk" not in packet_from_run(small, root)["clusters"][0]


def test_every_gate_refuses_its_own_failure(first_run) -> None:
    _, config, root, packet = first_run
    proposals = pd.read_csv(
        root / "evidence/summary/annotation_proposals.tsv", sep="\t", keep_default_na=False
    )
    aggregate = pd.read_csv(
        root / "evidence/summary/multimodal_label_scores.tsv", sep="\t", keep_default_na=False
    )
    from cell_curator.models import AnnotationDecision

    decisions = [
        AnnotationDecision(
            run_id=row["run_id"],
            scope=row["scope"],
            cluster_id=str(row["cluster_id"]),
            cell_type=row["cell_type"],
            confidence=row["confidence"],
            margin=row["margin"],
            normalized_entropy=row["normalized_entropy"],
            abstained=bool(row["abstained"]),
            rationale=row["rationale"],
            required_capabilities_complete=bool(row["required_capabilities_complete"]),
        )
        for row in proposals.to_dict("records")
    ]
    digest = packet_sha256(packet)

    def outcome(*answers: Proposal) -> dict:
        _, outcomes = apply_proposals(
            config,
            packet,
            ProposalSet(payload_sha256=digest, proposals=list(answers)),
            decisions,
            aggregate,
        )
        return outcomes[0]

    assert outcome(_answer(packet, "Type A", ("A1",)))["status"] == "accepted"
    assert "not one of the near-tied" in outcome(_answer(packet, "Type N", ("N1",)))["reason"]
    assert "not positive markers" in outcome(_answer(packet, "Type A", ("N1",)))["reason"]
    assert "not positive markers" in outcome(_answer(packet, "Type A", ()))["reason"]
    assert (
        "disagree"
        in outcome(
            _answer(packet, "Type A", ("A1",), "m1"), _answer(packet, "Type B", ("B1",), "m2")
        )["reason"]
    )
    with pytest.raises(ContractError, match="different packet"):
        apply_proposals(
            config,
            packet,
            ProposalSet(payload_sha256="0" * 64, proposals=[]),
            decisions,
            aggregate,
        )
    # A gate-passing candidate whose own negatives are violated is refused.
    violated = aggregate.copy()
    violated.loc[violated["label"].eq("Type A"), "negative_violation"] = 1.0
    _, outcomes = apply_proposals(
        config,
        packet,
        ProposalSet(payload_sha256=digest, proposals=[_answer(packet, "Type A", ("A1",))]),
        decisions,
        violated,
    )
    assert "own coverage or negative-marker gates" in outcomes[0]["reason"]


def test_a_tie_break_is_bound_at_plan_and_reaches_the_written_h5ad(
    first_run, tmp_path: Path
) -> None:
    config_path, _, _, packet = first_run
    answers = tmp_path / "proposals.json"
    answers.write_text(
        ProposalSet(
            payload_sha256=packet_sha256(packet),
            proposals=[_answer(packet, "Type A", ("A1", "A2"), "local-model")],
        ).model_dump_json()
    )
    value = yaml.safe_load(config_path.read_text())
    value["run"]["run_id"] = "with-llm"
    value["evidence"]["llm"] = {"proposals": str(answers)}
    second = tmp_path / "second.yaml"
    second.write_text(yaml.safe_dump(value, sort_keys=False))

    run_annotation(second)
    config = load_config(second)
    root = run_root(config)
    assert json.loads((root / "plan.json").read_text())["llm_proposals_sha256"]
    proposals = pd.read_csv(
        root / "evidence/summary/annotation_proposals.tsv", sep="\t", keep_default_na=False
    ).set_index("cluster_id")
    assert proposals.loc[0, "cell_type"] == "Type A"
    assert proposals.loc[0, "rationale"].startswith(PREFIX)
    # Below the confidence gate by construction: the model broke a measured tie.
    assert proposals.loc[0, "confidence"] < config.evidence.min_confidence

    _record_critics(config, root, tmp_path)
    output = tmp_path / "annotated.h5ad"
    write_back(config, root, output_path=output, confirm=_confirm)
    written = ad.read_h5ad(output).obs["cell_curator_type"].astype(str)
    assert set(written.iloc[:30]) == {"Type A"}

    # Editing the proposals after plan bound them is refused.
    answers.write_text(answers.read_text().replace("local-model", "other-model"))
    with pytest.raises(ContractError, match="changed after plan bound it"):
        from cell_curator.evidence import consolidate_lanes

        consolidate_lanes(config=config, run_root=root)


class _FakeModel(BaseHTTPRequestHandler):
    seen: list[dict] = []

    def do_POST(self) -> None:  # noqa: N802
        body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
        _FakeModel.seen.append({"auth": self.headers.get("Authorization"), "body": body})
        content = json.dumps({"label": "Type A", "cited_genes": ["A1"], "rationale": "A1 up"})
        payload = json.dumps({"choices": [{"message": {"content": content}}]}).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(payload)

    def log_message(self, *_: object) -> None:
        pass


def test_propose_asks_a_local_server_and_offline_refuses(
    first_run, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, config, root, packet = first_run
    server = HTTPServer(("127.0.0.1", 0), _FakeModel)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{server.server_port}/v1"
    monkeypatch.setenv("LOCAL_KEY", "sk-secret-value")
    institution = config.model_copy(
        update={"knowledge": config.knowledge.model_copy(update={"network": "institution"})}
    )
    try:
        answer = propose(
            institution,
            packet,
            base_url=url,
            model="llama",
            out_dir=tmp_path / "x",
            api_key_env="LOCAL_KEY",
        )
    finally:
        server.shutdown()
    assert answer.payload_sha256 == packet_sha256(packet)
    assert [(p.label, p.cited_genes) for p in answer.proposals] == [("Type A", ("A1",))]
    sent = _FakeModel.seen[-1]
    assert sent["auth"] == "Bearer sk-secret-value"
    assert sent["body"]["temperature"] == 0
    stored = b"".join(path.read_bytes() for path in (tmp_path / "x").iterdir())
    assert b"sk-secret-value" not in stored and b"barcode-" not in stored

    with pytest.raises(ContractError, match="offline"):
        propose(config, packet, base_url=url, model="llama", out_dir=tmp_path / "y")


def test_the_host_agent_packet_is_refused_unless_online_or_declared_local(first_run) -> None:
    config_path, _, root, _ = first_run
    args = build_parser().parse_args(
        ["llm", "packet", "--config", str(config_path), "--output", str(root / "p.json")]
    )
    with pytest.raises(ContractError, match="leaves the institution"):
        dispatch(args)
    value = yaml.safe_load(config_path.read_text())
    value["evidence"]["llm"] = {"host_agent_local": True}
    config_path.write_text(yaml.safe_dump(value, sort_keys=False))
    assert dispatch(args)["clusters"] == 1


def test_lane_bottom_up_reads_every_configured_lane(first_run) -> None:
    _, config, root, packet = first_run
    assert set(lane_bottom_up(config, root)["cluster_id"].astype(str)) >= {"0"}
    assert build_packet(config, [], pd.DataFrame(), pd.DataFrame(), {})["clusters"] == []
