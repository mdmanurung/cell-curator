from pathlib import Path

import pytest
from test_corpus_writeback import _completed

from cell_curator.pipeline import build_review_packet


def test_frozen_packet_keeps_corpus_definitions_and_state_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    config, root, _ = _completed(tmp_path, monkeypatch)
    provider = config.knowledge.providers[0]
    original = Path(provider.path).read_bytes()
    snapshot = root / "plan" / f"corpus-{provider.sha256}.json"
    assert snapshot.read_bytes() == original
    packet = build_review_packet(tmp_path / "config.yaml")
    assert (packet / "plan" / snapshot.name).read_bytes() == original
    assert (packet / "annotation" / "cell_labels.tsv").is_file()
    assert (packet / "annotation" / "levels" / "L3" / "states.tsv").is_file()
    # A frozen review packet is checked against its own inventory, not live sources.
    Path(provider.path).write_bytes(original + b"\n")
    assert build_review_packet(tmp_path / "config.yaml") == packet
    assert (packet / "plan" / snapshot.name).read_bytes() == original
