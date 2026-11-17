from pathlib import Path

import pytest
from consensus_fixture import write_consensus_run
from test_annotation_pipeline import annotation_fixture

from cell_curator.annotation import collect_level
from cell_curator.consensus_import import import_consensus_level
from cell_curator.guided import _next_command, load_state


def test_corpus_progress_moves_from_compiled_execution_through_review_to_critics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, expected = annotation_fixture(tmp_path, monkeypatch)
    state = load_state(root)
    assert "annotation collect-level" in _next_command(config_path, root, state)
    collect_level(config_path, "L1")
    assert "waiting for the consensus run" in _next_command(config_path, root, state)
    write_consensus_run(config_path, root, {c: v[0] for c, v in expected.items()})
    assert "annotation import-consensus" in _next_command(config_path, root, state)
    import_consensus_level(config_path, "L1")
    assert "annotation finalize" in _next_command(config_path, root, state)
    manifest = root / "annotation" / "annotation.manifest.json"
    manifest.parent.mkdir(exist_ok=True)
    manifest.write_text("{}")
    command = _next_command(config_path, root, state)
    assert "snakemake" in command and "annotation_review" in command
    report = root / "review" / "report.manifest.json"
    report.parent.mkdir(exist_ok=True)
    report.write_text("{}")
    assert "pipeline packet" in _next_command(config_path, root, state)
    packet = root / "review_packet" / "review_packet.manifest.json"
    packet.parent.mkdir()
    packet.write_text("{}")
    assert "critic import-findings" in _next_command(config_path, root, state)
