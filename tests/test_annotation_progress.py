from pathlib import Path

import pytest
from test_annotation_pipeline import annotation_fixture

from cell_curator.guided import _next_command, load_state


def test_corpus_progress_moves_from_compiled_execution_through_review_to_critics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path, root, _ = annotation_fixture(tmp_path, monkeypatch)
    state = load_state(root)
    assert "run execute" in _next_command(config_path, root, state)
    manifest = root / "annotation" / "annotation.manifest.json"
    manifest.parent.mkdir()
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
