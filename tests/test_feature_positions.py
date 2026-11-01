from cell_curator.contracts import feature_positions


def test_marker_genes_match_case_insensitively_but_never_ambiguously() -> None:
    positions = feature_positions(["Cd3e", "MS4A1", "Foo", "FOO"])
    assert positions["CD3E"] == 0 and "CD3E" in positions and positions.get("cd3e") == 0
    assert positions["MS4A1"] == 1
    assert "foo" not in positions and positions.get("fOO") is None
    assert positions["FOO"] == 3
