"""Register the `slow` marker and skip it unless the full run is asked for."""

import os


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: full AIDA/AIFI regression against the frozen benchmark numbers"
    )
    # Default is the fast run (about 6 min). Set ANNOTATION_CONSENSUS_FULL=1, or pass -m yourself,
    # to include the slow data regression.
    if not config.option.markexpr and os.environ.get("ANNOTATION_CONSENSUS_FULL") != "1":
        config.option.markexpr = "not slow"
