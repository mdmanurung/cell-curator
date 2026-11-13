"""Register the `slow` marker used by the annotation-consensus tests."""


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: full AIDA/AIFI regression against the frozen benchmark numbers"
    )
