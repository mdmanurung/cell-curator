# Contributing

Report bugs and ask questions in [GitHub Issues](https://github.com/mdmanurung/cell-curator/issues).
Include the command, the config file, and the full error.

To change code, open a pull request against `main`. It must pass:

```bash
uv sync --frozen --extra dev --extra workflow
uv run ruff check .
uv run mypy
uv run pytest -q
```

A change to a benchmark result follows the frozen protocol in
[benchmarks/README.md](benchmarks/README.md): declare a new protocol version before
running, never rewrite a published result.
