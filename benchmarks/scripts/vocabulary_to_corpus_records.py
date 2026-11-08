"""Turn the protocol v5/v6 marker vocabulary into corpus records, without changing a gene.

Purpose: protocol v7a runs the new engine on the same 30 Azimuth PBMC programs v6 used,
so the engine is the only change.
Run:     python benchmarks/scripts/vocabulary_to_corpus_records.py OUT_RECORDS.json
Inputs:  benchmarks/data/aida-pbmc-protocol-v5/pbmc_markers.json  [build_vocabulary.py]
Outputs: OUT_RECORDS.json  records for `cell-curator corpus import --source-kind generic`
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VOCABULARY = ROOT / "benchmarks/data/aida-pbmc-protocol-v5/pbmc_markers.json"
CONTEXT = {"species": "human", "tissue": "blood", "stage": "adult", "condition": "healthy"}


def main(out: Path) -> None:
    programs = json.loads(VOCABULARY.read_text())
    labels = [
        {
            "source_record_id": f"azimuth-l2-{index:02d}",
            "name": name,
            "source_label": item["azimuth_label"],
            "description": (
                f"Azimuth PBMC celltype.l2 type {item['azimuth_label']}: ten published markers; "
                "negatives are the nearest sibling's positives (protocol v5 rule)."
            ),
            "positive": item["pos"],
            "negative": item["negative"],
            "cl_id": item["cl_id"],
            "ontology_relation": "exact",  # as v6 used it; the ontology check runs at export
            "canonical_level": "L1",  # v6 ran one flat level
            "source_level": "celltype.l2",
        }
        for index, (name, item) in enumerate(programs.items())
    ]
    out.write_text(json.dumps({"context": CONTEXT, "labels": labels}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
