"""A missing identifier must never become a short string that reads as data.

`str(None)` yields "None" and `str(float("nan"))` yields "nan". Both are
truthy, so both survive an `if value` guard. A Cell Ontology identifier is
carried through a run as a string in several places, and each of those places
was a place a missing accession could be written into a durable, hashed artifact
as though the run had one.

The observed consequence: a run configured without an ontology produced a label
table whose `cl_id` column held the literal text "None", and a child evidence
summary that recorded "None" as a consensus accession -- which also meant the
`conflicting_ontology_accessions` gate could never fire for an absent
accession, because every lane agreed on the same fabricated value. A consumer
resolving it later failed with "ontology term is unknown: None", naming no
stage.

This file pins the normalizer itself. The places that use it are pinned where
they are reachable: `test_prediction_export.py` asserts a real run's artifacts
carry no stringified null, and `test_scientific_contracts.py` covers the
hierarchy and decision paths.
"""

from __future__ import annotations

from cell_curator.contracts import identifier_text


def test_the_normalizer_maps_every_shape_of_absence_to_absent() -> None:
    """A missing identifier has four spellings and all four must read as absent.

    `str(value)` turns `None` into "None" and `NaN` into "nan", and both strings
    are truthy, so they survive an `if value` guard and are written into durable
    artifacts as identifiers the run never had. The consequence was a label table
    reporting `cl_id` of "None" for a run configured without an ontology, and a
    consumer resolving it failing on a term it had never heard of with nothing
    naming the stage that produced it. An accession is never literally "none" or
    "nan", so this mapping is unambiguous.
    """

    for absent in (None, float("nan"), "", "   ", "None", "none", "NULL", "nan"):
        assert identifier_text(absent) == "", f"{absent!r} must read as absent"
    for present in ("CL:0000084", " PCL:0000123 ", "CL:0000625"):
        assert identifier_text(present) == present.strip()
