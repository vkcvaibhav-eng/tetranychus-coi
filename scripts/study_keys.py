"""Submission batches ("studies") for grouping records.

Records deposited together by one group are not independent: they share
specimens, laboratories, primers and naming decisions. The checks in
robustness.py, its2_analysis.py and trends.py therefore group records by
submission batch:

  GenBank records   first author (or consortium) of the "Direct Submission"
                    reference and the year of submission (for example
                    "li,g.-q.|2009")
  BOLD-only records the BOLD project code, i.e. the letters at the start of
                    the process identifier (for example "CNSID")

A record that cannot be placed forms its own batch.
"""
import glob
import re

from Bio import SeqIO

SUBMITTED = re.compile(r"Submitted \(\d{2}-[A-Z]{3}-(\d{4})\)")


def first_author(authors):
    """'Li,G.-Q., Xue,X.-F. and Hong,X.-Y.' -> 'li,g.-q.'"""
    authors = (authors or "").strip()
    if not authors:
        return ""
    head = re.split(r"\s+and\s+", authors)[0]
    parts = [p.strip() for p in head.split(",")]
    return (",".join(parts[:2]) if len(parts) > 1 else parts[0]).lower()


def submission_key(rec):
    """First author (or consortium) of the direct submission, with its year."""
    refs = rec.annotations.get("references", [])
    ds = next((r for r in refs if (r.title or "").strip() == "Direct Submission"), None)
    m = SUBMITTED.search(ds.journal or "") if ds else None
    year = m.group(1) if m else ""
    for ref in ([ds] if ds else []) + refs:
        who = first_author(ref.authors) or (getattr(ref, "consrtm", "") or "").strip().lower()
        if who:
            return f"{who}|{year}"
    return ""  # no author or consortium: the record is not grouped with others


def genbank_keys(pattern):
    """{accession without version: batch key} for every GenBank file matching pattern."""
    keys = {}
    for path in sorted(glob.glob(pattern)):
        for rec in SeqIO.parse(path, "genbank"):
            k = submission_key(rec)
            if k:
                keys.setdefault(rec.id.split(".")[0], k)
    return keys


def bold_project(processid):
    m = re.match(r"[A-Za-z]+", processid or "")
    return f"bold:{m.group(0).upper()}" if m else ""


def record_keys(merged_rows, gb_keys):
    """Batch key for each merged COI record."""
    out = {}
    for r in merged_rows:
        k = gb_keys.get(r.get("genbank_accession", "")) or bold_project(r.get("bold_processid", ""))
        out[r["record_id"]] = k or f"single:{r['record_id']}"
    return out
