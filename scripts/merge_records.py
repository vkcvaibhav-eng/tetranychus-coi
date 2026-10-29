"""Module 2: merge GenBank and BOLD COI records into one table, remove
duplicates, and score how complete each record's metadata is.

Reads:  data/ncbi_*.gb, data/bold_*.tsv, data/bold_*.jsonl
Writes: results/merged_records.csv, results/merged_sequences.fasta,
        results/merge_summary.txt

Usage: python scripts/merge_records.py [data_dir] [results_dir]
"""
import csv
import glob
import json
import os
import re
import sys
from collections import Counter

from Bio import SeqIO

FIELDS = ["record_id", "source", "genbank_accession", "bold_processid", "bold_bin",
          "species_label", "country", "latitude", "longitude", "host",
          "collection_date", "voucher", "identified_by", "institution",
          "marker_label", "publication", "journal", "metadata_score",
          "metadata_missing", "seq_length", "sequence"]

# Seven items a well-documented barcode record should carry.
METADATA_ITEMS = {
    "country": lambda r: bool(r["country"]),
    "coordinates": lambda r: r["latitude"] != "" and r["longitude"] != "",
    "host": lambda r: bool(r["host"]),
    "collection_date": lambda r: bool(r["collection_date"]),
    "voucher": lambda r: bool(r["voucher"]),
    "identifier": lambda r: bool(r["identified_by"]),
    "publication": lambda r: bool(r["publication"]),
}

LATLON = re.compile(r"([\d.]+)\s*([NS])\s+([\d.]+)\s*([EW])")
NUM = re.compile(r"-?\d+(?:\.\d+)?")


def clean_seq(s):
    return re.sub(r"[^A-Za-z]", "", str(s or "")).upper()


def strip_version(acc):
    return acc.split(".")[0].strip()


def empty_record():
    return {k: "" for k in FIELDS}


def first(q, *keys):
    for k in keys:
        if q.get(k):
            return q[k][0]
    return ""


def parse_genbank(path):
    for rec in SeqIO.parse(path, "genbank"):
        r = empty_record()
        src = next((f for f in rec.features if f.type == "source"), None)
        q = src.qualifiers if src else {}
        r["record_id"] = rec.id
        r["source"] = "GenBank"
        r["genbank_accession"] = strip_version(rec.id)
        r["species_label"] = first(q, "organism") or rec.annotations.get("organism", "")
        r["country"] = first(q, "geo_loc_name", "country")
        m = LATLON.match(first(q, "lat_lon"))
        if m:
            r["latitude"] = float(m.group(1)) * (1 if m.group(2) == "N" else -1)
            r["longitude"] = float(m.group(3)) * (1 if m.group(4) == "E" else -1)
        r["host"] = first(q, "host")
        r["collection_date"] = first(q, "collection_date")
        r["voucher"] = first(q, "specimen_voucher", "bio_material")
        r["identified_by"] = first(q, "identified_by")
        genes = [f.qualifiers.get("gene", f.qualifiers.get("product", [""]))[0]
                 for f in rec.features if f.type in ("gene", "CDS")]
        r["marker_label"] = ";".join(sorted(set(g for g in genes if g)))
        refs = rec.annotations.get("references", [])
        published = [ref for ref in refs if ref.title and ref.title != "Direct Submission"
                     and "Unpublished" not in (ref.journal or "")]
        if published:
            r["publication"] = published[0].title
            r["journal"] = published[0].journal
        r["sequence"] = clean_seq(rec.seq)
        yield r


def read_bold_rows(data_dir):
    rows = []
    for path in glob.glob(os.path.join(data_dir, "bold_*.tsv")):
        with open(path, newline="") as f:
            rows += list(csv.DictReader(f, delimiter="\t"))
    for path in glob.glob(os.path.join(data_dir, "bold_*.jsonl")):
        with open(path) as f:
            text = f.read().strip()
        rows += json.loads(text) if text.startswith("[") else \
            [json.loads(l) for l in text.splitlines() if l.strip()]
    return rows


def parse_bold_row(row):
    seq = clean_seq(row.get("nuc"))
    marker = row.get("marker_code") or ""
    if not seq or not marker.upper().startswith("COI"):
        return None
    r = empty_record()
    r["record_id"] = row.get("processid", "")
    r["source"] = "BOLD"
    r["bold_processid"] = row.get("processid", "")
    r["bold_bin"] = row.get("bin_uri", "") or ""
    accs = [strip_version(a) for a in re.split(r"[,;\s]+", str(row.get("insdc_acs") or "")) if a]
    r["genbank_accession"] = accs[0] if accs else ""
    r["species_label"] = row.get("species") or row.get("identification") or ""
    r["country"] = row.get("country/ocean") or ""
    coord = row.get("coord")
    if isinstance(coord, list) and len(coord) == 2:
        r["latitude"], r["longitude"] = float(coord[0]), float(coord[1])
    else:
        nums = NUM.findall(str(coord or ""))
        if len(nums) == 2:
            r["latitude"], r["longitude"] = float(nums[0]), float(nums[1])
    r["host"] = row.get("associated_taxa") or ""
    r["collection_date"] = row.get("collection_date_start") or row.get("collection_date") or ""
    r["voucher"] = row.get("museumid") or row.get("sampleid") or ""
    r["identified_by"] = row.get("identified_by") or ""
    r["institution"] = row.get("inst") or ""
    r["marker_label"] = marker
    r["sequence"] = seq
    return r


def merge_into(base, extra):
    """Fill empty fields of a GenBank record from its BOLD copy."""
    for k in ("bold_processid", "bold_bin", "latitude", "longitude", "host",
              "collection_date", "voucher", "identified_by", "institution", "country"):
        if base[k] == "" and extra[k] != "":
            base[k] = extra[k]
    base["source"] = "GenBank+BOLD"


def score(r):
    missing = [name for name, test in METADATA_ITEMS.items() if not test(r)]
    r["metadata_score"] = len(METADATA_ITEMS) - len(missing)
    r["metadata_missing"] = ";".join(missing)


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "data"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "results"
    os.makedirs(out_dir, exist_ok=True)

    log = Counter()
    by_acc = {}
    for path in sorted(glob.glob(os.path.join(data_dir, "ncbi_*.gb"))):
        for r in parse_genbank(path):
            log["genbank_read"] += 1
            if r["genbank_accession"] in by_acc:
                log["genbank_repeat_removed"] += 1
                continue
            by_acc[r["genbank_accession"]] = r

    bold_only, seen_pid = [], set()
    for row in read_bold_rows(data_dir):
        r = parse_bold_row(row)
        if r is None:
            log["bold_skipped_no_coi"] += 1
            continue
        log["bold_read"] += 1
        if r["bold_processid"] in seen_pid:
            log["bold_repeat_removed"] += 1
            continue
        seen_pid.add(r["bold_processid"])
        if r["genbank_accession"] and r["genbank_accession"] in by_acc:
            merge_into(by_acc[r["genbank_accession"]], r)
            log["bold_matched_to_genbank"] += 1
        else:
            bold_only.append(r)

    records = list(by_acc.values()) + bold_only
    for r in records:
        r["seq_length"] = len(r["sequence"])
        score(r)

    with open(os.path.join(out_dir, "merged_records.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(records)
    with open(os.path.join(out_dir, "merged_sequences.fasta"), "w") as f:
        for r in records:
            f.write(f">{r['record_id']}\n{r['sequence']}\n")

    src = Counter(r["source"] for r in records)
    scores = Counter(r["metadata_score"] for r in records)
    lines = [f"Total unique records: {len(records)}"]
    lines += [f"  {k}: {v}" for k, v in sorted(src.items())]
    lines += ["Processing:"] + [f"  {k}: {v}" for k, v in sorted(log.items())]
    lines += ["Metadata score (0-7, records):"]
    lines += [f"  {s}: {scores[s]}" for s in range(len(METADATA_ITEMS) + 1)]
    missing = Counter(m for r in records for m in r["metadata_missing"].split(";") if m)
    lines += ["Records missing each item:"]
    lines += [f"  {k}: {v} ({100 * v / len(records):.1f}%)" for k, v in missing.most_common()] \
        if records else []
    summary = "\n".join(lines)
    with open(os.path.join(out_dir, "merge_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
