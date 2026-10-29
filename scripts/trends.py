"""Module 9: how data quality changes with submission year and country.

Submission year: the GenBank "Direct Submission" date, or for BOLD-only
records the sequence upload date. For each period and for each country
the module reports the number of records and the share that
  - fail the sequence quality check,
  - carry no species name,
  - have coordinates, a voucher, and a linked publication,
  - disagree with the anchored species of their cluster (among records in
    anchored clusters),
together with the mean metadata score (0-7). A Spearman rank correlation
with a permutation p-value tests the year trend of the metadata score.
Because records deposited together are not independent, the trend is also
tested on submission batches (median year against mean score; see
study_keys.py), and a 95% interval for the record-level rho is taken from a
bootstrap that resamples batches. Countries are reported with the number
of batches behind their records.

Reads:  data/ncbi_*.gb, data/bold_*.tsv, results/genus/{merged_records,qc_records}.csv,
        results/anchors/record_anchor_check.csv
Writes: results/trends/by_year.csv, by_period.csv, by_country.csv, trends_summary.txt

Usage: python scripts/trends.py
"""
import argparse
import csv
import glob
import os
import re
from collections import defaultdict

import sys

import numpy as np
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from study_keys import genbank_keys, record_keys  # noqa: E402

PERIODS = [(0, 2005, "before 2006"), (2006, 2010, "2006-2010"), (2011, 2015, "2011-2015"),
           (2016, 2020, "2016-2020"), (2021, 2100, "2021 onwards")]
# BOLD and GenBank spell some countries differently; GenBank's form is kept.
COUNTRY_NAMES = {"United States": "USA", "United States of America": "USA", "Turkiye": "Turkey",
                 "Türkiye": "Turkey", "Russian Federation": "Russia", "Korea, Republic of": "South Korea",
                 "Republic of Korea": "South Korea", "Vietnam": "Viet Nam", "Iran, Islamic Republic of": "Iran",
                 "United Kingdom": "United Kingdom", "UK": "United Kingdom"}
SUBMITTED = re.compile(r"Submitted \(\d{2}-[A-Z]{3}-(\d{4})\)")


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def submission_years(data_dir):
    years = {}
    for path in glob.glob(os.path.join(data_dir, "ncbi_*.gb")):
        for rec in SeqIO.parse(path, "genbank"):
            for ref in rec.annotations.get("references", []):
                m = SUBMITTED.search(ref.journal or "")
                if m:
                    years[rec.id.split(".")[0]] = int(m.group(1))
                    break
    bold = {}
    for path in glob.glob(os.path.join(data_dir, "bold_*.tsv")):
        for r in read_tsv(path):
            y = (r.get("sequence_upload_date") or r.get("processid_minted_date") or "")[:4]
            if y.isdigit():
                bold[r["processid"]] = int(y)
    return years, bold


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def spearman(x, y, perms=2000, seed=1):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        return float("nan"), float("nan")
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    r = np.corrcoef(rx, ry)[0, 1]
    rng = np.random.default_rng(seed)
    null = [np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(perms)]
    p = (np.sum(np.abs(null) >= abs(r)) + 1) / (perms + 1)
    return float(r), float(p)


def summarise(group):
    n = len(group)

    def share(test):
        return f"{100 * sum(1 for g in group if test(g)) / n:.1f}" if n else ""

    anchored = [g for g in group if g["verdict"] in ("AGREES", "AGREES_SHARED_CLUSTER", "DISAGREES")]
    dis = [g for g in anchored if g["verdict"] == "DISAGREES"]
    return {"records": n, "batches": len({g["batch"] for g in group}),
            "qc_fail_pct": share(lambda g: g["qc"] == "FAIL"),
            "unnamed_pct": share(lambda g: not g["named"]),
            "coordinates_pct": share(lambda g: "coordinates" not in g["missing"]),
            "voucher_pct": share(lambda g: "voucher" not in g["missing"]),
            "publication_pct": share(lambda g: "publication" not in g["missing"]),
            "mean_metadata_score": f"{np.mean([g['score'] for g in group]):.2f}" if n else "",
            "in_anchored_clusters": len(anchored),
            "disagree_with_anchor_pct": f"{100 * sum(g['verdict'] == 'DISAGREES' for g in anchored) / len(anchored):.1f}"
            if anchored else "",
            "disagreeing_records": len(dis), "batches_with_disagreement": len({g["batch"] for g in dis})}


def main():
    p = argparse.ArgumentParser(description="Data quality by submission year and country")
    p.add_argument("--data", default="data")
    p.add_argument("--genus", default=os.path.join("results", "genus"))
    p.add_argument("--anchors-out", default=os.path.join("results", "anchors"))
    p.add_argument("--out", default=os.path.join("results", "trends"))
    p.add_argument("--min-country", type=int, default=20)
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    gb_years, bold_years = submission_years(args.data)
    merged = read(os.path.join(args.genus, "merged_records.csv"))
    qc = {r["record_id"]: r["qc_status"] for r in read(os.path.join(args.genus, "qc_records.csv"))}
    verdict_path = os.path.join(args.anchors_out, "record_anchor_check.csv")
    verdicts = {r["record_id"]: r["verdict"] for r in read(verdict_path)} if os.path.exists(verdict_path) else {}

    batch = record_keys(merged, genbank_keys(os.path.join(args.data, "ncbi_*.gb")))
    rows = []
    for r in merged:
        year = gb_years.get(r["genbank_accession"]) or bold_years.get(r["bold_processid"])
        words = r["species_label"].split()
        named = len(words) >= 2 and words[1].islower() and words[1] not in ("sp.", "sp", "cf.", "aff.")
        rows.append({"record_id": r["record_id"], "year": year,
                     "country": COUNTRY_NAMES.get(r["country"].split(":")[0].strip(),
                                                  r["country"].split(":")[0].strip()) or "(not given)",
                     "qc": qc.get(r["record_id"], ""), "named": named,
                     "missing": r["metadata_missing"].split(";"), "score": int(r["metadata_score"]),
                     "verdict": verdicts.get(r["record_id"], ""), "batch": batch[r["record_id"]]})

    with_year = [g for g in rows if g["year"]]
    by_year = defaultdict(list)
    for g in with_year:
        by_year[g["year"]].append(g)
    year_rows = [{"year": y, **summarise(gs)} for y, gs in sorted(by_year.items())]
    period_rows = []
    for lo, hi, label in PERIODS:
        gs = [g for g in with_year if lo <= g["year"] <= hi]
        if gs:
            period_rows.append({"period": label, **summarise(gs)})
    by_country = defaultdict(list)
    for g in rows:
        by_country[g["country"]].append(g)
    country_rows = sorted(({"country": c, **summarise(gs)} for c, gs in by_country.items()
                           if len(gs) >= args.min_country or c == "India"),
                          key=lambda r: -r["records"])

    for name, data in (("by_year.csv", year_rows), ("by_period.csv", period_rows), ("by_country.csv", country_rows)):
        with open(os.path.join(args.out, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0].keys()))
            w.writeheader()
            w.writerows(data)

    rho, pval = spearman([g["year"] for g in with_year], [g["score"] for g in with_year])
    groups = defaultdict(list)
    for g in with_year:
        groups[g["batch"]].append(g)
    b_year = [float(np.median([g["year"] for g in gs])) for gs in groups.values()]
    b_score = [float(np.mean([g["score"] for g in gs])) for gs in groups.values()]
    b_rho, b_p = spearman(b_year, b_score)
    rng = np.random.default_rng(20261030)
    keys = list(groups)
    boot = []
    for _ in range(2000):
        pick = [g for k in rng.choice(len(keys), size=len(keys), replace=True) for g in groups[keys[k]]]
        x = np.argsort(np.argsort([g["year"] for g in pick]))
        y = np.argsort(np.argsort([g["score"] for g in pick]))
        boot.append(np.corrcoef(x, y)[0, 1])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    lines = [f"Records: {len(rows)}; with a submission year: {len(with_year)} "
             f"({min(by_year)}-{max(by_year)})",
             f"Metadata score vs submission year: Spearman rho {rho:.3f}, permutation p {pval:.4f}; "
             f"95% interval resampling {len(keys)} submission batches {lo:.3f} to {hi:.3f}",
             f"Batch level (median year vs mean score, {len(keys)} batches): rho {b_rho:.3f}, permutation p {b_p:.4f}",
             "By period: records | QC fail % | unnamed % | coordinates % | voucher % | publication % | "
             "mean score | disagree with anchor %"]
    for r in period_rows:
        lines.append(f"  {r['period']}: {r['records']} | {r['qc_fail_pct']} | {r['unnamed_pct']} | "
                     f"{r['coordinates_pct']} | {r['voucher_pct']} | {r['publication_pct']} | "
                     f"{r['mean_metadata_score']} | {r['disagree_with_anchor_pct'] or '-'}"
                     f" (of {r['in_anchored_clusters']})")
    lines.append(f"By country (at least {args.min_country} records, plus India), same columns:")
    for r in country_rows:
        lines.append(f"  {r['country']}: {r['records']} | {r['qc_fail_pct']} | {r['unnamed_pct']} | "
                     f"{r['coordinates_pct']} | {r['voucher_pct']} | {r['publication_pct']} | "
                     f"{r['mean_metadata_score']} | {r['disagree_with_anchor_pct'] or '-'}"
                     f" (of {r['in_anchored_clusters']}) | batches {r['batches']}; disagreeing records "
                     f"{r['disagreeing_records']} from {r['batches_with_disagreement']} batches")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "trends_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
