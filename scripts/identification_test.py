"""Module 7: how well does the public COI database identify a specimen?

Queries are the anchor records (morphologically identified specimens,
evidence levels A and B). Each query is identified by its closest COI match
in a reference database, the way a user of GenBank or BOLD would:

  correct      the closest match(es) carry the query's species name
  wrong        the closest match(es) carry one other species name
  ambiguous    the closest matches carry different species names
  no match     nothing within --max-distance (K2P)

Two databases are compared:
  raw      every record that passed the quality check, with its public name
  cleaned  raw minus unnamed records and minus records whose name
           disagrees with the anchored species of their distance cluster

Leave-one-study-out: when the anchors of one study are queried, every record
from that study is removed from both databases, and that study's anchors are
not used to clean the database. Names are compared after the published
synonymies in data/anchors/synonyms.csv.

Reads:  results/genus/{merged_records.csv,qc_passed.fasta},
        results/species/record_checks.csv, data/anchors/*, reference COX1
Writes: results/identification/queries.csv, identification_summary.txt

Usage: python scripts/identification_test.py [--max-distance 0.03]
"""
import argparse
import csv
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_qc import make_aligner  # noqa: E402
from species_analysis import encode, k2p_row, project, species_name  # noqa: E402


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser(description="Identification accuracy of raw vs cleaned COI databases")
    p.add_argument("--genus", default=os.path.join("results", "genus"))
    p.add_argument("--species", default=os.path.join("results", "species"))
    p.add_argument("--anchors", default=os.path.join("data", "anchors"))
    p.add_argument("--reference", default=os.path.join("data", "reference", "coi_reference.fasta"))
    p.add_argument("--out", default=os.path.join("results", "identification"))
    p.add_argument("--levels", default="AB")
    p.add_argument("--min-overlap", type=int, default=300)
    p.add_argument("--max-distance", type=float, default=0.03)
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    syn_path = os.path.join(args.anchors, "synonyms.csv")
    synonyms = {r["name"]: r["accepted_name"] for r in read(syn_path)} if os.path.exists(syn_path) else {}

    def canon(label):
        n = species_name(label)
        return synonyms.get(n, n)

    merged = {r["record_id"]: r for r in read(os.path.join(args.genus, "merged_records.csv"))}
    checks = {r["record_id"]: r for r in read(os.path.join(args.species, "record_checks.csv"))}
    studies = {s["study_id"]: s for s in read(os.path.join(args.anchors, "anchor_studies.csv"))}
    by_acc = {r["genbank_accession"]: rid for rid, r in merged.items() if r["genbank_accession"]}
    anchor_study = {}
    for a in read(os.path.join(args.anchors, "anchor_records.csv")):
        rid = by_acc.get(a["accession"])
        if rid and rid in checks and studies[a["study_id"]]["evidence_level"] in args.levels:
            anchor_study[rid] = a["study_id"]

    # Records from the same study as an anchor study: matched by publication title
    # of that study's anchors, so unlisted records of the study are held out too.
    study_titles = defaultdict(set)
    for rid, sid in anchor_study.items():
        t = merged[rid]["publication"].strip().lower()
        if t:
            study_titles[sid].add(t)

    def study_of(rid):
        if rid in anchor_study:
            return anchor_study[rid]
        t = merged[rid]["publication"].strip().lower()
        for sid, titles in study_titles.items():
            if t and t in titles:
                return sid
        return ""

    # Haplotype distances (same placement and K2P as species_analysis.py).
    ref = str(SeqIO.read(args.reference, "fasta").seq).upper()
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(args.genus, "qc_passed.fasta"), "fasta")}
    aligner = make_aligner()
    hap_index, hap_seqs, hap_of = {}, [], {}
    for rid, s in seqs.items():
        proj = project(aligner, ref, s)
        if proj not in hap_index:
            hap_index[proj] = len(hap_seqs)
            hap_seqs.append(proj)
        hap_of[rid] = hap_index[proj]
    M = np.vstack([encode(s) for s in hap_seqs])
    D = np.vstack([k2p_row(M[i], M, args.min_overlap) for i in range(len(hap_seqs))])
    np.fill_diagonal(D, 0.0)
    print(f"{len(seqs)} sequences, {len(hap_seqs)} haplotypes, {len(anchor_study)} anchor queries")

    records = [rid for rid in seqs if rid in checks]
    name = {rid: canon(merged[rid]["species_label"]) for rid in records}
    cluster = {rid: checks[rid]["cluster"] for rid in records}
    study = {rid: study_of(rid) for rid in records}

    def cleaned_db(held_out):
        """Records kept after cleaning with anchors from all studies except held_out."""
        anchored = defaultdict(lambda: defaultdict(set))
        for rid, sid in anchor_study.items():
            if sid != held_out and name[rid]:
                anchored[cluster[rid]][studies[sid]["evidence_level"]].add(name[rid])
        keep = []
        for rid in records:
            if not name[rid]:
                continue
            levels = anchored.get(cluster[rid])
            if levels:
                best = levels[min(levels)]
                if name[rid] not in best:
                    continue
            keep.append(rid)
        return keep

    def identify(q, db):
        hq = hap_of[q]
        best_d, hits = math.inf, []
        for r in db:
            d = D[hq, hap_of[r]]
            if np.isnan(d):
                continue
            if d < best_d - 1e-12:
                best_d, hits = float(d), [r]
            elif abs(d - best_d) <= 1e-12:
                hits.append(r)
        if not hits or best_d > args.max_distance:
            return "no match", "", best_d
        names = Counter(name[r] for r in hits if name[r])
        if not names:
            return "no match", "", best_d
        if len(names) > 1:
            return "ambiguous", "; ".join(n for n, _ in names.most_common()), best_d
        top = next(iter(names))
        return ("correct" if top == name[q] else "wrong"), top, best_d

    rows = []
    cleaned_cache = {}
    for q, sid in sorted(anchor_study.items(), key=lambda kv: kv[1]):
        raw = [r for r in records if r != q and study[r] != sid]
        if sid not in cleaned_cache:
            cleaned_cache[sid] = cleaned_db(sid)
        clean = [r for r in cleaned_cache[sid] if r != q and study[r] != sid]
        res_raw, hit_raw, d_raw = identify(q, raw)
        res_cln, hit_cln, d_cln = identify(q, clean)
        rows.append({"query": q, "study_id": sid, "evidence_level": studies[sid]["evidence_level"],
                     "species": name[q], "raw_result": res_raw, "raw_match": hit_raw,
                     "raw_distance": "" if math.isinf(d_raw) else f"{max(d_raw, 0.0):.4f}",
                     "cleaned_result": res_cln, "cleaned_match": hit_cln,
                     "cleaned_distance": "" if math.isinf(d_cln) else f"{max(d_cln, 0.0):.4f}",
                     "raw_db_size": len(raw), "cleaned_db_size": len(clean)})

    with open(os.path.join(args.out, "queries.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["query"])
        w.writeheader()
        w.writerows(rows)

    order = ["correct", "wrong", "ambiguous", "no match"]
    n = len(rows)

    def line(key, subset):
        c = Counter(r[key] for r in subset)
        k = len(subset)
        return ", ".join(f"{o} {c[o]} ({100 * c[o] / k:.1f}%)" for o in order) if k else "none"

    lines = [f"Queries: {n} anchor records from {len({r['study_id'] for r in rows})} studies "
             f"(leave-one-study-out; closest match within K2P {args.max_distance})",
             f"Database size: raw about {rows[0]['raw_db_size'] if rows else 0}, "
             f"cleaned about {rows[0]['cleaned_db_size'] if rows else 0} records",
             f"Raw database:     {line('raw_result', rows)}",
             f"Cleaned database: {line('cleaned_result', rows)}"]
    excl = [r for r in rows if r["species"] not in ("Tetranychus turkestani",)]
    lines.append("Excluding T. turkestani (not separable from T. urticae by COI):")
    lines.append(f"  Raw:     {line('raw_result', excl)}")
    lines.append(f"  Cleaned: {line('cleaned_result', excl)}")
    lines.append("By study (raw correct -> cleaned correct, of queries):")
    per_study = defaultdict(list)
    for r in rows:
        per_study[(r["study_id"], r["evidence_level"])].append(r)
    for (sid, lvl), rs in sorted(per_study.items()):
        lines.append(f"  {sid} [{lvl}]: {sum(r['raw_result'] == 'correct' for r in rs)} -> "
                     f"{sum(r['cleaned_result'] == 'correct' for r in rs)} of {len(rs)}")
    lines.append("By species (raw correct -> cleaned correct, of queries):")
    per = defaultdict(list)
    for r in rows:
        per[r["species"]].append(r)
    for sp, rs in sorted(per.items(), key=lambda kv: -len(kv[1])):
        lines.append(f"  {sp}: {sum(r['raw_result'] == 'correct' for r in rs)} -> "
                     f"{sum(r['cleaned_result'] == 'correct' for r in rs)} of {len(rs)}")
    changed = Counter((r["raw_result"], r["cleaned_result"]) for r in rows if r["raw_result"] != r["cleaned_result"])
    lines.append("Changes from raw to cleaned:")
    lines += [f"  {a} -> {b}: {k}" for (a, b), k in changed.most_common()] or ["  none"]
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "identification_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
