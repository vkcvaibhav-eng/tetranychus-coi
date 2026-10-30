"""Module 7: how well does the public COI database identify a specimen?

Queries are the anchor records (morphologically identified specimens,
evidence levels A and B). Each query is identified by its closest COI match
in a reference database, the way a user of GenBank or BOLD would:

  correct      the closest match(es) carry the query's species name
  wrong        the closest match(es) carry one other species name
  ambiguous    the closest matches carry different species names
  no match     nothing within the match distance, or only unnamed records

Two databases are compared:
  raw      every record that passed the quality check, with its public name
  cleaned  raw minus unnamed records and minus records whose name
           disagrees with the anchored species of their distance cluster

Leave-one-study-out with fold isolation: when the anchors of one study are
queried, every record of that study is removed, and the clusters and the
cleaning are rebuilt from the remaining (training) records only, using the
anchors of the other studies. A record belongs to a study when it is one of
the study's anchors, cites the same publication title as one of them, or was
deposited in the same submission batch (study_keys.py) as one of them.

For comparison the earlier, less strict version is also run: clusters from
all records, and study membership by publication title only.

Names are compared after the published synonymies in data/anchors/synonyms.csv.

Reads:  results/genus/{merged_records.csv,qc_passed.fasta},
        results/species/record_checks.csv, data/anchors/*, data/ncbi_*.gb,
        reference COX1
Writes: results/identification/queries.csv, membership.csv,
        identification_summary.txt

Usage: python scripts/identification_test.py [--cluster 0.03] [--match 0.03]
"""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_qc import make_aligner  # noqa: E402
from species_analysis import UnionFind, encode, k2p_row, project, species_name  # noqa: E402
from study_keys import genbank_keys, record_keys  # noqa: E402

OUTCOMES = ["correct", "wrong", "ambiguous", "no match"]


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def load_context(genus, species_dir, anchors_dir, reference, data_dir, min_overlap=300):
    """Everything the identification test needs, computed once."""
    merged_rows = read(os.path.join(genus, "merged_records.csv"))
    merged = {r["record_id"]: r for r in merged_rows}
    checks = {r["record_id"]: r for r in read(os.path.join(species_dir, "record_checks.csv"))}
    studies = {s["study_id"]: s for s in read(os.path.join(anchors_dir, "anchor_studies.csv"))}
    synonyms = {r["name"]: r["accepted_name"] for r in read(os.path.join(anchors_dir, "synonyms.csv"))}
    by_acc = {r["genbank_accession"]: rid for rid, r in merged.items() if r["genbank_accession"]}
    anchor_all = {}
    for a in read(os.path.join(anchors_dir, "anchor_records.csv")):
        rid = by_acc.get(a["accession"])
        if rid and rid in checks:
            anchor_all[rid] = a["study_id"]

    ref = str(SeqIO.read(reference, "fasta").seq).upper()
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(genus, "qc_passed.fasta"), "fasta")}
    aligner = make_aligner()
    hap_index, hap_seqs, hap_of = {}, [], {}
    for rid, s in seqs.items():
        proj = project(aligner, ref, s)
        if proj not in hap_index:
            hap_index[proj] = len(hap_seqs)
            hap_seqs.append(proj)
        hap_of[rid] = hap_index[proj]
    M = np.vstack([encode(s) for s in hap_seqs])
    D = np.vstack([k2p_row(M[i], M, min_overlap) for i in range(len(hap_seqs))])
    np.fill_diagonal(D, 0.0)
    records = [rid for rid in seqs if rid in checks]
    return {"merged": merged, "merged_rows": merged_rows, "checks": checks, "studies": studies,
            "synonyms": synonyms, "anchor_all": anchor_all, "hap_seqs": hap_seqs, "hap_of": hap_of,
            "D": D, "Dfin": np.nan_to_num(D, nan=np.inf), "records": records,
            "index": {rid: i for i, rid in enumerate(records)},
            "rec_hap": np.array([hap_of[r] for r in records]),
            "raw_name": {rid: species_name(merged[rid]["species_label"]) for rid in records},
            "batch": record_keys(merged_rows, genbank_keys(os.path.join(data_dir, "ncbi_*.gb")))}


def cluster_ids(ctx, t, mask=None):
    """Single-linkage clusters at K2P t over the records selected by mask (all if None)."""
    rec_hap = ctx["rec_hap"]
    sel = np.ones(len(rec_hap), bool) if mask is None else mask
    haps = np.unique(rec_hap[sel])
    pos = {h: k for k, h in enumerate(haps)}
    uf = UnionFind(len(haps))
    sub = ctx["Dfin"][np.ix_(haps, haps)] <= t
    for i, j in np.argwhere(np.triu(sub, k=1)):
        uf.union(int(i), int(j))
    return {rid: int(haps[uf.find(pos[rec_hap[k]])])
            for k, rid in enumerate(ctx["records"]) if sel[k]}


def study_members(ctx, anchors, rule):
    """Records belonging to each anchor study: 'title' or 'title+batch'."""
    merged, batch = ctx["merged"], ctx["batch"]
    titles, batches = defaultdict(set), defaultdict(set)
    for rid, sid in anchors.items():
        t = merged[rid]["publication"].strip().lower()
        if t:
            titles[sid].add(t)
        batches[sid].add(batch.get(rid, ""))
    members = defaultdict(set)
    for rid in ctx["records"]:
        t = merged[rid]["publication"].strip().lower()
        for sid in set(anchors.values()):
            if anchors.get(rid) == sid or (t and t in titles[sid]) or \
                    (rule == "title+batch" and batch.get(rid, "") in batches[sid]):
                members[sid].add(rid)
    return members


def evaluate(ctx, t_cluster=0.03, t_match=0.03, levels="AB", synonyms=None, isolate=True, rule="title+batch"):
    """Leave-one-study-out identification of anchor specimens; one row per query."""
    syn = ctx["synonyms"] if synonyms is None else synonyms
    records, idx, rec_hap, D = ctx["records"], ctx["index"], ctx["rec_hap"], ctx["D"]
    studies = ctx["studies"]
    name = {rid: syn.get(n, n) for rid, n in ctx["raw_name"].items()}
    names_arr = np.array([name[r] for r in records], dtype=object)
    named = names_arr != ""
    anchors = {rid: sid for rid, sid in ctx["anchor_all"].items() if studies[sid]["evidence_level"] in levels}
    members = study_members(ctx, anchors, rule)
    shared_clusters = None if isolate else cluster_ids(ctx, t_cluster)

    def identify(q, mask):
        d = D[ctx["hap_of"][q], rec_hap]
        ok = mask & ~np.isnan(d)
        if not ok.any():
            return "no match", "", np.inf
        dmin = d[ok].min()
        if dmin > t_match:
            return "no match", "", float(dmin)
        hits = np.where(ok & (np.abs(d - dmin) <= 1e-12))[0]
        c = Counter(names_arr[i] for i in hits if names_arr[i])
        if not c:
            return "no match", "", float(dmin)
        if len(c) > 1:
            return "ambiguous", "; ".join(n for n, _ in c.most_common()), float(dmin)
        top = next(iter(c))
        return ("correct" if top == name[q] else "wrong"), top, float(dmin)

    rows = []
    for sid in sorted(set(anchors.values())):
        held = members[sid]
        train = np.array([r not in held for r in records])
        cl = cluster_ids(ctx, t_cluster, train) if isolate else shared_clusters
        best = defaultdict(lambda: defaultdict(set))
        for rid, s2 in anchors.items():
            if s2 != sid and rid not in held and name[rid]:
                best[cl[rid]][studies[s2]["evidence_level"]].add(name[rid])
        keep = train & named
        for k, rid in enumerate(records):
            if keep[k]:
                v = best.get(cl[rid])
                if v and name[rid] not in v[min(v)]:
                    keep[k] = False
        for q in sorted(r for r, s2 in anchors.items() if s2 == sid):
            rr, rm, rd = identify(q, train)
            cr, cm, cd = identify(q, keep)
            rows.append({"query": q, "study_id": sid, "evidence_level": studies[sid]["evidence_level"],
                         "species": name[q], "raw_result": rr, "raw_match": rm,
                         "raw_distance": "" if np.isinf(rd) else f"{abs(max(rd, 0.0)):.4f}",
                         "cleaned_result": cr, "cleaned_match": cm,
                         "cleaned_distance": "" if np.isinf(cd) else f"{abs(max(cd, 0.0)):.4f}",
                         "held_out_records": len(held), "raw_db_size": int(train.sum()),
                         "cleaned_db_size": int(keep.sum())})
    return rows


def tally(rows, key):
    c = Counter(r[key] for r in rows)
    return {o: c[o] for o in OUTCOMES}


def describe(rows, key):
    n = len(rows)
    c = tally(rows, key)
    return ", ".join(f"{o} {c[o]} ({100 * c[o] / n:.1f}%)" for o in OUTCOMES) if n else "none"


def main():
    p = argparse.ArgumentParser(description="Identification accuracy of raw vs cleaned COI databases")
    p.add_argument("--genus", default=os.path.join("results", "genus"))
    p.add_argument("--species", default=os.path.join("results", "species"))
    p.add_argument("--anchors", default=os.path.join("data", "anchors"))
    p.add_argument("--reference", default=os.path.join("data", "reference", "coi_reference.fasta"))
    p.add_argument("--data", default="data")
    p.add_argument("--out", default=os.path.join("results", "identification"))
    p.add_argument("--levels", default="AB")
    p.add_argument("--min-overlap", type=int, default=300)
    p.add_argument("--cluster", type=float, default=0.03, help="K2P threshold for the clusters used in cleaning")
    p.add_argument("--match", type=float, default=0.03, help="Largest K2P distance accepted as a match")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ctx = load_context(args.genus, args.species, args.anchors, args.reference, args.data, args.min_overlap)
    rows = evaluate(ctx, args.cluster, args.match, args.levels)
    legacy = evaluate(ctx, args.cluster, args.match, args.levels, isolate=False, rule="title")
    with open(os.path.join(args.out, "queries.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    anchors = {rid: sid for rid, sid in ctx["anchor_all"].items()
               if ctx["studies"][sid]["evidence_level"] in args.levels}
    m_title = study_members(ctx, anchors, "title")
    m_batch = study_members(ctx, anchors, "title+batch")
    mem_rows = [{"study_id": sid, "anchors": sum(1 for s in anchors.values() if s == sid),
                 "held_out_by_title": len(m_title[sid]), "held_out_by_title_and_batch": len(m_batch[sid]),
                 "added_by_batch": " ".join(sorted(m_batch[sid] - m_title[sid]))}
                for sid in sorted(set(anchors.values()))]
    with open(os.path.join(args.out, "membership.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mem_rows[0].keys()))
        w.writeheader()
        w.writerows(mem_rows)

    n = len(rows)
    lines = [f"Queries: {n} anchor records from {len({r['study_id'] for r in rows})} studies "
             f"(leave-one-study-out with clusters and cleaning rebuilt in each fold; clusters at K2P "
             f"{args.cluster}, matches within {args.match})",
             f"Raw database:     {describe(rows, 'raw_result')}",
             f"Cleaned database: {describe(rows, 'cleaned_result')}",
             "Earlier version (clusters from all records, study membership by publication title only):",
             f"  Raw:     {describe(legacy, 'raw_result')}",
             f"  Cleaned: {describe(legacy, 'cleaned_result')}",
             "Records held out per study (publication title -> title or submission batch):"]
    lines += [f"  {r['study_id']}: {r['held_out_by_title']} -> {r['held_out_by_title_and_batch']}" for r in mem_rows]
    changed = [(a["query"], a["study_id"], f"{b['raw_result']}/{b['cleaned_result']}",
                f"{a['raw_result']}/{a['cleaned_result']}")
               for a, b in zip(rows, legacy) if (a["raw_result"], a["cleaned_result"]) != (b["raw_result"], b["cleaned_result"])]
    lines.append(f"Queries whose outcome differs from the earlier version: {len(changed)}")
    lines += [f"  {q} ({s}): raw/cleaned {old} -> {new}" for q, s, old, new in changed]
    excl = [r for r in rows if r["species"] != "Tetranychus turkestani"]
    lines.append("Excluding T. turkestani (not separable from T. urticae by COI):")
    lines.append(f"  Raw:     {describe(excl, 'raw_result')}")
    lines.append(f"  Cleaned: {describe(excl, 'cleaned_result')}")
    lines.append("By study (raw correct -> cleaned correct; raw wrong -> cleaned wrong; of queries):")
    per_study = defaultdict(list)
    for r in rows:
        per_study[(r["study_id"], r["evidence_level"])].append(r)
    for (sid, lvl), rs in sorted(per_study.items()):
        a, b = tally(rs, "raw_result"), tally(rs, "cleaned_result")
        lines.append(f"  {sid} [{lvl}]: {a['correct']} -> {b['correct']}; {a['wrong']} -> {b['wrong']}; of {len(rs)}")
    lines.append("By species (raw correct -> cleaned correct, of queries):")
    per = defaultdict(list)
    for r in rows:
        per[r["species"]].append(r)
    for sp, rs in sorted(per.items(), key=lambda kv: -len(kv[1])):
        lines.append(f"  {sp}: {sum(r['raw_result'] == 'correct' for r in rs)} -> "
                     f"{sum(r['cleaned_result'] == 'correct' for r in rs)} of {len(rs)}")
    moves = Counter((r["raw_result"], r["cleaned_result"]) for r in rows if r["raw_result"] != r["cleaned_result"])
    lines.append("Changes from raw to cleaned:")
    lines += [f"  {a} -> {b}: {k}" for (a, b), k in moves.most_common()] or ["  none"]
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "identification_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
