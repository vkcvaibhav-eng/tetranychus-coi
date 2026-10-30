"""Module 7: how well does the public COI database identify a specimen?

Queries are the anchor records (morphologically identified specimens or
populations, evidence levels A and B, set per record when a study mixes
evidence). Each query is identified by its closest COI match within the
match distance, the way a user of GenBank or BOLD would. Four reference sets
are compared, all built from the records left after holding out the query's
study:

  all           every record that passed quality control
  named         records with a species name
  cleaned       named records minus those whose name differs from the
                anchored species (best evidence level present) of their cluster
  conservative  named records minus those whose name differs from every
                anchored species of their cluster, and only where one of those
                names is supported by anchors from at least two studies

Outcomes:
  correct          the closest match(es) carry the query's species name
  wrong            the closest match(es) carry one other species name
  ambiguous        the closest matches carry different species names
  nearest unnamed  the closest match(es) carry no species name ('all' only)
  no match         nothing within the match distance

For every query the coverage of its species in the named reference set is
recorded: no record of the species, records that do not overlap the query
by the minimum number of sites, or an overlapping record (with the distance
to the nearest one).

Leave-one-study-out with fold isolation: when the anchors of one study (or of
a group of related studies, data/anchors/study_groups.csv) are queried, every
record of that study is removed, and the clusters and the cleaning are rebuilt
from the remaining records, using the anchors of the other studies. A record
belongs to a study when it is one of the study's anchor records (any level),
cites the same publication title as one of them, or was deposited in the same
submission batch (study_keys.py). The query is judged by the name given in its
study (anchor_records.csv, column species), while reference records keep the
name in the database.

The earlier, less strict version (clusters from all records, membership by
publication title only) is also run to measure leakage.

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

OUTCOMES = ["correct", "wrong", "ambiguous", "nearest unnamed", "no match"]
SETS = ["all", "named", "cleaned", "conservative"]
COVERAGE = ["no record of the species", "no overlapping record", "overlapping record"]


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
    anchor_all, anchor_name, anchor_level = {}, {}, {}
    for a in read(os.path.join(anchors_dir, "anchor_records.csv")):
        rid = by_acc.get(a["accession"])
        if rid and rid in checks:
            anchor_all[rid] = a["study_id"]
            # evidence level of the record: its own when the study mixes evidence, else the study's
            anchor_level[rid] = (a.get("evidence_level") or "").strip() or studies[a["study_id"]]["evidence_level"]
            if (a.get("species") or "").strip():  # name given in the study, when the record's label differs
                anchor_name[rid] = a["species"].strip()

    gpath = os.path.join(anchors_dir, "study_groups.csv")  # optional: studies held out together
    groups = {r["study_id"]: r["group"] for r in read(gpath)} if os.path.exists(gpath) else {}

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
            "synonyms": synonyms, "anchor_all": anchor_all, "anchor_name": anchor_name,
            "anchor_level": anchor_level, "groups": groups, "hap_seqs": hap_seqs, "hap_of": hap_of,
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


def anchors_at(ctx, levels):
    """Anchor records whose own evidence level is in levels."""
    return {rid: sid for rid, sid in ctx["anchor_all"].items() if ctx["anchor_level"][rid] in levels}


def evaluate(ctx, t_cluster=0.03, t_match=0.03, levels="AB", synonyms=None, isolate=True, rule="title+batch",
             groups=None):
    """Leave-one-study-out identification of anchor records; one row per query."""
    syn = ctx["synonyms"] if synonyms is None else synonyms
    records, rec_hap, D = ctx["records"], ctx["rec_hap"], ctx["D"]
    lvl = ctx["anchor_level"]
    db_name = {rid: syn.get(n, n) for rid, n in ctx["raw_name"].items()}  # names as the database holds them
    names_arr = np.array([db_name[r] for r in records], dtype=object)
    named = names_arr != ""
    truth = dict(db_name)  # anchors are judged by the name given in their study
    for rid, n in ctx.get("anchor_name", {}).items():
        truth[rid] = syn.get(n, n)
    anchors = anchors_at(ctx, levels)
    used = set(anchors.values())
    members = study_members(ctx, {r: s for r, s in ctx["anchor_all"].items() if s in used}, rule)
    groups = {} if groups is None else groups  # studies held out together (sensitivity run)
    folds = defaultdict(set)
    for sid in used:
        folds[groups.get(sid, sid)].add(sid)
    shared_clusters = None if isolate else cluster_ids(ctx, t_cluster)

    def identify(q, mask, allow_unnamed):
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
            return ("nearest unnamed" if allow_unnamed else "no match"), "", float(dmin)
        if len(c) > 1:
            return "ambiguous", "; ".join(n for n, _ in c.most_common()), float(dmin)
        top = next(iter(c))
        return ("correct" if top == truth[q] else "wrong"), top, float(dmin)

    def fmt(x):
        return "" if np.isinf(x) else f"{abs(max(x, 0.0)):.4f}"

    rows = []
    for fold, sids in sorted(folds.items()):
        held = set().union(*(members[s] for s in sids))
        train = np.array([r not in held for r in records])
        cl = cluster_ids(ctx, t_cluster, train) if isolate else shared_clusters
        best = defaultdict(lambda: defaultdict(set))   # cluster -> level -> names
        support = defaultdict(lambda: defaultdict(set))  # cluster -> name -> studies
        for rid, s2 in anchors.items():
            if s2 not in sids and rid not in held and truth[rid]:
                best[cl[rid]][lvl[rid]].add(truth[rid])
                support[cl[rid]][truth[rid]].add(s2)
        keep_c = train & named
        keep_k = train & named
        for k, rid in enumerate(records):
            if not keep_c[k]:
                continue
            v = best.get(cl[rid])
            if v and db_name[rid] not in v[min(v)]:
                keep_c[k] = False
            sup = support.get(cl[rid])
            if sup and db_name[rid] not in sup and any(len(st) >= 2 for st in sup.values()):
                keep_k[k] = False
        masks = {"all": train, "named": train & named, "cleaned": keep_c, "conservative": keep_k}
        for q in sorted(r for r, s2 in anchors.items() if s2 in sids):
            row = {"query": q, "study_id": anchors[q], "fold": fold, "evidence_level": lvl[q], "species": truth[q]}
            for name_, m in masks.items():
                res, match, dist = identify(q, m, name_ == "all")
                row.update({f"{name_}_result": res, f"{name_}_match": match, f"{name_}_distance": fmt(dist)})
            for name_ in ("named", "cleaned"):
                same = masks[name_] & (names_arr == truth[q])
                d = D[ctx["hap_of"][q], rec_hap][same]
                if not same.any():
                    cov, near = COVERAGE[0], np.inf
                elif np.all(np.isnan(d)):
                    cov, near = COVERAGE[1], np.inf
                else:
                    cov, near = COVERAGE[2], float(np.nanmin(d))
                row.update({f"coverage_{name_}": cov, f"nearest_same_species_{name_}": fmt(near)})
            row.update({"held_out_records": len(held), **{f"{n}_db_size": int(m.sum()) for n, m in masks.items()},
                        **{f"{n}_db_species": len(set(names_arr[m]) - {""}) for n, m in masks.items()}})
            rows.append(row)
    return rows


def tally(rows, key):
    c = Counter(r[key] for r in rows)
    return {o: c[o] for o in OUTCOMES}


def describe(rows, key):
    n = len(rows)
    c = tally(rows, key)
    return ", ".join(f"{o} {c[o]} ({100 * c[o] / n:.1f}%)" for o in OUTCOMES if c[o] or o != "nearest unnamed") \
        if n else "none"


def main():
    p = argparse.ArgumentParser(description="Identification accuracy against four reference sets")
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

    anchors = anchors_at(ctx, args.levels)
    used = set(anchors.values())
    allstudy = {r: s for r, s in ctx["anchor_all"].items() if s in used}
    m_title = study_members(ctx, allstudy, "title")
    m_batch = study_members(ctx, allstudy, "title+batch")
    mem_rows = [{"study_id": sid, "anchor_queries": sum(1 for s in anchors.values() if s == sid),
                 "held_out_by_title": len(m_title[sid]), "held_out_by_title_and_batch": len(m_batch[sid]),
                 "added_by_batch": " ".join(sorted(m_batch[sid] - m_title[sid]))} for sid in sorted(used)]
    with open(os.path.join(args.out, "membership.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mem_rows[0].keys()))
        w.writeheader()
        w.writerows(mem_rows)

    n = len(rows)
    lines = [f"Queries: {n} anchor records from {len(used)} studies (leave-one-study-out with clusters and "
             f"cleaning rebuilt in each fold; clusters at K2P {args.cluster}, matches within {args.match})"]
    for s_ in SETS:
        lines.append(f"  {s_:>12}: {describe(rows, s_ + '_result')}")
    cov = Counter(r["coverage_named"] for r in rows)
    lines.append("Coverage of the query's species in the named reference set: "
                 + "; ".join(f"{c} {cov[c]}" for c in COVERAGE))
    for c in COVERAGE:
        rs = [r for r in rows if r["coverage_named"] == c]
        if rs:
            lines.append(f"  {c} ({len(rs)}): named -> {describe(rs, 'named_result')}")
    lines.append("Earlier version (clusters from all records, study membership by publication title only):")
    for s_ in SETS:
        lines.append(f"  {s_:>12}: {describe(legacy, s_ + '_result')}")
    changed = sum(1 for a, b in zip(rows, legacy) if any(a[f"{s_}_result"] != b[f"{s_}_result"] for s_ in SETS))
    lines.append(f"Queries whose outcome differs from the earlier version: {changed}")
    lines.append("Records held out per study (publication title -> title or submission batch):")
    lines += [f"  {r['study_id']}: {r['held_out_by_title']} -> {r['held_out_by_title_and_batch']}" for r in mem_rows]
    lines.append("By study (correct: all / named / cleaned / conservative; wrong likewise; of queries):")
    per_study = defaultdict(list)
    for r in rows:
        per_study[r["study_id"]].append(r)
    for sid, rs in sorted(per_study.items()):
        c = [tally(rs, s_ + "_result") for s_ in SETS]
        lines.append(f"  {sid}: {' / '.join(str(x['correct']) for x in c)}; {' / '.join(str(x['wrong']) for x in c)};"
                     f" of {len(rs)}")
    lines.append("By species (correct: all / named / cleaned / conservative, of queries):")
    per = defaultdict(list)
    for r in rows:
        per[r["species"]].append(r)
    for sp, rs in sorted(per.items(), key=lambda kv: -len(kv[1])):
        lines.append(f"  {sp}: " + " / ".join(str(sum(r[s_ + '_result'] == 'correct' for r in rs)) for s_ in SETS)
                     + f" of {len(rs)}")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "identification_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
