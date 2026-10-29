"""Module 10: sensitivity and uncertainty checks of the COI results.

1. Thresholds. Clusters, anchor verdicts and the identification test are
   repeated at several K2P thresholds (the same value is used for clusters
   and for the largest distance accepted as a match).
2. Anchors and taxonomy. At K2P 0.03 the identification test is repeated
   with level A anchors only, and with T. malaysiensis treated as a junior
   synonym of T. macfarlanei (Zeity et al. 2017) in addition to the two
   synonymies used throughout.
3. Uncertainty. For the main identification test: the full table of
   outcomes raw -> cleaned, an exact McNemar test on the records whose
   correctness changed, and 95% intervals from a bootstrap that resamples
   anchor studies (not records), since anchors of one study are not
   independent.
4. Cluster structure. Largest K2P distance within each cluster (single
   linkage can chain records that are far apart), and the clusters whose 5'
   and 3' fragments are joined only through records long enough to overlap
   both.
5. Label agreement comparable with ITS2. For every named COI record, the
   name of its nearest other record, (a) over all records and (b) over
   records from other submission batches only. ITS2 is scored the same way
   in its2_analysis.py; for specimens sequenced for both genes the two
   markers are compared on the same specimens.
6. Records that disagree with the anchored species of their cluster, by
   pair of names and by submission batch.

Reads:  results/genus/*, results/species/record_checks.csv, data/anchors/*,
        data/ncbi_*.gb, reference COX1, results/its2/* (optional)
Writes: results/robustness/*.csv, results/robustness/robustness_summary.txt

Usage: python scripts/robustness.py [--boot 10000]
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
from species_analysis import UnionFind, encode, k2p_row, project, species_name  # noqa: E402
from study_keys import genbank_keys, record_keys  # noqa: E402

OUTCOMES = ["correct", "wrong", "ambiguous", "no match"]
THRESHOLDS = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05]
EXTRA_SYNONYMS = {"Tetranychus malaysiensis": "Tetranychus macfarlanei"}


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["none"])
        w.writeheader()
        w.writerows(rows)


def binom_two_sided(k, n):
    """Exact two-sided binomial p-value for k of n at p = 0.5."""
    if n == 0:
        return 1.0
    probs = [math.comb(n, i) / 2 ** n for i in range(n + 1)]
    return min(1.0, sum(p for p in probs if p <= probs[k] + 1e-12))


def main():
    p = argparse.ArgumentParser(description="Sensitivity and uncertainty checks")
    p.add_argument("--genus", default=os.path.join("results", "genus"))
    p.add_argument("--species", default=os.path.join("results", "species"))
    p.add_argument("--anchors", default=os.path.join("data", "anchors"))
    p.add_argument("--its2", default=os.path.join("results", "its2"))
    p.add_argument("--data", default="data")
    p.add_argument("--reference", default=os.path.join("data", "reference", "coi_reference.fasta"))
    p.add_argument("--out", default=os.path.join("results", "robustness"))
    p.add_argument("--min-overlap", type=int, default=300)
    p.add_argument("--boot", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20261030)
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    base_syn = {r["name"]: r["accepted_name"] for r in read(os.path.join(args.anchors, "synonyms.csv"))}
    merged_rows = read(os.path.join(args.genus, "merged_records.csv"))
    merged = {r["record_id"]: r for r in merged_rows}
    qc = {r["record_id"]: r for r in read(os.path.join(args.genus, "qc_records.csv"))}
    checks = {r["record_id"]: r for r in read(os.path.join(args.species, "record_checks.csv"))}
    studies = {s["study_id"]: s for s in read(os.path.join(args.anchors, "anchor_studies.csv"))}
    by_acc = {r["genbank_accession"]: rid for rid, r in merged.items() if r["genbank_accession"]}
    anchor_all = {}
    for a in read(os.path.join(args.anchors, "anchor_records.csv")):
        rid = by_acc.get(a["accession"])
        if rid and rid in checks:
            anchor_all[rid] = a["study_id"]
    batch = record_keys(merged_rows, genbank_keys(os.path.join(args.data, "ncbi_*.gb")))

    # Haplotypes on reference coordinates and K2P distances (as species_analysis.py).
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
    H = len(hap_seqs)
    M = np.vstack([encode(s) for s in hap_seqs])
    D = np.vstack([k2p_row(M[i], M, args.min_overlap) for i in range(H)])
    np.fill_diagonal(D, 0.0)
    Dfin = np.nan_to_num(D, nan=np.inf)
    records = [rid for rid in seqs if rid in checks]
    print(f"{len(records)} records, {H} haplotypes")

    raw_name = {rid: species_name(merged[rid]["species_label"]) for rid in records}

    def clusters_at(t):
        uf = UnionFind(H)
        for i, j in np.argwhere(np.triu(Dfin <= t, k=1)):
            uf.union(int(i), int(j))
        return {rid: uf.find(hap_of[rid]) for rid in records}

    def run(t, levels, syn):
        """Clusters, anchor verdicts and identification test at threshold t."""
        name = {rid: syn.get(raw_name[rid], raw_name[rid]) for rid in records}
        cl = clusters_at(t)
        anchors = {rid: sid for rid, sid in anchor_all.items() if studies[sid]["evidence_level"] in levels}
        # cluster composition (submitted names, as in species_analysis.py)
        comp = defaultdict(Counter)
        for rid in records:
            comp[cl[rid]][raw_name[rid] or "(unnamed)"] += 1
        n_cl = len(comp)
        mixed = sum(1 for c in comp.values() if len([n for n in c if n != "(unnamed)"]) > 1)
        sp_clusters = defaultdict(set)
        for rid in records:
            if raw_name[rid]:
                sp_clusters[raw_name[rid]].add(cl[rid])
        split = sum(1 for v in sp_clusters.values() if len(v) > 1)
        cinna = len(sp_clusters.get("Tetranychus cinnabarinus", ()))
        urt = Counter(cl[r] for r in records if raw_name[r] == "Tetranychus urticae")
        urt_share = 100 * urt.most_common(1)[0][1] / sum(urt.values()) if urt else float("nan")

        # anchor verdicts for non-anchor records
        best = defaultdict(lambda: defaultdict(set))
        for rid, sid in anchors.items():
            if name[rid]:
                best[cl[rid]][studies[sid]["evidence_level"]].add(name[rid])
        anchored = {c: sorted(v[min(v)]) for c, v in best.items()}
        verdict = {}
        for rid in records:
            if rid in anchor_all:  # anchors of every level are left out, as in anchor_check.py
                continue
            names = anchored.get(cl[rid], [])
            own = name[rid]
            if not names:
                verdict[rid] = "NO_ANCHOR"
            elif not own:
                verdict[rid] = "UNNAMED_SUGGESTED" if len(names) == 1 else "UNNAMED_SEVERAL_POSSIBLE"
            elif names == [own]:
                verdict[rid] = "AGREES"
            elif own in names:
                verdict[rid] = "AGREES_SHARED_CLUSTER"
            else:
                verdict[rid] = "DISAGREES"

        # identification test (same rules as identification_test.py)
        titles = defaultdict(set)
        for rid, sid in anchors.items():
            tt = merged[rid]["publication"].strip().lower()
            if tt:
                titles[sid].add(tt)

        def study_of(rid):
            if rid in anchors:
                return anchors[rid]
            tt = merged[rid]["publication"].strip().lower()
            for sid, ts in titles.items():
                if tt and tt in ts:
                    return sid
            return ""
        study = {rid: study_of(rid) for rid in records}

        def cleaned(held_out):
            lv = defaultdict(lambda: defaultdict(set))
            for rid, sid in anchors.items():
                if sid != held_out and name[rid]:
                    lv[cl[rid]][studies[sid]["evidence_level"]].add(name[rid])
            keep = []
            for rid in records:
                if not name[rid]:
                    continue
                v = lv.get(cl[rid])
                if v and name[rid] not in v[min(v)]:
                    continue
                keep.append(rid)
            return keep

        def identify(q, db):
            hq = hap_of[q]
            ds = np.array([D[hq, hap_of[r]] for r in db]) if db else np.array([])
            ok = ~np.isnan(ds)
            if not ok.any():
                return "no match", ""
            dmin = ds[ok].min()
            if dmin > t:
                return "no match", ""
            hits = [db[i] for i in np.where(ok & (np.abs(ds - dmin) <= 1e-12))[0]]
            names = Counter(name[r] for r in hits if name[r])
            if not names:
                return "no match", ""
            if len(names) > 1:
                return "ambiguous", "; ".join(n for n, _ in names.most_common())
            top = next(iter(names))
            return ("correct" if top == name[q] else "wrong"), top

        rows, cache = [], {}
        for q, sid in sorted(anchors.items(), key=lambda kv: kv[1]):
            raw_db = [r for r in records if r != q and study[r] != sid]
            if sid not in cache:
                cache[sid] = cleaned(sid)
            cln_db = [r for r in cache[sid] if r != q and study[r] != sid]
            rr, rm = identify(q, raw_db)
            cr, cm = identify(q, cln_db)
            rows.append({"query": q, "study_id": sid, "evidence_level": studies[sid]["evidence_level"],
                         "species": name[q], "raw_result": rr, "raw_match": rm,
                         "cleaned_result": cr, "cleaned_match": cm})
        return {"clusters": n_cl, "mixed": mixed, "split": split, "cinnabarinus_clusters": cinna,
                "urticae_main_pct": urt_share, "verdict": verdict, "anchored": anchored,
                "cluster": cl, "queries": rows, "name": name}

    def acc(rows, key):
        c = Counter(r[key] for r in rows)
        return {o: c[o] for o in OUTCOMES}

    lines = []
    # 1. thresholds
    thr_rows, base = [], None
    for t in THRESHOLDS:
        res = run(t, "AB", base_syn)
        if abs(t - 0.03) < 1e-9:
            base = res
        n = len(res["queries"])
        a_raw, a_cln = acc(res["queries"], "raw_result"), acc(res["queries"], "cleaned_result")
        v = Counter(res["verdict"].values())
        nv = sum(v.values())
        thr_rows.append({
            "threshold": t, "clusters": res["clusters"], "clusters_several_names": res["mixed"],
            "species_split": res["split"], "cinnabarinus_clusters": res["cinnabarinus_clusters"],
            "urticae_in_largest_cluster_pct": f"{res['urticae_main_pct']:.1f}",
            "anchored_clusters": len(res["anchored"]),
            "records_disagree_pct": f"{100 * v['DISAGREES'] / nv:.1f}",
            "records_no_anchor_pct": f"{100 * v['NO_ANCHOR'] / nv:.1f}",
            "queries": n,
            **{f"raw_{o.replace(' ', '_')}": a_raw[o] for o in OUTCOMES},
            **{f"cleaned_{o.replace(' ', '_')}": a_cln[o] for o in OUTCOMES},
            "raw_correct_pct": f"{100 * a_raw['correct'] / n:.1f}",
            "cleaned_correct_pct": f"{100 * a_cln['correct'] / n:.1f}"})
    write(os.path.join(args.out, "thresholds.csv"), thr_rows)
    lines.append("1. Thresholds (clusters and largest accepted match distance):")
    lines.append("   K2P | clusters | several names | species split | cinnabarinus clusters | urticae largest % | "
                 "disagree % | correct raw -> cleaned (wrong raw -> cleaned)")
    for r in thr_rows:
        lines.append(f"   {r['threshold']:.3f} | {r['clusters']} | {r['clusters_several_names']} | "
                     f"{r['species_split']} | {r['cinnabarinus_clusters']} | {r['urticae_in_largest_cluster_pct']} | "
                     f"{r['records_disagree_pct']} | {r['raw_correct']}/{r['queries']} -> {r['cleaned_correct']}/"
                     f"{r['queries']} ({r['raw_wrong']} -> {r['cleaned_wrong']})")

    # 2. anchors and taxonomy
    ext_path = os.path.join(args.anchors, "synonyms_extended.csv")
    syn_plus = ({r["name"]: r["accepted_name"] for r in read(ext_path)} if os.path.exists(ext_path)
                else {**base_syn, **EXTRA_SYNONYMS})
    variants = [("main (levels A+B; 2 synonymies)", base),
                ("level A anchors only", run(0.03, "A", base_syn)),
                ("plus malaysiensis = macfarlanei", run(0.03, "AB", syn_plus)),
                ("level A only, plus malaysiensis = macfarlanei", run(0.03, "A", syn_plus))]
    var_rows = []
    lines.append("2. Anchor evidence and taxonomic treatment (K2P 0.03):")
    for label, res in variants:
        qs = res["queries"]
        n = len(qs)
        a_raw, a_cln = acc(qs, "raw_result"), acc(qs, "cleaned_result")
        v = Counter(res["verdict"].values())
        var_rows.append({"variant": label, "queries": n, "studies": len({q["study_id"] for q in qs}),
                         **{f"raw_{o.replace(' ', '_')}": a_raw[o] for o in OUTCOMES},
                         **{f"cleaned_{o.replace(' ', '_')}": a_cln[o] for o in OUTCOMES},
                         "anchored_clusters": len(res["anchored"]),
                         "records_disagree": v["DISAGREES"], "non_anchor_records": sum(v.values())})
        lines.append(f"   {label}: {n} queries from {var_rows[-1]['studies']} studies; raw "
                     + ", ".join(f"{o} {a_raw[o]}" for o in OUTCOMES) + "; cleaned "
                     + ", ".join(f"{o} {a_cln[o]}" for o in OUTCOMES)
                     + f"; disagreeing records {v['DISAGREES']} of {sum(v.values())}")
    write(os.path.join(args.out, "variants.csv"), var_rows)

    # 3. uncertainty of the main test
    qs = base["queries"]
    write(os.path.join(args.out, "queries_main.csv"), qs)
    trans = Counter((q["raw_result"], q["cleaned_result"]) for q in qs)
    trans_rows = [{"raw": a, **{b: trans[(a, b)] for b in OUTCOMES}, "total": sum(trans[(a, b)] for b in OUTCOMES)}
                  for a in OUTCOMES]
    trans_rows.append({"raw": "total", **{b: sum(trans[(a, b)] for a in OUTCOMES) for b in OUTCOMES},
                       "total": len(qs)})
    write(os.path.join(args.out, "transitions.csv"), trans_rows)
    gained = sum(1 for q in qs if q["raw_result"] != "correct" and q["cleaned_result"] == "correct")
    lost = sum(1 for q in qs if q["raw_result"] == "correct" and q["cleaned_result"] != "correct")
    p_mc = binom_two_sided(min(gained, lost), gained + lost)

    by_study = defaultdict(list)
    for q in qs:
        by_study[q["study_id"]].append(q)
    sids = sorted(by_study)
    rng = np.random.default_rng(args.seed)
    stat = defaultdict(list)
    for _ in range(args.boot):
        pick = [q for s in rng.choice(sids, size=len(sids), replace=True) for q in by_study[s]]
        n = len(pick)
        for key in ("raw", "cleaned"):
            for o in ("correct", "wrong"):
                stat[f"{key}_{o}"].append(100 * sum(q[f"{key}_result"] == o for q in pick) / n)
        stat["diff_correct"].append(stat["cleaned_correct"][-1] - stat["raw_correct"][-1])
        stat["diff_wrong"].append(stat["cleaned_wrong"][-1] - stat["raw_wrong"][-1])
    n = len(qs)
    point = {f"{k}_{o}": 100 * sum(q[f"{k}_result"] == o for q in qs) / n
             for k in ("raw", "cleaned") for o in ("correct", "wrong")}
    point["diff_correct"] = point["cleaned_correct"] - point["raw_correct"]
    point["diff_wrong"] = point["cleaned_wrong"] - point["raw_wrong"]
    ci_rows = [{"measure": k, "estimate_pct": f"{point[k]:.1f}",
                "ci_low": f"{np.percentile(v, 2.5):.1f}", "ci_high": f"{np.percentile(v, 97.5):.1f}"}
               for k, v in stat.items()]
    write(os.path.join(args.out, "identification_intervals.csv"), ci_rows)
    lines.append("3. Main identification test: outcomes raw (rows) -> cleaned (columns):")
    lines.append("   " + " | ".join(["raw \\ cleaned"] + OUTCOMES + ["total"]))
    for r in trans_rows:
        lines.append("   " + " | ".join(str(r[k]) for k in ["raw"] + OUTCOMES + ["total"]))
    lines.append(f"   Became correct {gained}, stopped being correct {lost}; exact McNemar p = {p_mc:.3f}")
    lines.append(f"   95% intervals, bootstrap over {len(sids)} anchor studies ({args.boot} replicates, seed {args.seed}):")
    for r in ci_rows:
        lines.append(f"     {r['measure']}: {r['estimate_pct']}% ({r['ci_low']} to {r['ci_high']})")

    sp_rows = []
    per_sp = defaultdict(list)
    for q in qs:
        per_sp[q["species"]].append(q)
    for sp, rs in sorted(per_sp.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        sp_rows.append({"species": sp, "queries": len(rs), "studies": len({r["study_id"] for r in rs}),
                        **{f"raw_{o.replace(' ', '_')}": sum(r["raw_result"] == o for r in rs) for o in OUTCOMES},
                        **{f"cleaned_{o.replace(' ', '_')}": sum(r["cleaned_result"] == o for r in rs)
                           for o in OUTCOMES},
                        "wrong_names_cleaned": "; ".join(f"{k} ({v})" for k, v in Counter(
                            r["cleaned_match"] for r in rs if r["cleaned_result"] == "wrong").most_common())})
    write(os.path.join(args.out, "identification_by_species.csv"), sp_rows)

    # 4. cluster structure at 0.03
    cl = base["cluster"]
    members = defaultdict(set)
    for rid in records:
        members[cl[rid]].add(hap_of[rid])
    diam_rows = []
    cov5 = np.array([sum(c != "-" for c in s[:700]) for s in hap_seqs])
    cov3 = np.array([sum(c != "-" for c in s[700:]) for s in hap_seqs])
    bridge = (cov5 >= args.min_overlap) & (cov3 >= args.min_overlap)
    only5 = (cov5 >= args.min_overlap) & ~bridge
    only3 = (cov3 >= args.min_overlap) & ~bridge
    cid_names = {}
    for rid in records:
        cid_names.setdefault(cl[rid], Counter())[base["name"][rid] or "(unnamed)"] += 1
    joined_only_by_bridges = 0
    for c, hs in members.items():
        hs = sorted(hs)
        sub = D[np.ix_(hs, hs)]
        fin = sub[np.isfinite(sub)]
        diam = float(fin.max()) if fin.size else 0.0
        has5 = any(only5[h] for h in hs)
        has3 = any(only3[h] for h in hs)
        split_without_bridges = ""
        if has5 and has3:
            keep = [h for h in hs if not bridge[h]]
            uf = UnionFind(len(keep))
            for a in range(len(keep)):
                for b in range(a + 1, len(keep)):
                    if Dfin[keep[a], keep[b]] <= 0.03:
                        uf.union(a, b)
            comps = defaultdict(set)
            for a, h in enumerate(keep):
                comps[uf.find(a)].add(h)
            mixed_comp = any(any(only5[h] for h in s) and any(only3[h] for h in s) for s in comps.values())
            split_without_bridges = "no" if mixed_comp else "yes"
            joined_only_by_bridges += split_without_bridges == "yes"
        n_rec = sum(1 for rid in records if cl[rid] == c)
        diam_rows.append({"cluster": checks[next(r for r in records if cl[r] == c)]["cluster"],
                          "records": n_rec, "haplotypes": len(hs),
                          "largest_k2p_within": f"{diam:.4f}",
                          "has_5prime_only": "yes" if has5 else "no", "has_3prime_only": "yes" if has3 else "no",
                          "bridging_haplotypes": int(sum(bridge[h] for h in hs)),
                          "halves_joined_only_through_bridges": split_without_bridges,
                          "composition": "; ".join(f"{k} ({v})" for k, v in cid_names[c].most_common())})
    diam_rows.sort(key=lambda r: -r["records"])
    write(os.path.join(args.out, "cluster_structure.csv"), diam_rows)
    wide = [r for r in diam_rows if float(r["largest_k2p_within"]) > 0.03]
    both_halves = [r for r in diam_rows if r["has_5prime_only"] == "yes" and r["has_3prime_only"] == "yes"]
    n_bridge_rec = sum(1 for rid in records if bridge[hap_of[rid]])
    lines.append("4. Cluster structure at K2P 0.03:")
    lines.append(f"   Clusters whose largest internal K2P exceeds 0.03 (chaining): {len(wide)} of {len(diam_rows)}")
    for r in wide[:10]:
        lines.append(f"     {r['cluster']} ({r['records']} records): largest {r['largest_k2p_within']}")
    lines.append(f"   Records overlapping both halves by >= {args.min_overlap} sites (bridges): {n_bridge_rec}")
    lines.append(f"   Clusters holding both 5'-only and 3'-only records: {len(both_halves)}; "
                 f"of these, joined only through bridges: {joined_only_by_bridges}")
    for r in both_halves:
        lines.append(f"     {r['cluster']} ({r['records']} records, {r['bridging_haplotypes']} bridging haplotypes, "
                     f"joined only through bridges: {r['halves_joined_only_through_bridges']})")

    # 5. nearest-neighbour label agreement (comparable with ITS2)
    name = base["name"]
    named = [rid for rid in records if name[rid]]
    nn_rows = []
    rec_hap = np.array([hap_of[r] for r in records])
    rec_batch = np.array([batch.get(r, "") for r in records])
    idx = np.arange(len(records))
    for qi, q in enumerate(records):
        hq = hap_of[q]
        out = {"record_id": q, "accepted_name": name[q], "batch": batch.get(q, "")}
        d_all = D[hq, rec_hap]
        for label, mask in (("any", idx != qi),
                            ("other_batch", (idx != qi) & (rec_batch != rec_batch[qi]))):
            ok = mask & ~np.isnan(d_all)
            nm, dist = "", ""
            if ok.any():
                dmin = d_all[ok].min()
                hits = [records[i] for i in np.where(ok & (np.abs(d_all - dmin) <= 1e-12))[0]]
                c = Counter(name[r] for r in hits if name[r])
                if c:
                    nm, dist = c.most_common(1)[0][0], f"{max(float(dmin), 0.0):.4f}"
            out[f"nn_name_{label}"] = nm
            out[f"nn_distance_{label}"] = dist
            out[f"nn_match_{label}"] = ("yes" if nm == name[q] else "no") if (nm and name[q]) else ""
        nn_rows.append(out)
    write(os.path.join(args.out, "coi_nearest_neighbour.csv"), nn_rows)
    nn_by = {r["record_id"]: r for r in nn_rows}

    def agree(rows, key):
        v = [r[key] for r in rows if r[key]]
        return sum(x == "yes" for x in v), len(v)

    lines.append("5. Nearest-neighbour label agreement:")
    for label in ("any", "other_batch"):
        k, m = agree(nn_rows, f"nn_match_{label}")
        lines.append(f"   COI, all named records, {label.replace('_', ' ')}: {k} of {m} ({100 * k / m:.1f}%)")
    its_path = os.path.join(args.its2, "its2_records.csv")
    pair_path = os.path.join(args.its2, "coi_its2_pairs.csv")
    if os.path.exists(its_path) and os.path.exists(pair_path):
        its = {r["record_id"]: r for r in read(its_path)}
        for label in ("any", "other_batch"):
            key = "nn_match" if label == "any" else "nn_match_other_batch"
            if key in next(iter(its.values())):
                k, m = agree(its.values(), key)
                lines.append(f"   ITS2, all named records, {label.replace('_', ' ')}: {k} of {m} ({100 * k / m:.1f}%)")
        pairs = read(pair_path)
        pair_rows = []
        for pr in pairs:
            c, i = nn_by.get(pr["coi_record"]), its.get(pr["its2_record"])
            if not c or not i:
                continue
            pair_rows.append({"specimen": pr["specimen"], "coi_record": pr["coi_record"],
                              "its2_record": pr["its2_record"], "label": c["accepted_name"] or i["accepted_name"],
                              "coi_nn_other_batch": c["nn_name_other_batch"],
                              "its2_nn_other_batch": i.get("nn_name_other_batch", ""),
                              "coi_match": c["nn_match_other_batch"],
                              "its2_match": i.get("nn_match_other_batch", "")})
        write(os.path.join(args.out, "paired_label_agreement.csv"), pair_rows)
        both = [r for r in pair_rows if r["coi_match"] and r["its2_match"]]
        if both:
            kc = sum(r["coi_match"] == "yes" for r in both)
            ki = sum(r["its2_match"] == "yes" for r in both)
            only_c = sum(r["coi_match"] == "yes" and r["its2_match"] == "no" for r in both)
            only_i = sum(r["coi_match"] == "no" and r["its2_match"] == "yes" for r in both)
            lines.append(f"   Same specimens ({len(both)} with both genes, other batches only): COI {kc} "
                         f"({100 * kc / len(both):.1f}%), ITS2 {ki} ({100 * ki / len(both):.1f}%); "
                         f"COI only {only_c}, ITS2 only {only_i}, exact McNemar p = "
                         f"{binom_two_sided(min(only_c, only_i), only_c + only_i):.3f}")
            per = defaultdict(lambda: [0, 0, 0])
            for r in both:
                per[r["label"]][0] += 1
                per[r["label"]][1] += r["coi_match"] == "yes"
                per[r["label"]][2] += r["its2_match"] == "yes"
            for sp, (m, kc_, ki_) in sorted(per.items(), key=lambda kv: -kv[1][0]):
                lines.append(f"     {sp}: {m} specimens; COI {kc_}, ITS2 {ki_}")

    # 6. disagreements
    dis = [rid for rid, v in base["verdict"].items() if v == "DISAGREES"]
    pair_c = Counter((name[r], "; ".join(base["anchored"][cl[r]])) for r in dis)
    dis_rows = []
    for (own, anc), k in pair_c.most_common():
        rs = [r for r in dis if name[r] == own and "; ".join(base["anchored"][cl[r]]) == anc]
        dis_rows.append({"record_name": own, "submitted_names": "; ".join(sorted({raw_name[r] for r in rs})),
                         "cluster_anchored_species": anc, "records": k,
                         "batches": len({batch.get(r) for r in rs}),
                         "countries": "; ".join(f"{c} ({v})" for c, v in Counter(
                             merged[r]["country"].split(":")[0].strip() or "(not given)" for r in rs).most_common())})
    write(os.path.join(args.out, "disagreements.csv"), dis_rows)
    lines.append(f"6. Records disagreeing with the anchored species of their cluster: {len(dis)} "
                 f"from {len({batch.get(r) for r in dis})} submission batches")
    for r in dis_rows:
        lines.append(f"   {r['records']} x {r['submitted_names']} in cluster anchored as {r['cluster_anchored_species']}"
                     f" ({r['batches']} batches; {r['countries']})")

    summary = "\n".join(lines)
    with open(os.path.join(args.out, "robustness_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
