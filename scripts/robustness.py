"""Module 10: sensitivity and uncertainty checks of the COI results.

1. Thresholds. Clusters, anchor verdicts and the identification test are
   repeated at several K2P thresholds, first with the cluster threshold and
   the largest accepted match distance equal, then varied separately.
2. Anchors, taxonomy and fold isolation. At K2P 0.03 the identification
   test is repeated with level A anchors only, with T. malaysiensis treated
   as a junior synonym of T. macfarlanei (Zeity et al. 2017), and in the
   earlier, less strict form (clusters from all records, study membership by
   publication title only) to measure any leakage from held-out records.
3. Uncertainty. For the main identification test: the full table of
   outcomes raw -> cleaned, an exact McNemar test on the queries whose
   correctness changed, and 95% intervals from a bootstrap that resamples
   anchor studies, since anchors of one study are not independent.
4. Cluster structure. Largest K2P distance within each cluster (single
   linkage can chain records that are far apart), and the clusters whose 5'
   and 3' fragments are joined only through records long enough to overlap
   both.
5. Label agreement comparable with ITS2. For every named COI record, the
   name of its nearest other record, over all records and over records from
   other submission batches only. For specimens linked one-to-one to an ITS2
   record (its2_analysis.py), the two markers are compared on the same
   specimens, with a bootstrap over submission batches and a sign-flip test
   on batch totals, because specimens of one batch are not independent.
6. Records that disagree with the anchored species of their cluster, by
   pair of names and by submission batch.

Reads:  results/genus/*, results/species/record_checks.csv, data/anchors/*,
        data/ncbi_*.gb, reference COX1, results/its2/* (optional)
Writes: results/robustness/*.csv, results/robustness/robustness_summary.txt

Usage: python scripts/robustness.py [--boot 10000]
"""
import argparse
import csv
import itertools
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from identification_test import OUTCOMES, cluster_ids, evaluate, load_context, tally  # noqa: E402
from species_analysis import UnionFind  # noqa: E402

THRESHOLDS = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05]
MATCH_DISTANCES = [0.01, 0.02, 0.03, 0.04, 0.05]
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


def sign_flip_p(totals, draws=200000, seed=1):
    """Two-sided sign-flip test of sum(totals) = 0, flipping whole groups."""
    t = np.array([x for x in totals if x != 0], dtype=float)
    if t.size == 0:
        return 1.0
    obs = abs(t.sum())
    if t.size <= 20:
        sums = [abs(np.dot(signs, t)) for signs in itertools.product((1, -1), repeat=t.size)]
        return float(np.mean(np.array(sums) >= obs - 1e-9))
    rng = np.random.default_rng(seed)
    signs = rng.choice((1, -1), size=(draws, t.size))
    return float((np.sum(np.abs(signs @ t) >= obs - 1e-9) + 1) / (draws + 1))


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

    ctx = load_context(args.genus, args.species, args.anchors, args.reference, args.data, args.min_overlap)
    records, D, Dfin, hap_of = ctx["records"], ctx["D"], ctx["Dfin"], ctx["hap_of"]
    merged, studies, anchor_all, batch = ctx["merged"], ctx["studies"], ctx["anchor_all"], ctx["batch"]
    raw_name, base_syn = ctx["raw_name"], ctx["synonyms"]
    ext_path = os.path.join(args.anchors, "synonyms_extended.csv")
    syn_plus = ({r["name"]: r["accepted_name"] for r in read(ext_path)} if os.path.exists(ext_path)
                else {**base_syn, **EXTRA_SYNONYMS})
    print(f"{len(records)} records, {len(ctx['hap_seqs'])} haplotypes")

    def cluster_stats(t, levels="AB", syn=None):
        """Clusters and anchor verdicts over all records at threshold t."""
        syn = base_syn if syn is None else syn
        name = {rid: syn.get(raw_name[rid], raw_name[rid]) for rid in records}
        study_name = {rid: syn.get(n, n) for rid, n in ctx.get("anchor_name", {}).items()}
        cl = cluster_ids(ctx, t)
        comp = defaultdict(Counter)
        for rid in records:
            comp[cl[rid]][raw_name[rid] or "(unnamed)"] += 1
        sp_clusters = defaultdict(set)
        for rid in records:
            if raw_name[rid]:
                sp_clusters[raw_name[rid]].add(cl[rid])
        urt = Counter(cl[r] for r in records if raw_name[r] == "Tetranychus urticae")
        anchors = {rid: sid for rid, sid in anchor_all.items() if studies[sid]["evidence_level"] in levels}
        best = defaultdict(lambda: defaultdict(set))
        for rid, sid in anchors.items():
            n = study_name.get(rid, name[rid])
            if n:
                best[cl[rid]][studies[sid]["evidence_level"]].add(n)
        anchored = {c: sorted(v[min(v)]) for c, v in best.items()}
        verdict = {}
        for rid in records:
            if rid in anchor_all:  # anchors of every level are left out, as in anchor_check.py
                continue
            names, own = anchored.get(cl[rid], []), name[rid]
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
        return {"clusters": len(comp),
                "mixed": sum(1 for c in comp.values() if len([n for n in c if n != "(unnamed)"]) > 1),
                "split": sum(1 for v in sp_clusters.values() if len(v) > 1),
                "cinnabarinus_clusters": len(sp_clusters.get("Tetranychus cinnabarinus", ())),
                "urticae_main_pct": 100 * urt.most_common(1)[0][1] / sum(urt.values()) if urt else float("nan"),
                "verdict": verdict, "anchored": anchored, "cluster": cl, "name": name}

    def outcome_cols(qs):
        a, b = tally(qs, "raw_result"), tally(qs, "cleaned_result")
        return {**{f"raw_{o.replace(' ', '_')}": a[o] for o in OUTCOMES},
                **{f"cleaned_{o.replace(' ', '_')}": b[o] for o in OUTCOMES}}

    lines = []
    # 1a. thresholds, cluster and match distance equal
    thr_rows, base, base_q = [], None, None
    for t in THRESHOLDS:
        st = cluster_stats(t)
        qs = evaluate(ctx, t, t)
        if abs(t - 0.03) < 1e-9:
            base, base_q = st, qs
        v = Counter(st["verdict"].values())
        nv = sum(v.values())
        n = len(qs)
        oc = outcome_cols(qs)
        thr_rows.append({"threshold": t, "clusters": st["clusters"], "clusters_several_names": st["mixed"],
                         "species_split": st["split"], "cinnabarinus_clusters": st["cinnabarinus_clusters"],
                         "urticae_in_largest_cluster_pct": f"{st['urticae_main_pct']:.1f}",
                         "anchored_clusters": len(st["anchored"]),
                         "records_disagree_pct": f"{100 * v['DISAGREES'] / nv:.1f}",
                         "records_no_anchor_pct": f"{100 * v['NO_ANCHOR'] / nv:.1f}",
                         "queries": n, **oc,
                         "raw_correct_pct": f"{100 * oc['raw_correct'] / n:.1f}",
                         "cleaned_correct_pct": f"{100 * oc['cleaned_correct'] / n:.1f}"})
    write(os.path.join(args.out, "thresholds.csv"), thr_rows)
    lines.append("1a. Thresholds, cluster threshold = largest match distance:")
    tetra = any(raw_name[r] == "Tetranychus urticae" for r in records)
    lines.append("   K2P | clusters | several names | species split | "
                 + ("cinnabarinus clusters | urticae largest % | " if tetra else "")
                 + "disagree % | correct raw -> cleaned (wrong raw -> cleaned)")
    for r in thr_rows:
        extra = f"{r['cinnabarinus_clusters']} | {r['urticae_in_largest_cluster_pct']} | " if tetra else ""
        lines.append(f"   {r['threshold']:.3f} | {r['clusters']} | {r['clusters_several_names']} | "
                     f"{r['species_split']} | {extra}"
                     f"{r['records_disagree_pct']} | {r['raw_correct']}/{r['queries']} -> {r['cleaned_correct']}/"
                     f"{r['queries']} ({r['raw_wrong']} -> {r['cleaned_wrong']})")

    # 1b. cluster threshold and match distance varied separately
    grid_rows = []
    for tc in THRESHOLDS:
        for tm in MATCH_DISTANCES:
            qs = base_q if (abs(tc - 0.03) < 1e-9 and abs(tm - 0.03) < 1e-9) else evaluate(ctx, tc, tm)
            grid_rows.append({"cluster_threshold": tc, "match_distance": tm, "queries": len(qs), **outcome_cols(qs)})
    write(os.path.join(args.out, "threshold_grid.csv"), grid_rows)
    lines.append("1b. Correct answers, cleaned database (raw in brackets), by cluster threshold (rows) "
                 "and largest match distance (columns):")
    lines.append("   cluster \\ match | " + " | ".join(f"{tm:.2f}" for tm in MATCH_DISTANCES))
    for tc in THRESHOLDS:
        cells = [r for r in grid_rows if r["cluster_threshold"] == tc]
        lines.append(f"   {tc:.3f} | " + " | ".join(f"{r['cleaned_correct']} ({r['raw_correct']})" for r in cells))

    # 2. anchors, taxonomy, isolation
    n_syn = {0: "no synonymies", 1: "one synonymy", 2: "two synonymies"}.get(len(base_syn), f"{len(base_syn)} synonymies")
    added = sorted(set(syn_plus.items()) - set(base_syn.items()))
    plus_label = "plus " + "; ".join(f"{a.split()[-1]} = {b.split()[-1]}" for a, b in added) if added else "extended synonymy (no change)"
    variants = [(f"main (levels A+B; {n_syn}; folds isolated)", base_q, base),
                ("level A anchors only", evaluate(ctx, 0.03, 0.03, "A"), cluster_stats(0.03, "A")),
                (plus_label, evaluate(ctx, 0.03, 0.03, "AB", syn_plus),
                 cluster_stats(0.03, "AB", syn_plus)),
                ("earlier version (clusters from all records; title membership)",
                 evaluate(ctx, 0.03, 0.03, isolate=False, rule="title"), base)]
    var_rows = []
    lines.append("2. Anchor evidence, taxonomic treatment and fold isolation (K2P 0.03):")
    for label, qs, st in variants:
        oc = outcome_cols(qs)
        v = Counter(st["verdict"].values())
        changed = sum(1 for a, b in zip(qs, base_q)
                      if (a["raw_result"], a["cleaned_result"]) != (b["raw_result"], b["cleaned_result"])) \
            if len(qs) == len(base_q) else ""
        var_rows.append({"variant": label, "queries": len(qs), "studies": len({q["study_id"] for q in qs}), **oc,
                         "queries_changed_vs_main": changed, "anchored_clusters": len(st["anchored"]),
                         "records_disagree": v["DISAGREES"], "non_anchor_records": sum(v.values())})
        lines.append(f"   {label}: {len(qs)} queries from {var_rows[-1]['studies']} studies; raw "
                     + ", ".join(f"{o} {oc['raw_' + o.replace(' ', '_')]}" for o in OUTCOMES) + "; cleaned "
                     + ", ".join(f"{o} {oc['cleaned_' + o.replace(' ', '_')]}" for o in OUTCOMES)
                     + (f"; queries changed vs main {changed}" if changed != "" else ""))
    write(os.path.join(args.out, "variants.csv"), var_rows)

    # 3. uncertainty of the main test
    qs = base_q
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
                        **outcome_cols(rs),
                        "wrong_names_cleaned": "; ".join(f"{k} ({v})" for k, v in Counter(
                            r["cleaned_match"] for r in rs if r["cleaned_result"] == "wrong").most_common())})
    write(os.path.join(args.out, "identification_by_species.csv"), sp_rows)

    # 4. cluster structure at 0.03
    cl = base["cluster"]
    hap_seqs = ctx["hap_seqs"]
    members = defaultdict(set)
    for rid in records:
        members[cl[rid]].add(hap_of[rid])
    cov5 = np.array([sum(c != "-" for c in s[:700]) for s in hap_seqs])
    cov3 = np.array([sum(c != "-" for c in s[700:]) for s in hap_seqs])
    bridge = (cov5 >= args.min_overlap) & (cov3 >= args.min_overlap)
    only5 = (cov5 >= args.min_overlap) & ~bridge
    only3 = (cov3 >= args.min_overlap) & ~bridge
    cid_names = defaultdict(Counter)
    for rid in records:
        cid_names[cl[rid]][base["name"][rid] or "(unnamed)"] += 1
    checks = ctx["checks"]
    diam_rows, joined_only_by_bridges = [], 0
    for c, hs in members.items():
        hs = sorted(hs)
        sub = D[np.ix_(hs, hs)]
        fin = sub[np.isfinite(sub)]
        diam = float(fin.max()) if fin.size else 0.0
        has5, has3 = any(only5[h] for h in hs), any(only3[h] for h in hs)
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
        in_c = [r for r in records if cl[r] == c]
        diam_rows.append({"cluster": checks[in_c[0]]["cluster"], "records": len(in_c), "haplotypes": len(hs),
                          "largest_k2p_within": f"{diam:.4f}",
                          "has_5prime_only": "yes" if has5 else "no", "has_3prime_only": "yes" if has3 else "no",
                          "bridging_haplotypes": int(sum(bridge[h] for h in hs)),
                          "halves_joined_only_through_bridges": split_without_bridges,
                          "composition": "; ".join(f"{k} ({v})" for k, v in cid_names[c].most_common())})
    diam_rows.sort(key=lambda r: -r["records"])
    write(os.path.join(args.out, "cluster_structure.csv"), diam_rows)
    wide = [r for r in diam_rows if float(r["largest_k2p_within"]) > 0.03]
    both_halves = [r for r in diam_rows if r["has_5prime_only"] == "yes" and r["has_3prime_only"] == "yes"]
    lines.append("4. Cluster structure at K2P 0.03:")
    lines.append(f"   Clusters whose largest internal K2P exceeds 0.03 (chaining): {len(wide)} of {len(diam_rows)}")
    for r in wide[:10]:
        lines.append(f"     {r['cluster']} ({r['records']} records): largest {r['largest_k2p_within']}")
    lines.append(f"   Records overlapping both halves by >= {args.min_overlap} sites (bridges): "
                 f"{sum(1 for rid in records if bridge[hap_of[rid]])}")
    lines.append(f"   Clusters holding both 5'-only and 3'-only records: {len(both_halves)}; "
                 f"of these, joined only through bridges: {joined_only_by_bridges}")

    # 5. nearest-neighbour label agreement (comparable with ITS2)
    name = base["name"]
    rec_hap, rec_batch = ctx["rec_hap"], np.array([batch.get(r, "") for r in records])
    idx = np.arange(len(records))
    nn_rows = []
    for qi, q in enumerate(records):
        out = {"record_id": q, "accepted_name": name[q], "batch": batch.get(q, "")}
        d_all = D[hap_of[q], rec_hap]
        for label, mask in (("any", idx != qi), ("other_batch", (idx != qi) & (rec_batch != rec_batch[qi]))):
            ok = mask & ~np.isnan(d_all)
            nm, dist = "", ""
            if ok.any():
                dmin = d_all[ok].min()
                hits = [records[i] for i in np.where(ok & (np.abs(d_all - dmin) <= 1e-12))[0]]
                c = Counter(name[r] for r in hits if name[r])
                if c:
                    nm, dist = c.most_common(1)[0][0], f"{abs(max(float(dmin), 0.0)):.4f}"
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
        for label, key in (("any", "nn_match"), ("other batch", "nn_match_other_batch")):
            if key in next(iter(its.values())):
                k, m = agree(its.values(), key)
                lines.append(f"   ITS2, all named records, {label}: {k} of {m} ({100 * k / m:.1f}%)")
        pairs = read(pair_path)
        pair_rows = []
        for pr in pairs:
            c, i = nn_by.get(pr["coi_record"]), its.get(pr["its2_record"])
            if not c or not i:
                continue
            pair_rows.append({"specimen": pr["specimen"], "coi_record": pr["coi_record"],
                              "its2_record": pr["its2_record"], "link": pr.get("link", ""),
                              "confident": pr.get("confident", "yes"),
                              "label": c["accepted_name"] or i["accepted_name"],
                              "coi_batch": c["batch"],
                              "coi_nn_other_batch": c["nn_name_other_batch"],
                              "its2_nn_other_batch": i.get("nn_name_other_batch", ""),
                              "coi_match": c["nn_match_other_batch"], "its2_match": i.get("nn_match_other_batch", "")})
        write(os.path.join(args.out, "paired_label_agreement.csv"), pair_rows)
        lines.append(f"   Specimen links: {len(pair_rows)} pairs, {sum(r['confident'] == 'yes' for r in pair_rows)} "
                     f"one-to-one without conflict")
        both = [r for r in pair_rows if r["confident"] == "yes" and r["coi_match"] and r["its2_match"]]
        if both:
            kc = sum(r["coi_match"] == "yes" for r in both)
            ki = sum(r["its2_match"] == "yes" for r in both)
            only_c = sum(r["coi_match"] == "yes" and r["its2_match"] == "no" for r in both)
            only_i = sum(r["coi_match"] == "no" and r["its2_match"] == "yes" for r in both)
            groups = defaultdict(list)
            for r in both:
                groups[r["coi_batch"]].append(r)
            totals = [sum((r["coi_match"] == "yes") - (r["its2_match"] == "yes") for r in g) for g in groups.values()]
            keys = list(groups)
            rngp = np.random.default_rng(args.seed)
            diffs = []
            for _ in range(args.boot):
                pick = [r for k_ in rngp.choice(len(keys), size=len(keys), replace=True) for r in groups[keys[k_]]]
                diffs.append(100 * (sum(r["coi_match"] == "yes" for r in pick)
                                    - sum(r["its2_match"] == "yes" for r in pick)) / len(pick))
            lo, hi = np.percentile(diffs, [2.5, 97.5])
            p_flip = sign_flip_p(totals)
            lines.append(f"   Same specimens ({len(both)} one-to-one, other batches only, {len(groups)} COI submission "
                         f"batches): COI {kc} ({100 * kc / len(both):.1f}%), ITS2 {ki} ({100 * ki / len(both):.1f}%); "
                         f"COI only {only_c}, ITS2 only {only_i}")
            lines.append(f"     difference {100 * (kc - ki) / len(both):.1f} points; 95% interval over batches "
                         f"{lo:.1f} to {hi:.1f}; sign-flip test on batch totals p = {p_flip:.3f}; "
                         f"McNemar ignoring batches p = {binom_two_sided(min(only_c, only_i), only_c + only_i):.3f}")
            per = defaultdict(lambda: [0, 0, 0])
            for r in both:
                per[r["label"]][0] += 1
                per[r["label"]][1] += r["coi_match"] == "yes"
                per[r["label"]][2] += r["its2_match"] == "yes"
            for sp, (m, kc_, ki_) in sorted(per.items(), key=lambda kv: -kv[1][0]):
                lines.append(f"     {sp}: {m} specimens; COI {kc_}, ITS2 {ki_}")
            write(os.path.join(args.out, "paired_by_batch.csv"),
                  [{"coi_batch": k_, "specimens": len(g), "coi_agree": sum(r["coi_match"] == "yes" for r in g),
                    "its2_agree": sum(r["its2_match"] == "yes" for r in g)} for k_, g in groups.items()])

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
