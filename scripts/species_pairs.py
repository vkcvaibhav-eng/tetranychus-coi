"""Module 12: what kind of COI overlap lies behind each pair of species names.

A single-linkage cluster holding two names does not show that the two
species share mitochondrial haplotypes: chaining can join records much
further apart than the threshold. For every pair of accepted names that share
a cluster, or come within the cluster threshold of each other, this module
reports, over sequences that share at least the minimum number of sites:

  identical         at least one pair of sequences at K2P 0
  close             smallest K2P above 0 and within the threshold
  chained only      in one cluster, but no pair within the threshold
  separate          not in one cluster (listed only when within the threshold)

once for all records under each name (database labels) and once for anchor
records only (identifications from the anchor studies, levels A and B, names
as given in the study), with the largest distance within each species for
comparison.

Reads:  as identification_test.py
Writes: <out>/species_pairs.csv, <out>/species_pairs_summary.txt

Usage: python scripts/species_pairs.py [--genus results/genus --species results/species --anchors data/anchors]
"""
import argparse
import csv
import itertools
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from identification_test import anchors_at, cluster_ids, load_context  # noqa: E402


def main():
    p = argparse.ArgumentParser(description="COI overlap behind pairs of species names")
    p.add_argument("--genus", default=os.path.join("results", "genus"))
    p.add_argument("--species", default=os.path.join("results", "species"))
    p.add_argument("--anchors", default=os.path.join("data", "anchors"))
    p.add_argument("--reference", default=os.path.join("data", "reference", "coi_reference.fasta"))
    p.add_argument("--data", default="data")
    p.add_argument("--out", default=os.path.join("results", "pairs"))
    p.add_argument("--threshold", type=float, default=0.03)
    p.add_argument("--min-overlap", type=int, default=300)
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ctx = load_context(args.genus, args.species, args.anchors, args.reference, args.data, args.min_overlap)
    syn, D, hap = ctx["synonyms"], ctx["D"], ctx["hap_of"]
    records = ctx["records"]
    name = {r: syn.get(n, n) for r, n in ctx["raw_name"].items()}
    anchors = anchors_at(ctx, "AB")
    a_name = {r: syn.get(ctx["anchor_name"].get(r, name[r]), ctx["anchor_name"].get(r, name[r])) for r in anchors}
    cl = cluster_ids(ctx, args.threshold)

    by_name, by_anchor = defaultdict(list), defaultdict(list)
    for r in records:
        if name[r]:
            by_name[name[r]].append(r)
    for r, n in a_name.items():
        if n:
            by_anchor[n].append(r)

    def haps(rs):
        return sorted({hap[r] for r in rs})

    def min_between(ra, rb):
        sub = D[np.ix_(haps(ra), haps(rb))]
        return float(np.nanmin(sub)) if np.isfinite(sub).any() else np.nan

    def max_within(rs):
        h = haps(rs)
        if len(h) < 2:
            return 0.0
        sub = D[np.ix_(h, h)]
        return float(np.nanmax(sub)) if np.isfinite(sub).any() else np.nan

    def kind(dmin, same_cluster):
        if np.isnan(dmin):
            return "chained only (no overlapping pair)" if same_cluster else "no overlapping pair"
        if dmin <= 1e-9:
            return "identical"
        if dmin <= args.threshold:
            return "close" if same_cluster else "close, separate clusters"
        return "chained only" if same_cluster else "separate"

    clusters_of = {n: {cl[r] for r in rs} for n, rs in by_name.items()}
    a_clusters = {n: {cl[r] for r in rs} for n, rs in by_anchor.items()}
    rows = []
    for a, b in itertools.combinations(sorted(set(by_name) | set(by_anchor)), 2):
        both_db = a in by_name and b in by_name
        both_an = a in by_anchor and b in by_anchor
        shared = clusters_of.get(a, set()) & clusters_of.get(b, set())
        ash = a_clusters.get(a, set()) & a_clusters.get(b, set())
        dmin = min_between(by_name[a], by_name[b]) if both_db else np.nan
        ad = min_between(by_anchor[a], by_anchor[b]) if both_an else np.nan
        if not (shared or ash or dmin <= args.threshold or ad <= args.threshold):
            continue
        rows.append({"species_a": a, "species_b": b, "records_a": len(by_name.get(a, [])),
                     "records_b": len(by_name.get(b, [])), "shared_clusters": len(shared),
                     "min_k2p_all_records": f"{abs(dmin):.4f}" if not np.isnan(dmin) else "",
                     "kind_all_records": kind(dmin, bool(shared)) if both_db else "not both named in the database",
                     "max_within_a": f"{max_within(by_name[a]):.4f}" if a in by_name else "",
                     "max_within_b": f"{max_within(by_name[b]):.4f}" if b in by_name else "",
                     "anchors_a": len(by_anchor.get(a, [])), "anchors_b": len(by_anchor.get(b, [])),
                     "min_k2p_anchors": f"{abs(ad):.4f}" if not np.isnan(ad) else "",
                     "kind_anchors": kind(ad, bool(ash)) if both_an else "no anchors of both"})
    rows.sort(key=lambda r: (r["kind_anchors"] == "no anchors of both", r["species_a"], r["species_b"]))
    with open(os.path.join(args.out, "species_pairs.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    lines = [f"Pairs of accepted names sharing a cluster (K2P {args.threshold}) or within it: {len(rows)}",
             "species pair | all records: kind (min K2P) | anchors: kind (min K2P)"]
    for r in rows:
        lines.append(f"  {r['species_a']} / {r['species_b']} | {r['kind_all_records']} ({r['min_k2p_all_records']}) | "
                     f"{r['kind_anchors']} ({r['min_k2p_anchors']})")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "species_pairs_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
