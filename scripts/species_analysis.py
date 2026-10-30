"""Module 4: species-level analysis of quality-checked COI sequences.

Every sequence that passed the quality check is placed on the reference
COX1 coordinates (insertions relative to the reference are dropped), so
sequences from the 5' and 3' halves of the gene are only compared where
they actually overlap. On that basis the module computes:

- K2P distances between unique haplotypes over shared sites (minimum
  overlap set by --min-overlap);
- the barcode gap of each named species (largest distance within the
  species against the smallest distance to any other named species);
- haplotypes carried by records under more than one species name;
- records whose nearest neighbour belongs to a different named species;
- single-linkage clusters at --threshold and the names in each.

Reads:  <results>/merged_records.csv, <results>/qc_records.csv,
        <results>/qc_passed.fasta, the reference COX1 fasta
Writes: <out>/haplotypes.csv, <out>/record_checks.csv,
        <out>/barcode_gap.csv, <out>/clusters.csv, <out>/species_summary.txt

Usage: python scripts/species_analysis.py --results results/genus --out results/species
"""
import argparse
import csv
import math
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_qc import best_alignment, make_aligner  # noqa: E402

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
UNNAMED_WORDS = {"sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "nr", "nr."}


def species_name(label):
    """Two-word binomial if the label has a real species epithet, else ''."""
    words = (label or "").split()
    if len(words) >= 2 and re.fullmatch(r"[a-z][a-z-]+", words[1]) \
            and words[1] not in UNNAMED_WORDS:
        return f"{words[0]} {words[1]}"
    return ""


def project(aligner, ref, seq):
    """Sequence placed on reference coordinates; '-' where not covered."""
    best = best_alignment(aligner, ref, seq)
    out = ["-"] * len(ref)
    if best is None:
        return "".join(out)
    _, aln, oriented = best
    for (ts, te), (qs, _) in zip(*aln.aligned):
        out[ts:te] = oriented[qs:qs + (te - ts)]
    return "".join(out)


def encode(projected):
    return np.array([CODE.get(c, -1) for c in projected], dtype=np.int8)


def k2p_row(a, M, min_overlap):
    """K2P distance from encoded sequence a to every row of M (NaN if overlap too short)."""
    valid = (a >= 0) & (M >= 0)
    n = valid.sum(axis=1)
    diff = valid & (a != M)
    transition = diff & ((a % 2) == (M % 2))  # A<->G or C<->T
    ts = transition.sum(axis=1)
    tv = diff.sum(axis=1) - ts
    with np.errstate(divide="ignore", invalid="ignore"):
        P, Q = ts / n, tv / n
        w1, w2 = 1 - 2 * P - Q, 1 - 2 * Q
        d = -0.5 * np.log(w1) - 0.25 * np.log(w2)
    d = np.where((w1 > 0) & (w2 > 0), d, np.inf)  # saturated
    d[n < min_overlap] = np.nan
    return d


class UnionFind:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, i):
        while self.p[i] != i:
            self.p[i] = self.p[self.p[i]]
            i = self.p[i]
        return i

    def union(self, i, j):
        self.p[self.find(i)] = self.find(j)


def fmt(x):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    return f"{abs(max(float(x), 0.0)):.4f}"  # abs() avoids printing -0.0000


def main():
    p = argparse.ArgumentParser(description="Species-level COI analysis")
    p.add_argument("--results", default=os.path.join("results", "genus"))
    p.add_argument("--out", default=os.path.join("results", "species"))
    p.add_argument("--reference", default=os.path.join("data", "reference", "coi_reference.fasta"))
    p.add_argument("--min-overlap", type=int, default=300)
    p.add_argument("--threshold", type=float, default=0.03,
                   help="K2P distance for clustering (default 0.03)")
    p.add_argument("--synonyms", default="",
                   help="CSV of published synonymies (name, accepted_name); applied to names when given. "
                        "Without it the submitted names are used, as needed to describe labelling.")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ref = str(SeqIO.read(args.reference, "fasta").seq).upper()
    with open(os.path.join(args.results, "merged_records.csv"), newline="") as f:
        meta = {r["record_id"]: r for r in csv.DictReader(f)}
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(args.results, "qc_passed.fasta"), "fasta")}

    # Place sequences on reference coordinates and collapse identical ones.
    aligner = make_aligner()
    hap_index, hap_seqs, members = {}, [], defaultdict(list)
    record_hap = {}
    for rid, s in seqs.items():
        proj = project(aligner, ref, s)
        if proj not in hap_index:
            hap_index[proj] = len(hap_seqs)
            hap_seqs.append(proj)
        h = hap_index[proj]
        members[h].append(rid)
        record_hap[rid] = h
    H = len(hap_seqs)
    M = np.vstack([encode(s) for s in hap_seqs]) if H else np.zeros((0, len(ref)), np.int8)
    print(f"{len(seqs)} sequences, {H} unique haplotypes on reference coordinates")

    D = np.full((H, H), np.nan, dtype=np.float32)
    for i in range(H):
        D[i] = k2p_row(M[i], M, args.min_overlap)
    np.fill_diagonal(D, 0.0)

    synonyms = {}
    if args.synonyms:
        with open(args.synonyms, newline="") as f:
            synonyms = {r["name"]: r["accepted_name"] for r in csv.DictReader(f)}
    names = {}
    for rid in seqs:
        n = species_name(meta.get(rid, {}).get("species_label", ""))
        names[rid] = synonyms.get(n, n)
    hap_names = {h: Counter(names[r] for r in members[h] if names[r]) for h in range(H)}

    # Haplotypes shared by more than one species name.
    shared = {h for h in range(H) if len(hap_names[h]) > 1}

    # Nearest same-name and other-name distances per haplotype and name.
    def nearest(h, name, same):
        best_d, best_name = math.inf, ""
        for j in range(H):
            d = D[h, j]
            if np.isnan(d):
                continue
            if j == h:
                if same and hap_names[h][name] > 1:
                    best_d, best_name = 0.0, name
                elif not same and any(n != name for n in hap_names[h]):
                    other = next(n for n in hap_names[h] if n != name)
                    return 0.0, other
                continue
            names_j = hap_names[j]
            if same and name in names_j and d < best_d:
                best_d, best_name = float(d), name
            if not same:
                others = [n for n in names_j if n != name]
                if others and d < best_d:
                    best_d, best_name = float(d), others[0]
        return best_d, best_name

    # Clusters.
    uf = UnionFind(H)
    close = np.argwhere(np.triu(np.nan_to_num(D, nan=np.inf) <= args.threshold, k=1))
    for i, j in close:
        uf.union(int(i), int(j))
    cluster_of = {h: uf.find(h) for h in range(H)}
    cluster_ids = {c: f"C{k + 1:03d}" for k, c in enumerate(sorted(set(cluster_of.values())))}
    cluster_records = defaultdict(list)
    for rid, h in record_hap.items():
        cluster_records[cluster_ids[cluster_of[h]]].append(rid)

    # Record-level checks.
    rows, nn_cache = [], {}
    for rid in seqs:
        h, name = record_hap[rid], names[rid]
        row = {"record_id": rid,
               "species_label": meta.get(rid, {}).get("species_label", ""),
               "species_name": name, "haplotype": f"H{h + 1:04d}",
               "cluster": cluster_ids[cluster_of[h]],
               "nearest_same_species": "", "nearest_other_species": "",
               "nearest_other_species_name": "", "check": ""}
        flags = []
        if not name:
            flags.append("UNNAMED")
        else:
            key = (h, name)
            if key not in nn_cache:
                nn_cache[key] = (nearest(h, name, True), nearest(h, name, False))
            (ds, _), (do, on) = nn_cache[key]
            row["nearest_same_species"] = fmt(ds)
            row["nearest_other_species"] = fmt(do)
            row["nearest_other_species_name"] = on
            if h in shared:
                flags.append("SHARED_HAPLOTYPE")
            if do < ds:
                flags.append("NEAREST_IS_OTHER_SPECIES")
            if math.isinf(ds):
                flags.append("NO_CONSPECIFIC_FOR_COMPARISON")
        row["check"] = ";".join(flags) or "OK"
        rows.append(row)

    with open(os.path.join(args.out, "record_checks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["record_id"])
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(args.out, "haplotypes.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["haplotype", "cluster", "n_records", "species_names", "ref_start", "ref_end", "records"])
        for h in range(H):
            covered = [i for i, c in enumerate(hap_seqs[h]) if c != "-"]
            w.writerow([f"H{h + 1:04d}", cluster_ids[cluster_of[h]], len(members[h]),
                        "; ".join(f"{n} ({c})" for n, c in hap_names[h].most_common()),
                        covered[0] + 1 if covered else "", covered[-1] + 1 if covered else "",
                        " ".join(members[h])])

    # Barcode gap per named species.
    by_name = defaultdict(set)
    for rid, n in names.items():
        if n:
            by_name[n].add(record_hap[rid])
    gap_rows = []
    for n, haps in sorted(by_name.items()):
        haps = sorted(haps)
        n_rec = sum(1 for r in names if names[r] == n)
        intra = [D[i, j] for a, i in enumerate(haps) for j in haps[a + 1:] if not np.isnan(D[i, j])]
        max_intra = float(max(intra)) if intra else (0.0 if n_rec > 1 else math.nan)
        inter, nn_name = math.inf, ""
        for i in haps:
            for j in range(H):
                others = [x for x in hap_names[j] if x != n]
                if others and not np.isnan(D[i, j]) and D[i, j] < inter:
                    inter, nn_name = float(D[i, j]), others[0]
        gap = inter - max_intra if not (math.isnan(max_intra) or math.isinf(inter)) else math.nan
        gap_rows.append({"species": n, "records": n_rec, "haplotypes": len(haps),
                         "max_intraspecific": fmt(max_intra), "min_interspecific": fmt(inter),
                         "nearest_species": nn_name, "barcode_gap": fmt(gap),
                         "gap_present": "" if math.isnan(gap) else ("yes" if gap > 0 else "no")})
    with open(os.path.join(args.out, "barcode_gap.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(gap_rows[0].keys()) if gap_rows else ["species"])
        w.writeheader()
        w.writerows(gap_rows)

    # Cluster composition.
    cl_rows = []
    for cid, recs in sorted(cluster_records.items(), key=lambda kv: -len(kv[1])):
        comp = Counter(names[r] or "(unnamed)" for r in recs)
        named = [n for n in comp if n != "(unnamed)"]
        cl_rows.append({"cluster": cid, "records": len(recs),
                        "named_species": len(named),
                        "composition": "; ".join(f"{n} ({c})" for n, c in comp.most_common())})
    with open(os.path.join(args.out, "clusters.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cluster", "records", "named_species", "composition"])
        w.writeheader()
        w.writerows(cl_rows)

    # Summary.
    checks = Counter(fl for r in rows for fl in r["check"].split(";"))
    split = [n for n in by_name if len({cluster_of[h] for h in by_name[n]}) > 1]
    mixed = [c for c in cl_rows if c["named_species"] > 1]
    with_gap = [g for g in gap_rows if g["gap_present"]]
    lines = [f"Names: {'accepted (synonymies from ' + args.synonyms + ')' if synonyms else 'as submitted'}",
             f"Sequences: {len(seqs)}  Haplotypes: {H}  Named species: {len(by_name)}",
             f"Distance: K2P on shared sites (min overlap {args.min_overlap} bp); "
             f"clusters at {args.threshold:.3f}",
             f"Records without a species name: {checks.get('UNNAMED', 0)}",
             f"Records on haplotypes shared by more than one species name: {checks.get('SHARED_HAPLOTYPE', 0)}",
             f"Records whose nearest neighbour is another species: {checks.get('NEAREST_IS_OTHER_SPECIES', 0)}",
             f"Species with barcode gap (of {len(with_gap)} testable): "
             f"{sum(g['gap_present'] == 'yes' for g in with_gap)}",
             f"Species split over more than one cluster: {len(split)}"
             + (f" ({', '.join(sorted(split))})" if split else ""),
             f"Clusters containing more than one species name: {len(mixed)}"]
    for c in mixed[:15]:
        lines.append(f"  {c['cluster']} ({c['records']} records): {c['composition']}")
    lines.append("Species without a barcode gap:")
    for g in with_gap:
        if g["gap_present"] == "no":
            lines.append(f"  {g['species']}: max within {g['max_intraspecific']}, "
                         f"nearest {g['nearest_species']} at {g['min_interspecific']}")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "species_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
