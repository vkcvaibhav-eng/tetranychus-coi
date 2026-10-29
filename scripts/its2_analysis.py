"""Module 6: nuclear ITS2 check of the COI results.

ITS2 (nuclear ribosomal internal transcribed spacer 2) is inherited
independently of mitochondrial COI. Where the two genes disagree for the
same specimen, the COI result may reflect shared or introgressed
mitochondria rather than a wrong species name.

Steps:
1. ITS2 is cut from each GenBank record: an annotated ITS2 feature if
   present; otherwise the stretch between annotated 5.8S and 28S genes;
   otherwise the part that aligns to a reference ITS2 taken from an
   annotated T. urticae record.
2. Distances between unique ITS2 sequences from local alignments
   (1 - identity over aligned columns, gaps included), with a minimum
   aligned length; single-linkage clusters at --threshold.
3. Distances between species (after published synonymies).
4. Specimens with both genes are paired by first author and isolate/strain/
   voucher name, and their COI and ITS2 clusters are compared.

Reads:  data/its2/*.gb, data/ncbi_*.gb, results/species/record_checks.csv,
        results/species/clusters.csv, data/anchors/synonyms.csv
Writes: results/its2/its2_records.csv, its2_clusters.csv,
        species_distances.csv, coi_its2_pairs.csv, its2_summary.txt

Usage: python scripts/its2_analysis.py [--threshold 0.015]
"""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import glob
import os
import re
from collections import Counter, defaultdict
from statistics import median

from Bio import SeqIO
from Bio.Align import PairwiseAligner

ITS2_PAT = re.compile(r"^internal transcribed spacer 2(, ITS2)?$|^ITS2$|^ITS 2$", re.I)


def make_aligner():
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 2
    a.mismatch_score = -3
    a.open_gap_score = -5
    a.extend_gap_score = -2
    return a


def feature_label(f):
    q = f.qualifiers
    return (q.get("product") or q.get("note") or q.get("gene") or [""])[0].strip()


def specimen_key(rec):
    src = next((f for f in rec.features if f.type == "source"), None)
    q = src.qualifiers if src else {}
    sid = next((q[k][0] for k in ("isolate", "strain", "specimen_voucher") if k in q), "")
    refs = rec.annotations.get("references", [])
    author = refs[0].authors.split(",")[0].strip().lower() if refs and refs[0].authors else ""
    sid = re.sub(r"\W", "", sid.lower())
    return (author, sid) if author and sid else None


def country(rec):
    src = next((f for f in rec.features if f.type == "source"), None)
    q = src.qualifiers if src else {}
    return (q.get("geo_loc_name") or q.get("country") or [""])[0].split(":")[0]


def extract_its2(rec, aligner, ref):
    seq = str(rec.seq).upper()
    for f in rec.features:
        if f.type in ("misc_RNA", "misc_feature") and ITS2_PAT.match(feature_label(f)):
            return str(f.extract(rec.seq)).upper(), "ITS2 feature"
    s58 = [f for f in rec.features if f.type == "rRNA" and "5.8S" in feature_label(f)]
    s28 = [f for f in rec.features if f.type == "rRNA" and "28S" in feature_label(f)]
    if s58 and s28 and int(s58[0].location.end) < int(s28[0].location.start):
        return seq[int(s58[0].location.end):int(s28[0].location.start)], "between 5.8S and 28S"
    if ref:
        aln = next(iter(aligner.align(ref, seq)), None)
        if aln is not None and aln.score > 100:
            qb = aln.aligned[1]
            return seq[int(qb[0][0]):int(qb[-1][1])], "aligned to reference ITS2"
    return "", "not found"


def distance(aligner, a, b, min_len):
    aln = next(iter(aligner.align(a, b)), None)
    if aln is None:
        return None
    tb, qb = aln.aligned
    matches = columns = 0
    for i, ((ts, te), (qs, qe)) in enumerate(zip(tb, qb)):
        matches += sum(x == y for x, y in zip(a[ts:te], b[qs:qe]))
        columns += te - ts
        if i + 1 < len(tb):
            columns += (tb[i + 1][0] - te) + (qb[i + 1][0] - qe)
    if tb[-1][1] - tb[0][0] < min_len:
        return None
    return 1 - matches / columns


_HAPS, _MIN = [], 0


def _row(i):
    """Distances from haplotype i to all later haplotypes (run in worker processes)."""
    aligner = make_aligner()
    out = []
    for j in range(i + 1, len(_HAPS)):
        d = distance(aligner, _HAPS[i], _HAPS[j], _MIN)
        if d is not None:
            out.append((i, j, d))
    return out


def main():
    p = argparse.ArgumentParser(description="ITS2 check of the COI results")
    p.add_argument("--its2-dir", default=os.path.join("data", "its2"))
    p.add_argument("--coi-glob", default=os.path.join("data", "ncbi_*.gb"))
    p.add_argument("--species", default=os.path.join("results", "species"))
    p.add_argument("--anchors", default=os.path.join("data", "anchors"))
    p.add_argument("--out", default=os.path.join("results", "its2"))
    p.add_argument("--threshold", type=float, default=0.015)
    p.add_argument("--min-length", type=int, default=150)
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    syn_path = os.path.join(args.anchors, "synonyms.csv")
    synonyms = {}
    if os.path.exists(syn_path):
        with open(syn_path, newline="") as f:
            synonyms = {r["name"]: r["accepted_name"] for r in csv.DictReader(f)}

    def canon(label):
        words = (label or "").split()
        name = " ".join(words[:2]) if len(words) >= 2 and words[1].islower() \
            and words[1] not in ("sp.", "sp", "cf.", "aff.") else ""
        return synonyms.get(name, name)

    records = {}
    for path in sorted(glob.glob(os.path.join(args.its2_dir, "*.gb"))):
        for rec in SeqIO.parse(path, "genbank"):
            records.setdefault(rec.id, rec)

    # Reference ITS2: the most common length among annotated T. urticae ITS2.
    aligner = make_aligner()
    annotated = []
    for rec in records.values():
        if canon(rec.annotations.get("organism", "")) != "Tetranychus urticae":
            continue
        for f in rec.features:
            if f.type in ("misc_RNA", "misc_feature") and ITS2_PAT.match(feature_label(f)):
                annotated.append(str(f.extract(rec.seq)).upper())
    ref = Counter(annotated).most_common(1)[0][0] if annotated else ""
    print(f"{len(records)} ITS2 records; reference ITS2 {len(ref)} bp from {len(annotated)} annotated T. urticae")

    rows, haps, hap_of = [], {}, {}
    for rid, rec in records.items():
        its2, method = extract_its2(rec, aligner, ref)
        its2 = re.sub(r"[^ACGTRYKMSWBDHVN]", "", its2)
        row = {"record_id": rid, "species_label": rec.annotations.get("organism", ""),
               "accepted_name": canon(rec.annotations.get("organism", "")),
               "country": country(rec), "its2_length": len(its2), "extraction": method,
               "haplotype": "", "cluster": ""}
        if len(its2) >= args.min_length:
            if its2 not in haps:
                haps[its2] = len(haps)
            hap_of[rid] = haps[its2]
            row["haplotype"] = f"I{haps[its2] + 1:04d}"
        rows.append(row)
    hap_seqs = sorted(haps, key=haps.get)
    H = len(hap_seqs)
    print(f"{len(hap_of)} usable ITS2 sequences, {H} unique")

    global _HAPS, _MIN
    _HAPS, _MIN = hap_seqs, args.min_length
    D = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for chunk in pool.map(_row, range(H), chunksize=4):
            for i, j, d in chunk:
                D[(i, j)] = d

    parent = list(range(H))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for (i, j), d in D.items():
        if d <= args.threshold:
            parent[find(i)] = find(j)
    roots = sorted({find(i) for i in range(H)})
    cid = {r: f"I{k + 1:03d}" for k, r in enumerate(roots)}
    for row in rows:
        if row["record_id"] in hap_of:
            row["cluster"] = cid[find(hap_of[row["record_id"]])]

    with open(os.path.join(args.out, "its2_records.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    comp = defaultdict(Counter)
    for row in rows:
        if row["cluster"]:
            comp[row["cluster"]][row["accepted_name"] or "(unnamed)"] += 1
    cl_rows = sorted(({"cluster": c, "records": sum(v.values()),
                       "named_species": len([n for n in v if n != "(unnamed)"]),
                       "composition": "; ".join(f"{n} ({k})" for n, k in v.most_common())}
                      for c, v in comp.items()), key=lambda r: -r["records"])
    with open(os.path.join(args.out, "its2_clusters.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cluster", "records", "named_species", "composition"])
        w.writeheader()
        w.writerows(cl_rows)
    dominant = {c["cluster"]: c["composition"].split(" (")[0] for c in cl_rows}

    # Species-pair distances.
    names_of_hap = defaultdict(set)
    for row in rows:
        if row["record_id"] in hap_of and row["accepted_name"]:
            names_of_hap[hap_of[row["record_id"]]].add(row["accepted_name"])
    pair_d = defaultdict(list)
    for (i, j), d in D.items():
        for a in names_of_hap[i]:
            for b in names_of_hap[j]:
                pair_d[tuple(sorted((a, b)))].append(d)
    sp_rows = [{"species_a": a, "species_b": b, "comparisons": len(v),
                "min": f"{min(v):.4f}", "median": f"{median(v):.4f}", "max": f"{max(v):.4f}"}
               for (a, b), v in sorted(pair_d.items())]
    with open(os.path.join(args.out, "species_distances.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["species_a", "species_b", "comparisons", "min", "median", "max"])
        w.writeheader()
        w.writerows(sp_rows)

    # Pair specimens across genes.
    with open(os.path.join(args.species, "record_checks.csv"), newline="") as f:
        coi_checks = {r["record_id"]: r for r in csv.DictReader(f)}
    with open(os.path.join(args.species, "clusters.csv"), newline="") as f:
        coi_dom = {r["cluster"]: r["composition"].split(" (")[0] for r in csv.DictReader(f)}
    coi_by_key = {}
    for path in sorted(glob.glob(args.coi_glob)):
        for rec in SeqIO.parse(path, "genbank"):
            k = specimen_key(rec)
            if k and rec.id in coi_checks:
                coi_by_key.setdefault(k, rec.id)
    its_rows = {r["record_id"]: r for r in rows}
    pairs = []
    for rid, rec in records.items():
        k = specimen_key(rec)
        if not k or k not in coi_by_key or not its_rows[rid]["cluster"]:
            continue
        crid = coi_by_key[k]
        c = coi_checks[crid]
        coi_name = canon(c["species_label"])
        its_name = its_rows[rid]["accepted_name"]
        coi_cluster_name = canon(coi_dom.get(c["cluster"], ""))
        its_cluster_name = canon(dominant.get(its_rows[rid]["cluster"], ""))
        pairs.append({"specimen": f"{k[0]}:{k[1]}", "coi_record": crid, "its2_record": rid,
                      "coi_label": coi_name, "its2_label": its_name,
                      "coi_cluster": c["cluster"], "coi_cluster_dominant": coi_cluster_name,
                      "its2_cluster": its_rows[rid]["cluster"], "its2_cluster_dominant": its_cluster_name,
                      "country": its_rows[rid]["country"],
                      "genes_agree": "yes" if coi_cluster_name == its_cluster_name else "no"})
    with open(os.path.join(args.out, "coi_its2_pairs.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pairs[0].keys()) if pairs else ["specimen"])
        w.writeheader()
        w.writerows(pairs)

    # Summary.
    lines = [f"ITS2 records: {len(records)}; usable: {len(hap_of)}; unique: {H}",
             "Extraction: " + ", ".join(f"{k} {v}" for k, v in Counter(r["extraction"] for r in rows).most_common()),
             f"Distance: local alignment, min aligned length {args.min_length} bp; clusters at {args.threshold}",
             f"ITS2 clusters: {len(cl_rows)}; with more than one species name: "
             f"{sum(c['named_species'] > 1 for c in cl_rows)}"]
    for c in cl_rows[:12]:
        lines.append(f"  {c['cluster']} ({c['records']}): {c['composition']}")
    lines.append("Key species pairs (ITS2 distance min / median):")
    for a, b in [("Tetranychus turkestani", "Tetranychus urticae"),
                 ("Tetranychus truncatus", "Tetranychus urticae"),
                 ("Tetranychus kanzawai", "Tetranychus parakanzawai"),
                 ("Tetranychus macfarlanei", "Tetranychus malaysiensis"),
                 ("Tetranychus ludeni", "Tetranychus neocaledonicus")]:
        v = pair_d.get(tuple(sorted((a, b))))
        between = f"{min(v):.4f} / {median(v):.4f}" if v else "no comparison"
        within = []
        for n in (a, b):
            w_ = pair_d.get((n, n))
            within.append(f"within {n.split()[1]} {median(w_):.4f}" if w_ else f"within {n.split()[1]} n/a")
        lines.append(f"  {a.split()[1]} vs {b.split()[1]}: {between}  ({'; '.join(within)})")
    agree = Counter(p["genes_agree"] for p in pairs)
    lines.append(f"Specimens with both genes: {len(pairs)}; COI and ITS2 clusters agree: "
                 f"{agree.get('yes', 0)}, disagree: {agree.get('no', 0)}")
    dis = Counter((p["coi_label"], p["coi_cluster_dominant"], p["its2_cluster_dominant"], p["country"])
                  for p in pairs if p["genes_agree"] == "no")
    for (lab, cdom, idom, ctry), n in dis.most_common(12):
        lines.append(f"  {n} x labelled {lab or '(unnamed)'} ({ctry}): COI groups with {cdom}, ITS2 with {idom}")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "its2_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
