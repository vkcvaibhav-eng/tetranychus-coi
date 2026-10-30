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
2. Distances between unique ITS2 sequences from semi-global alignments
   (free end gaps; 1 - identity over the overlapping columns, internal gaps
   included), with a minimum overlap. Each record is assigned the name
   of its nearest other record (nearest-neighbour). Single-linkage clusters
   at --threshold are reported too, but they chain closely related species
   together and are not used for the COI comparison.
3. Distances between species (after published synonymies).
4. Specimens with both genes are paired by first author and isolate/strain/
   voucher name; the species of the specimen's COI cluster is compared with
   its ITS2 nearest-neighbour name.

Reads:  data/its2/*.gb, data/ncbi_*.gb, results/species/record_checks.csv,
        results/species/clusters.csv, data/anchors/synonyms.csv
Writes: results/its2/its2_records.csv, its2_clusters.csv,
        species_distances.csv, coi_its2_pairs.csv, its2_summary.txt

Usage: python scripts/its2_analysis.py [--threshold 0.015]
"""
import argparse
import csv
import hashlib
from concurrent.futures import ProcessPoolExecutor
import glob
import os
import re
import sys
from collections import Counter, defaultdict
from statistics import median

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from study_keys import submission_key  # noqa: E402

ITS2_PAT = re.compile(r"^internal transcribed spacer 2(, ITS2)?$|^ITS2$|^ITS 2$", re.I)


def make_aligner(local=False):
    """Semi-global by default: free end gaps, so the whole overlap is compared."""
    a = PairwiseAligner()
    a.mode = "local" if local else "global"
    if not local:
        a.end_gap_score = 0
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
    # aligner here is local: it finds the ITS2 stretch inside longer records
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
    p.add_argument("--distance-cache", default="",
                   help="Optional .npz file to store or reuse the ITS2 distances (checked against the haplotypes)")
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
    aligner = make_aligner(local=True)
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
    digest = hashlib.sha256(("\n".join(hap_seqs) + f"|{args.min_length}").encode()).hexdigest()
    if args.distance_cache and os.path.exists(args.distance_cache):
        z = np.load(args.distance_cache)
        if str(z["digest"]) == digest:
            D = {(int(i), int(j)): float(d) for i, j, d in zip(z["i"], z["j"], z["d"])}
            print(f"ITS2 distances read from {args.distance_cache}")
    if not D:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for chunk in pool.map(_row, range(H), chunksize=4):
                for i, j, d in chunk:
                    D[(i, j)] = d
        if args.distance_cache:
            keys = list(D)
            np.savez_compressed(args.distance_cache, digest=digest, i=np.array([k[0] for k in keys]),
                                j=np.array([k[1] for k in keys]), d=np.array([D[k] for k in keys]))

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

    # Nearest-neighbour name for every record (excluding the record itself).
    recs_of_hap = defaultdict(list)
    for row in rows:
        if row["record_id"] in hap_of:
            recs_of_hap[hap_of[row["record_id"]]].append(row)
    neigh = defaultdict(list)
    for (i, j), d in D.items():
        neigh[i].append((d, j))
        neigh[j].append((d, i))
    for row in rows:
        row["nn_distance"] = row["nn_name"] = row["nn_match"] = ""
        if row["record_id"] not in hap_of:
            continue
        h = hap_of[row["record_id"]]
        cands = [(0.0, r) for r in recs_of_hap[h] if r is not row]
        if not cands and neigh[h]:
            dmin = min(d for d, _ in neigh[h])
            cands = [(d, r) for d, j in neigh[h] if d == dmin for r in recs_of_hap[j]]
        named = Counter(r["accepted_name"] for _, r in cands if r["accepted_name"])
        if named:
            row["nn_distance"] = f"{cands[0][0]:.4f}"
            row["nn_name"] = named.most_common(1)[0][0]
            if row["accepted_name"]:
                row["nn_match"] = "yes" if row["nn_name"] == row["accepted_name"] else "no"

    # The same, but only records from other submission batches count as neighbours,
    # so that identical sequences deposited together cannot confirm each other.
    for row in rows:
        row["batch"] = submission_key(records[row["record_id"]]) or f"single:{row['record_id']}"
    sorted_neigh = {h: sorted(v, key=lambda x: x[0]) for h, v in neigh.items()}
    for row in rows:
        row["nn_name_other_batch"] = row["nn_distance_other_batch"] = row["nn_match_other_batch"] = ""
        if row["record_id"] not in hap_of:
            continue
        h, b = hap_of[row["record_id"]], row["batch"]
        cands = [(0.0, r) for r in recs_of_hap[h] if r is not row and r["batch"] != b]
        if not cands:
            dmin = None
            for d, j in sorted_neigh.get(h, []):
                if dmin is not None and d > dmin:
                    break
                found = [(d, r) for r in recs_of_hap[j] if r["batch"] != b]
                if found:
                    dmin = d
                    cands += found
        named = Counter(r["accepted_name"] for _, r in cands if r["accepted_name"])
        if named:
            row["nn_distance_other_batch"] = f"{cands[0][0]:.4f}"
            row["nn_name_other_batch"] = named.most_common(1)[0][0]
            if row["accepted_name"]:
                row["nn_match_other_batch"] = "yes" if row["nn_name_other_batch"] == row["accepted_name"] else "no"

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
    # Specimens are linked by first author and isolate, strain or voucher name. A link is
    # accepted only when exactly one COI record and one ITS2 record share the key and
    # their names (after synonymies) and countries do not conflict; other links are kept
    # in the table but flagged and left out of the comparison.
    coi_by_key, coi_country = defaultdict(list), {}
    for path in sorted(glob.glob(args.coi_glob)):
        for rec in SeqIO.parse(path, "genbank"):
            k = specimen_key(rec)
            if k and rec.id in coi_checks and rec.id not in coi_by_key[k]:
                coi_by_key[k].append(rec.id)
                coi_country[rec.id] = country(rec)
    its_rows = {r["record_id"]: r for r in rows}
    its_by_key = defaultdict(list)
    for rid, rec in records.items():
        k = specimen_key(rec)
        if k and its_rows[rid]["cluster"]:
            its_by_key[k].append(rid)
    pairs = []
    for k in sorted(set(coi_by_key) & set(its_by_key)):
        cois, itss = coi_by_key[k], its_by_key[k]
        link = "one-to-one" if len(cois) == 1 and len(itss) == 1 else \
            f"ambiguous ({len(cois)} COI, {len(itss)} ITS2 records)"
        for crid in cois:
            for rid in itss:
                c = coi_checks[crid]
                coi_name = canon(c["species_label"])
                its_name = its_rows[rid]["accepted_name"]
                coi_cluster_name = canon(coi_dom.get(c["cluster"], ""))
                its_cluster_name = its_rows[rid]["nn_name"]
                conflicts = []
                if coi_name and its_name and coi_name != its_name:
                    conflicts.append("names differ")
                cc, ic = coi_country.get(crid, ""), its_rows[rid]["country"]
                if cc and ic and cc != ic:
                    conflicts.append("countries differ")
                pairs.append({"specimen": f"{k[0]}:{k[1]}", "coi_record": crid, "its2_record": rid,
                              "coi_label": coi_name, "its2_label": its_name,
                              "coi_cluster": c["cluster"], "coi_cluster_dominant": coi_cluster_name,
                              "its2_cluster": its_rows[rid]["cluster"], "its2_nearest_species": its_cluster_name,
                              "country": ic or cc, "link": link, "link_conflict": "; ".join(conflicts),
                              "confident": "yes" if link == "one-to-one" and not conflicts else "no",
                              "genes_agree": "yes" if coi_cluster_name == its_cluster_name else "no"})
    with open(os.path.join(args.out, "coi_its2_pairs.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pairs[0].keys()) if pairs else ["specimen"])
        w.writeheader()
        w.writerows(pairs)

    # Summary.
    lines = [f"ITS2 records: {len(records)}; usable: {len(hap_of)}; unique: {H}",
             "Single-linkage clusters chain related species; species assignment uses nearest neighbours.",
             "Extraction: " + ", ".join(f"{k} {v}" for k, v in Counter(r["extraction"] for r in rows).most_common()),
             f"Distance: semi-global alignment (free end gaps), min overlap {args.min_length} bp; clusters at {args.threshold}",
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
    nn = [r for r in rows if r["nn_match"]]
    lines.append(f"ITS2 nearest-neighbour name matches own name: "
                 f"{sum(r['nn_match'] == 'yes' for r in nn)} of {len(nn)} named records "
                 f"({100 * sum(r['nn_match'] == 'yes' for r in nn) / max(len(nn), 1):.1f}%)")
    per_sp = defaultdict(Counter)
    for r in nn:
        per_sp[r["accepted_name"]][r["nn_match"]] += 1
    for sp, c in sorted(per_sp.items(), key=lambda kv: -sum(kv[1].values()))[:12]:
        lines.append(f"  {sp}: {c['yes']} of {c['yes'] + c['no']} match")
    nb = [r for r in rows if r["nn_match_other_batch"]]
    lines.append(f"Nearest neighbour from another submission batch matches own name: "
                 f"{sum(r['nn_match_other_batch'] == 'yes' for r in nb)} of {len(nb)} named records "
                 f"({100 * sum(r['nn_match_other_batch'] == 'yes' for r in nb) / max(len(nb), 1):.1f}%)")
    per_sp = defaultdict(Counter)
    for r in nb:
        per_sp[r["accepted_name"]][r["nn_match_other_batch"]] += 1
    for sp, c in sorted(per_sp.items(), key=lambda kv: -sum(kv[1].values()))[:12]:
        lines.append(f"  {sp}: {c['yes']} of {c['yes'] + c['no']} match")
    conf = [p for p in pairs if p["confident"] == "yes"]
    lines.append(f"Specimen links between genes: {len(pairs)} pairs for {len({p['specimen'] for p in pairs})} keys; "
                 f"one-to-one without conflict: {len(conf)}; ambiguous: "
                 f"{sum(p['link'] != 'one-to-one' for p in pairs)} pairs; one-to-one with conflicting names or "
                 f"countries: {sum(p['link'] == 'one-to-one' and p['confident'] == 'no' for p in pairs)}")
    agree = Counter(p["genes_agree"] for p in conf)
    lines.append("COI vs ITS2 for one-to-one specimens (COI cluster species vs ITS2 nearest species):")
    lines.append(f"  agree {agree.get('yes', 0)}, disagree {agree.get('no', 0)} of {len(conf)}")
    dis = Counter((p["coi_label"], p["coi_cluster_dominant"], p["its2_nearest_species"], p["country"])
                  for p in conf if p["genes_agree"] == "no")
    for (lab, cdom, idom, ctry), n in dis.most_common(12):
        lines.append(f"  {n} x labelled {lab or '(unnamed)'} ({ctry}): COI groups with {cdom or '-'}, "
                     f"ITS2 nearest {idom or '-'}")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "its2_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
