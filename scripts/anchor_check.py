"""Module 5: check distance clusters against anchor records.

Anchor records are sequences from published studies whose specimens were
identified morphologically (data/anchors/). Names are compared after applying published synonymies listed in
data/anchors/synonyms.csv. Each study has an evidence level:
  A  morphology with males examined and vouchers in a named collection
  B  morphology stated, vouchers not confirmed
  C  lead only; listed but not used for decisions

For every cluster from species_analysis.py the module lists the species
names carried by anchors of the best evidence level present in that cluster
(A outranks B; a B anchor naming a different species is a conflict), and then
classifies each record:
  AGREES                 its name matches the only anchored species in the cluster
  AGREES_SHARED_CLUSTER  its name is one of several anchored species in the cluster
                         (COI does not separate them)
  DISAGREES              its name differs from every anchored species in the cluster
  UNNAMED_SUGGESTED      no species name; the cluster has one anchored species
  NO_ANCHOR              the cluster has no anchor

Reads:  data/anchors/anchor_studies.csv, data/anchors/anchor_records.csv (accession, study_id and,
        optionally, species: the name given in the study when the record carries another label),
        <genus>/merged_records.csv, <species>/record_checks.csv, <species>/clusters.csv
Writes: <out>/anchor_status.csv, <out>/cluster_anchors.csv,
        <out>/record_anchor_check.csv, <out>/anchor_summary.txt

Usage: python scripts/anchor_check.py [--levels AB]
"""
import argparse
import csv
import os
from collections import Counter, defaultdict


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser(description="Check clusters against anchor records")
    p.add_argument("--anchors", default=os.path.join("data", "anchors"))
    p.add_argument("--genus", default=os.path.join("results", "genus"))
    p.add_argument("--species", default=os.path.join("results", "species"))
    p.add_argument("--out", default=os.path.join("results", "anchors"))
    p.add_argument("--levels", default="AB", help="Evidence levels used for decisions")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    studies = {s["study_id"]: s for s in read(os.path.join(args.anchors, "anchor_studies.csv"))}
    anchor_rows = read(os.path.join(args.anchors, "anchor_records.csv"))
    merged = read(os.path.join(args.genus, "merged_records.csv"))
    checks = {r["record_id"]: r for r in read(os.path.join(args.species, "record_checks.csv"))}
    clusters = {c["cluster"]: c for c in read(os.path.join(args.species, "clusters.csv"))}
    by_acc = {r["genbank_accession"]: r for r in merged if r["genbank_accession"]}
    syn_path = os.path.join(args.anchors, "synonyms.csv")
    synonyms = {r["name"]: r["accepted_name"] for r in read(syn_path)} if os.path.exists(syn_path) else {}

    def canon(name):
        return synonyms.get(name, name)

    # Anchor status and anchored species per cluster.
    status, cluster_anchor = [], defaultdict(Counter)
    anchor_record_ids = {}
    for a in anchor_rows:
        st = studies[a["study_id"]]
        rec = by_acc.get(a["accession"])
        row = {"accession": a["accession"], "study_id": a["study_id"],
               "evidence_level": st["evidence_level"], "found": "no", "record_id": "",
               "species_name": "", "passed_qc": "", "cluster": ""}
        if rec:
            rid = rec["record_id"]
            chk = checks.get(rid)
            row.update(found="yes", record_id=rid,
                       species_name=canon((a.get("species") or "").strip()
                                          or " ".join(rec["species_label"].split()[:2])),
                       passed_qc="yes" if chk else "no", cluster=chk["cluster"] if chk else "")
            if chk:
                anchor_record_ids[rid] = row
                own = row["species_name"] if (a.get("species") or "").strip() else canon(chk["species_name"])
                if st["evidence_level"] in args.levels and own:
                    cluster_anchor[chk["cluster"]][(own, st["evidence_level"])] += 1
        status.append(row)

    def anchored_names(cluster):
        """Names at the best evidence level present in the cluster."""
        entries = cluster_anchor.get(cluster, {})
        if not entries:
            return []
        best = min(lvl for _, lvl in entries)  # "A" < "B"
        return sorted({name for name, lvl in entries if lvl == best})

    for row in status:
        names = anchored_names(row["cluster"]) if row["cluster"] else []
        row["cluster_anchored_species"] = "; ".join(names)
        if not row["cluster"]:
            row["anchor_consistency"] = "not in clusters"
        elif row["evidence_level"] not in args.levels:
            row["anchor_consistency"] = "lead only"
        elif names == [row["species_name"]]:
            row["anchor_consistency"] = "consistent"
        elif row["species_name"] in names:
            row["anchor_consistency"] = "cluster shared with other anchored species"
        else:
            row["anchor_consistency"] = "conflict"

    with open(os.path.join(args.out, "anchor_status.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(status[0].keys()))
        w.writeheader()
        w.writerows(status)

    # Cluster table.
    cl_rows = []
    for cid, c in clusters.items():
        names = anchored_names(cid)
        detail = "; ".join(f"{n} [{lvl}] x{k}" for (n, lvl), k in sorted(cluster_anchor.get(cid, {}).items()))
        cl_rows.append({"cluster": cid, "records": c["records"], "composition": c["composition"],
                        "anchored_species": "; ".join(names), "anchor_detail": detail,
                        "status": "no anchor" if not names else
                        ("one anchored species" if len(names) == 1 else "several anchored species")})
    with open(os.path.join(args.out, "cluster_anchors.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cl_rows[0].keys()))
        w.writeheader()
        w.writerows(cl_rows)

    # Record-level verdicts.
    rec_rows = []
    for rid, chk in checks.items():
        names = anchored_names(chk["cluster"])
        own = canon(chk["species_name"]) if chk["species_name"] else ""
        if not names:
            verdict = "NO_ANCHOR"
        elif not own:
            verdict = "UNNAMED_SUGGESTED" if len(names) == 1 else "UNNAMED_SEVERAL_POSSIBLE"
        elif names == [own]:
            verdict = "AGREES"
        elif own in names:
            verdict = "AGREES_SHARED_CLUSTER"
        else:
            verdict = "DISAGREES"
        rec_rows.append({"record_id": rid, "species_label": chk["species_label"], "accepted_name": own,
                         "cluster": chk["cluster"], "anchored_species": "; ".join(names),
                         "is_anchor": "yes" if rid in anchor_record_ids else "no", "verdict": verdict})
    with open(os.path.join(args.out, "record_anchor_check.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rec_rows[0].keys()))
        w.writeheader()
        w.writerows(rec_rows)

    # Summary.
    used = [s for s in status if s["evidence_level"] in args.levels]
    found = [s for s in used if s["found"] == "yes"]
    in_cl = [s for s in used if s["cluster"]]
    cons = Counter(s["anchor_consistency"] for s in in_cl)
    verdicts = Counter(r["verdict"] for r in rec_rows if r["is_anchor"] == "no")
    n_other = sum(verdicts.values())
    anchored_clusters = [c for c in cl_rows if c["anchored_species"]]
    lines = [f"Anchor studies: {len(studies)} (levels used for decisions: {args.levels})",
             f"Synonymies applied: {', '.join(f'{k} = {v}' for k, v in synonyms.items()) or 'none'}",
             f"Anchor accessions (levels {args.levels}): {len(used)}; found in data: {len(found)}; "
             f"passed QC and clustered: {len(in_cl)}",
             "Anchor consistency:"] + [f"  {k}: {v}" for k, v in cons.most_common()]
    conflicts = [s for s in in_cl if s["anchor_consistency"] == "conflict"]
    for s in conflicts:
        lines.append(f"    conflict: {s['accession']} ({s['study_id']}, {s['species_name']}) "
                     f"in cluster with {s['cluster_anchored_species']}")
    lines += [f"Clusters with anchors: {len(anchored_clusters)} of {len(cl_rows)} "
              f"(several anchored species: {sum(c['status'] == 'several anchored species' for c in cl_rows)})",
              f"Non-anchor records: {n_other}"]
    lines += [f"  {k}: {v} ({100 * v / n_other:.1f}%)" for k, v in verdicts.most_common()] if n_other else []
    lines.append("Clusters holding several anchored species:")
    for c in cl_rows:
        if c["status"] == "several anchored species":
            lines.append(f"  {c['cluster']} ({c['records']} records): {c['anchor_detail']}")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "anchor_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
