"""Supplementary tables S1-S6 as one Excel workbook.

S1  Anchor studies and anchor accessions
S2  Metadata and quality control for every record (accession versions)
S3  Haplotype, cluster, nearest neighbours and anchor verdict for every
    record that passed quality control
S4  Sensitivity: K2P thresholds, anchor evidence and taxonomic treatment
S5  Identification test: every query, outcomes by species, raw -> cleaned
S6  Records disagreeing with the anchored species of their cluster, and
    COI-ITS2 comparison for specimens sequenced for both genes

Usage: python scripts/supplementary_tables.py [--out results/supplementary]
"""
import argparse
import csv
import os
import re

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter


def read(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def clean(text):
    """Drop notes on where a detail was found (they belong in the data file, not the table)."""
    text = re.sub(r"\s*\((?:per|claimed by)[^)]*\)", "", text or "")
    text = re.sub(r"\s*claimed by [A-Za-z ,]+ only", "", text).strip()
    text = re.sub(r"^uncertain", "not stated", text)
    return text[:1].upper() + text[1:]


def add_sheet(wb, title, caption, blocks):
    ws = wb.create_sheet(title)
    ws.append([caption])
    ws["A1"].font = Font(bold=True)
    for heading, header, rows in blocks:
        ws.append([])
        if heading:
            ws.append([heading])
            ws.cell(ws.max_row, 1).font = Font(italic=True)
        ws.append(header)
        for c in range(1, len(header) + 1):
            ws.cell(ws.max_row, c).font = Font(bold=True)
        for r in rows:
            ws.append([r.get(h, "") if isinstance(r, dict) else r[i] for i, h in enumerate(header)])
    for col in range(1, ws.max_column + 1):
        width = max((len(str(ws.cell(r, col).value or "")) for r in range(2, ws.max_row + 1)), default=8)
        ws.column_dimensions[get_column_letter(col)].width = min(max(width + 2, 8), 60)
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(vertical="top")
    ws.freeze_panes = None
    return ws


def main():
    p = argparse.ArgumentParser(description="Supplementary tables")
    p.add_argument("--out", default=os.path.join("results", "supplementary"))
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    studies = read("data/anchors/anchor_studies.csv")
    anchors = read("data/anchors/anchor_records.csv")
    status = {r["accession"]: r for r in read("results/anchors/anchor_status.csv")}
    merged = {r["record_id"]: r for r in read("results/genus/merged_records.csv")}
    qc = read("results/genus/qc_records.csv")
    checks = {r["record_id"]: r for r in read("results/species/record_checks.csv")}
    verdicts = {r["record_id"]: r for r in read("results/anchors/record_anchor_check.csv")}
    nn = {r["record_id"]: r for r in read("results/robustness/coi_nearest_neighbour.csv")}

    wb = Workbook()
    wb.remove(wb.active)

    # S1
    def detail(s, key):
        return clean(s[key]) if s["evidence_level"] in "AB" else "Not verified"

    s1_studies = [{"Study": s["citation"],
                   "DOI": s["doi"] if s["doi"].startswith("10.") and "check" not in s["doi"] else "",
                   "Evidence level": s["evidence_level"],
                   "Used for decisions": "yes" if s["evidence_level"] in "AB" else "no (lead only)",
                   "Morphological identification": detail(s, "morphology_verified"),
                   "Males examined": detail(s, "males_examined"), "Identifier": detail(s, "identifier_named"),
                   "Voucher repository": detail(s, "voucher_linked"), "Crossing tests": detail(s, "crossing_tests"),
                   "Material sequenced": "", "Supporting passage": "",
                   "Anchor accessions": sum(a["study_id"] == s["study_id"] for a in anchors)}
                  for s in studies]
    by_id = {s["study_id"]: s for s in studies}
    s1_acc = []
    for a in anchors:
        st = status.get(a["accession"], {})
        rid = st.get("record_id", "")
        s1_acc.append({"Accession": rid or a["accession"], "Study": by_id[a["study_id"]]["citation"].split(" (")[0]
                       + " (" + by_id[a["study_id"]]["citation"].split(" (")[1].split(")")[0] + ")",
                       "Evidence level": by_id[a["study_id"]]["evidence_level"],
                       "Name in record": merged.get(rid, {}).get("species_label", ""),
                       "Accepted name": st.get("species_name", ""),
                       "Passed QC": st.get("passed_qc", ""), "Cluster": st.get("cluster", ""),
                       "Anchored species of cluster": st.get("cluster_anchored_species", ""),
                       "Consistency": st.get("anchor_consistency", "")})
    h1 = list(s1_studies[0].keys())
    h1b = list(s1_acc[0].keys())
    add_sheet(wb, "S1 Anchors", "Table S1. Anchor studies, evidence levels and COI accessions.",
              [("Studies", h1, s1_studies), ("Anchor accessions", h1b, s1_acc)])

    # S2
    s2 = []
    for r in qc:
        m = merged.get(r["record_id"], {})
        s2.append({"Record": r["record_id"], "GenBank accession": m.get("genbank_accession", ""),
                   "BOLD process ID": m.get("bold_processid", ""), "Source": r["source"],
                   "Name in record": r["species_label"], "Country": m.get("country", ""),
                   "Metadata score (0-7)": m.get("metadata_score", ""),
                   "Missing items": m.get("metadata_missing", "").replace(";", "; "),
                   "Automated BOLD identification": m.get("identification_automated", ""),
                   "Coordinate note": m.get("coordinate_source", ""),
                   "Length (bp)": r["seq_length"], "COX1 region": r["cox1_half"],
                   "Reference start": r["ref_start"], "Reference end": r["ref_end"],
                   "Identity to reference": r["identity"], "Reading frame intact": r["reading_frame_intact"],
                   "QC flags": r["qc_flags"].replace(";", "; "), "QC result": r["qc_status"]})
    add_sheet(wb, "S2 Records and QC", "Table S2. Metadata and quality-control result for every record. "
              "Records are identified by accession and version, or by BOLD process ID for BOLD-only records.",
              [("", list(s2[0].keys()), s2)])

    # S3
    s3 = []
    for rid, c in checks.items():
        v, n = verdicts.get(rid, {}), nn.get(rid, {})
        s3.append({"Record": rid, "Name in record": c["species_label"], "Accepted name": v.get("accepted_name", ""),
                   "Haplotype": c["haplotype"], "Cluster (K2P 0.03)": c["cluster"],
                   "Nearest same name (K2P)": c["nearest_same_species"],
                   "Nearest other name (K2P)": c["nearest_other_species"],
                   "Nearest other name": c["nearest_other_species_name"], "Checks": c["check"].replace(";", "; "),
                   "Anchored species of cluster": v.get("anchored_species", ""), "Anchor": v.get("is_anchor", ""),
                   "Verdict": v.get("verdict", ""), "Submission batch": n.get("batch", ""),
                   "Nearest record name": n.get("nn_name_any", ""),
                   "Nearest record name, other batches": n.get("nn_name_other_batch", "")})
    add_sheet(wb, "S3 Clusters", "Table S3. Haplotype, cluster, nearest neighbours and anchor verdict for every "
              "record that passed quality control.", [("", list(s3[0].keys()), s3)])

    # S4
    thr = read("results/robustness/thresholds.csv")
    var = read("results/robustness/variants.csv")
    add_sheet(wb, "S4 Sensitivity", "Table S4. Sensitivity of clusters, anchor verdicts and identification "
              "accuracy to the K2P threshold, the evidence level of anchors and the taxonomic treatment.",
              [("K2P threshold (clusters and largest accepted match distance)", list(thr[0].keys()), thr),
               ("Anchor evidence and taxonomic treatment, K2P 0.03", list(var[0].keys()), var)])

    # S5
    qs = read("results/robustness/queries_main.csv")
    sp = read("results/robustness/identification_by_species.csv")
    tr = read("results/robustness/transitions.csv")
    ci = read("results/robustness/identification_intervals.csv")
    add_sheet(wb, "S5 Identification", "Table S5. Identification of anchor specimens by closest COI match "
              "(K2P 0.03, leaving out each anchor study in turn).",
              [("Outcomes, raw database (rows) against cleaned database (columns)", list(tr[0].keys()), tr),
               ("Percentages with 95% intervals from a bootstrap over anchor studies", list(ci[0].keys()), ci),
               ("By species", list(sp[0].keys()), sp),
               ("Every query", list(qs[0].keys()), qs)])

    # S6
    dis = read("results/robustness/disagreements.csv")
    pairs = read("results/its2/coi_its2_pairs.csv")
    plab = {r["specimen"]: r for r in read("results/robustness/paired_label_agreement.csv")}
    for r in pairs:
        x = plab.get(r["specimen"], {})
        r["coi_nn_other_batch"] = x.get("coi_nn_other_batch", "")
        r["its2_nn_other_batch"] = x.get("its2_nn_other_batch", "")
    blocks = [("Records disagreeing with the anchored species of their cluster", list(dis[0].keys()), dis)]
    if pairs:
        blocks.append(("Specimens sequenced for COI and ITS2", list(pairs[0].keys()), pairs))
    add_sheet(wb, "S6 Disagreements and ITS2", "Table S6. Records that disagree with the anchored species of "
              "their cluster, and COI-ITS2 comparison for specimens sequenced for both genes.", blocks)

    path = os.path.join(args.out, "Supplementary_Tables_S1-S6.xlsx")
    wb.save(path)
    print(f"Written {path}")


if __name__ == "__main__":
    main()
