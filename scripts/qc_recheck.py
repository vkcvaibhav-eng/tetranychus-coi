"""Module 13: recheck of failed and frame-flagged sequences against a panel of references.

The quality screen (sequence_qc.py) aligns every sequence to one reference
COX1, so a frameshift or stop-codon flag can be an artefact of aligning a
divergent species. Here every sequence that failed the screen, and every
passing sequence with alignment gaps that are not a multiple of three, is
screened again against a panel of references: the annotated COX1 gene of one
mitochondrial genome per species in the downloaded GenBank files, plus the
main reference. The reference with the highest alignment score is kept and
the screen is repeated with it. Independently of any alignment, each sequence
is translated in all six reading frames (invertebrate mitochondrial code),
and the frame with the fewest internal stop codons is recorded.

Reads:  <results>/merged_records.csv, <results>/qc_records.csv, <data>/ncbi_*.gb, reference
Writes: <results>/qc_recheck.csv, <results>/qc_recheck_summary.txt

Usage: python scripts/qc_recheck.py [--results results/genus --data data --reference FILE]
"""
import argparse
import csv
import glob
import os
import sys
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_qc import GENETIC_CODE, check, make_aligner  # noqa: E402


def panel(data_dir, reference):
    """One COX1 gene per species from mitochondrial genomes, plus the main reference."""
    refs = {"main reference": str(SeqIO.read(reference, "fasta").seq).upper()}
    seen = set()
    for path in sorted(glob.glob(os.path.join(data_dir, "ncbi_*.gb"))):
        for rec in SeqIO.parse(path, "genbank"):
            if len(rec) < 2500:
                continue
            sp = " ".join(rec.annotations.get("organism", "").split()[:2])
            if sp in seen:
                continue
            for f in rec.features:
                if f.type in ("CDS", "gene") and [g.lower() for g in f.qualifiers.get("gene", [])] in (["cox1"], ["coi"],
                                                                                                    ["cox-1"]):
                    s = str(f.extract(rec.seq)).upper()
                    if len(s) > 1400:
                        refs[f"{sp} ({rec.id})"] = s
                        seen.add(sp)
                        break
    return refs


def six_frames(seq):
    """Fewest internal stop codons over the six frames, and that frame's length in codons."""
    best = None
    for strand, s in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        for f in range(3):
            coding = s[f:]
            coding = coding[:len(coding) - len(coding) % 3].replace("N", "A")
            prot = str(Seq(coding).translate(table=GENETIC_CODE))
            stops = prot[:-1].count("*")  # a stop in the last codon may be the gene's own
            if best is None or stops < best[0]:
                best = (stops, f"{strand}{f + 1}", len(prot))
    return best


def main():
    p = argparse.ArgumentParser(description="Recheck flagged sequences against a panel of references")
    p.add_argument("--results", default=os.path.join("results", "genus"))
    p.add_argument("--data", default="data")
    p.add_argument("--reference", default=os.path.join("data", "reference", "coi_reference.fasta"))
    p.add_argument("--min-length", type=int, default=300)
    p.add_argument("--max-ambiguous", type=float, default=0.01)
    p.add_argument("--min-identity", type=float, default=0.75)
    args = p.parse_args()

    refs = panel(args.data, args.reference)
    merged = {r["record_id"]: r for r in csv.DictReader(open(os.path.join(args.results, "merged_records.csv")))}
    qc = list(csv.DictReader(open(os.path.join(args.results, "qc_records.csv"))))
    targets = [r for r in qc if r["qc_status"] == "FAIL"
               or (r["alignment_frameshift_gaps"] not in ("", "0") and r["qc_status"] == "PASS")]
    aligner = make_aligner()
    rows = []
    for r in targets:
        rec = merged[r["record_id"]]
        best = None
        for label, ref in refs.items():
            out, _ = check(aligner, ref, rec, args)
            score = (out.get("identity") or 0, -len(out.get("qc_flags", "")))
            if best is None or score > best[0]:
                best = (score, label, out)
        _, label, out = best
        stops, frame, codons = six_frames(rec["sequence"].upper())
        rows.append({"record_id": r["record_id"], "species_label": r["species_label"],
                     "flags_main_reference": r["qc_flags"], "status_main_reference": r["qc_status"],
                     "odd_gaps_main_reference": r["alignment_frameshift_gaps"],
                     "best_panel_reference": label, "identity_best": out.get("identity", ""),
                     "flags_best_reference": out.get("qc_flags", ""), "status_best_reference": out.get("qc_status", ""),
                     "odd_gaps_best_reference": out.get("alignment_frameshift_gaps", ""),
                     "six_frame_min_internal_stops": stops, "six_frame_best": frame})
    with open(os.path.join(args.results, "qc_recheck.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    fails = [r for r in rows if r["status_main_reference"] == "FAIL"]
    odd = [r for r in rows if r["status_main_reference"] == "PASS"]
    flag_main = Counter(f for r in fails for f in r["flags_main_reference"].split(";") if f)
    flag_best = Counter(f for r in fails for f in r["flags_best_reference"].split(";") if f)
    lines = [f"Panel: {len(refs)} references ({', '.join(refs)})",
             f"Failed sequences rechecked: {len(fails)}; still failing with the best panel reference: "
             f"{sum(r['status_best_reference'] == 'FAIL' for r in fails)}",
             "Flags, main reference -> best panel reference:"]
    for fl in sorted(set(flag_main) | set(flag_best)):
        lines.append(f"  {fl}: {flag_main[fl]} -> {flag_best[fl]}")
    fs = [r for r in fails if "FRAMESHIFT" in r["flags_main_reference"] or "STOP_CODON" in r["flags_main_reference"]]
    lines.append(f"Sequences flagged for frameshift or stop codon: {len(fs)}; with an open frame (no internal stop) "
                 f"in some reading frame: {sum(r['six_frame_min_internal_stops'] == 0 for r in fs)}")
    lines.append(f"Passing sequences with odd-length alignment gaps: {len(odd)}; still with such gaps against the best "
                 f"panel reference: {sum(r['odd_gaps_best_reference'] not in ('', '0', 0) for r in odd)}; "
                 f"failing against it: {sum(r['status_best_reference'] == 'FAIL' for r in odd)}")
    summary = "\n".join(lines)
    with open(os.path.join(args.results, "qc_recheck_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
