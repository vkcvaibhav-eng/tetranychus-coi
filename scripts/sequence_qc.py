"""Module 3: quality check of every COI sequence against a reference COX1
coding sequence (default: data/reference/coi_reference.fasta).

For each sequence it records orientation, the part of COI covered, identity
to the reference, ambiguous bases, internal stop codons (invertebrate
mitochondrial code, table 5), and indels that break the reading frame.
Stop codons and frameshifts point to sequencing errors or nuclear
pseudogenes (NUMTs). A sequence flagged in the reference frame is checked
again by translation in all six frames: if one frame is open and its protein
matches the reference protein (identity >= 0.70 over >= 90% of its length),
the flag is treated as an alignment artefact and dropped.

Reads:  results/merged_records.csv
Writes: results/qc_records.csv, results/qc_passed.fasta, results/qc_summary.txt

Usage: python scripts/sequence_qc.py [--reference FILE] [--results DIR]
       [--min-length 300] [--max-ambiguous 0.01] [--min-identity 0.75]
"""
import argparse
import csv
import os
import sys
from collections import Counter

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

GENETIC_CODE = 5  # invertebrate mitochondrial
QC_FIELDS = ["record_id", "source", "species_label", "seq_length", "orientation",
             "ref_start", "ref_end", "ref_coverage", "cox1_half", "identity",
             "query_aligned_fraction", "ambiguous_bases", "ambiguous_fraction",
             "reading_frame_intact", "internal_stops", "alignment_frameshift_gaps",
             "translation_check", "qc_flags", "qc_status"]
MIN_PROTEIN_IDENTITY = 0.70   # translation-guided check: open frame must match the reference protein
MIN_PROTEIN_COVER = 0.90      # over this share of the translated sequence
_PROT_ALIGNER = None


def translation_check(seq, ref):
    """Translation-guided check of a sequence flagged for a frameshift or stop codon.

    The flag comes from reading the nucleotide alignment in the reference
    frame, and a local alignment against a divergent reference can place
    compensating gaps that misplace that frame. Here the sequence is read in
    all six frames; if one frame has no internal stop codon and its protein
    matches the reference protein over most of its length, the sequence
    encodes an intact COX1 stretch and the flag is an alignment artefact.
    Returns (verdict, protein identity, frame).
    """
    global _PROT_ALIGNER
    if _PROT_ALIGNER is None:
        _PROT_ALIGNER = PairwiseAligner()
        _PROT_ALIGNER.mode = "local"
        _PROT_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
        _PROT_ALIGNER.open_gap_score, _PROT_ALIGNER.extend_gap_score = -10, -1
    ref_prot = str(Seq(ref[:len(ref) - len(ref) % 3]).translate(table=GENETIC_CODE)).rstrip("*")
    best = ("no open frame", 0.0, "")
    for strand, s in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        for f in range(3):
            coding = s[f:]
            coding = coding[:len(coding) - len(coding) % 3]
            if len(coding) < 90:
                continue
            prot = str(Seq(coding).translate(table=GENETIC_CODE))
            if "*" in prot[:-1]:
                continue
            prot = prot.rstrip("*")
            aln = _PROT_ALIGNER.align(ref_prot, prot)[0]
            span = sum(qe - qs for (qs, qe) in aln.aligned[1])
            same = sum(ref_prot[ts + i] == prot[qs + i]
                       for (ts, te), (qs, qe) in zip(*aln.aligned) for i in range(te - ts))
            ident = same / span if span else 0.0
            if span >= MIN_PROTEIN_COVER * len(prot) and ident >= MIN_PROTEIN_IDENTITY and ident > best[1]:
                best = ("open COX1 frame", ident, f"{strand}{f + 1}")
    return best


def make_aligner():
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 2
    a.mismatch_score = -3
    a.open_gap_score = -6
    a.extend_gap_score = -2
    # Ambiguous bases (N, R, Y ...) are scored as mismatches by default.
    return a


def best_alignment(aligner, ref, seq):
    """Align the sequence and its reverse complement; keep the better one."""
    best = None
    for orient, s in (("forward", seq), ("reverse", str(Seq(seq).reverse_complement()))):
        try:
            aln = next(iter(aligner.align(ref, s)))
        except StopIteration:
            continue
        if best is None or aln.score > best[1].score:
            best = (orient, aln, s)
    return best


def block_stats(ref, query, blocks_t, blocks_q, ref_len):
    matches = columns = stops = shifts = 0
    for i, ((ts, te), (qs, qe)) in enumerate(zip(blocks_t, blocks_q)):
        seg_r, seg_q = ref[ts:te], query[qs:qe]
        matches += sum(a == b for a, b in zip(seg_r, seg_q))
        columns += te - ts
        # Translate this block in the reference reading frame.
        offset = (3 - ts % 3) % 3
        coding = seg_q[offset:]
        coding = coding[:len(coding) - len(coding) % 3]
        if coding:
            protein = str(Seq(coding).translate(table=GENETIC_CODE))
            if te >= ref_len - 3 and protein.endswith("*"):
                protein = protein[:-1]  # terminal stop codon of the gene
            stops += protein.count("*")
        if i + 1 < len(blocks_t):
            rgap = blocks_t[i + 1][0] - te
            qgap = blocks_q[i + 1][0] - qe
            columns += rgap + qgap
            if (qgap - rgap) % 3 != 0:
                shifts += 1
    return matches / columns if columns else 0.0, stops, shifts


def frame_stops(oriented, ref_start, ref_end, q_start, q_end, ref_len):
    """Stop codons when the whole aligned stretch is read in one frame.

    A real indel that is not a multiple of three breaks this frame and
    produces stop codons downstream; an alignment artefact does not.
    """
    offset = (3 - ref_start % 3) % 3
    coding = oriented[q_start + offset:q_end]
    coding = coding[:len(coding) - len(coding) % 3]
    if not coding:
        return 0
    protein = str(Seq(coding).translate(table=GENETIC_CODE))
    if ref_end >= ref_len - 3 and protein.endswith("*"):
        protein = protein[:-1]
    return protein.count("*")


def check(aligner, ref, rec, args):
    seq = rec["sequence"].upper()
    out = {k: rec.get(k, "") for k in ("record_id", "source", "species_label")}
    out["seq_length"] = len(seq)
    ambiguous = sum(1 for c in seq if c not in "ACGT")
    out["ambiguous_bases"] = ambiguous
    out["ambiguous_fraction"] = round(ambiguous / len(seq), 4) if seq else ""
    flags = []
    if len(seq) < args.min_length:
        flags.append("SHORT")
    if seq and ambiguous / len(seq) > args.max_ambiguous:
        flags.append("AMBIGUOUS")

    best = best_alignment(aligner, ref, seq) if seq else None
    if best is None or best[1].score < 100:
        flags.append("NOT_ALIGNED")
        out.update(orientation="", qc_flags=";".join(flags), qc_status="FAIL")
        return out, None

    orient, aln, oriented = best
    blocks_t, blocks_q = aln.aligned
    ref_start, ref_end = int(blocks_t[0][0]), int(blocks_t[-1][1])
    q_start, q_end = int(blocks_q[0][0]), int(blocks_q[-1][1])
    identity, block_stops, shifts = block_stats(ref, oriented, blocks_t, blocks_q, len(ref))
    whole_stops = frame_stops(oriented, ref_start, ref_end, q_start, q_end, len(ref))

    out["orientation"] = orient
    out["ref_start"] = ref_start + 1
    out["ref_end"] = ref_end
    out["ref_coverage"] = round((ref_end - ref_start) / len(ref), 3)
    if (ref_end - ref_start) / len(ref) >= 0.9:
        out["cox1_half"] = "full"
    else:
        out["cox1_half"] = "5prime" if (ref_start + ref_end) / 2 < len(ref) / 2 else "3prime"
    out["identity"] = round(identity, 4)
    out["query_aligned_fraction"] = round((q_end - q_start) / len(seq), 3)
    out["reading_frame_intact"] = "yes" if whole_stops == 0 else "no"
    out["internal_stops"] = whole_stops
    out["alignment_frameshift_gaps"] = shifts

    if identity < args.min_identity:
        flags.append("LOW_IDENTITY")
    if (q_end - q_start) / len(seq) < 0.9:
        flags.append("PARTIAL_ALIGNMENT")
    # Flag only when the reading frame is really broken. Gaps of odd length
    # with an intact frame are alignment artefacts against a divergent reference.
    out["translation_check"] = ""
    if whole_stops:
        verdict, ident, frame = translation_check(oriented, ref)
        out["translation_check"] = f"{verdict} ({frame}, protein identity {ident:.2f})" if frame else verdict
        if verdict != "open COX1 frame":
            if shifts:
                flags.append("FRAMESHIFT")
            if block_stops or not shifts:
                flags.append("STOP_CODON")
    out["qc_flags"] = ";".join(flags)
    out["qc_status"] = "PASS" if not flags else "FAIL"
    return out, oriented


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--reference", default=os.path.join("data", "reference", "coi_reference.fasta"))
    p.add_argument("--results", default="results")
    p.add_argument("--min-length", type=int, default=300)
    p.add_argument("--max-ambiguous", type=float, default=0.01)
    p.add_argument("--min-identity", type=float, default=0.75)
    args = p.parse_args()

    if not os.path.exists(args.reference):
        sys.exit(f"Reference not found: {args.reference}. Run scripts/fetch_reference.py first.")
    ref = str(SeqIO.read(args.reference, "fasta").seq).upper()
    ref_protein = str(Seq(ref[:len(ref) - len(ref) % 3]).translate(table=GENETIC_CODE))
    if "*" in ref_protein.rstrip("*"):
        sys.exit("The reference itself has internal stop codons; check the reference and genetic code.")

    with open(os.path.join(args.results, "merged_records.csv"), newline="") as f:
        records = list(csv.DictReader(f))

    aligner = make_aligner()
    rows, passed = [], []
    for rec in records:
        row, oriented = check(aligner, ref, rec, args)
        rows.append(row)
        if row["qc_status"] == "PASS":
            passed.append((row["record_id"], oriented))

    with open(os.path.join(args.results, "qc_records.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=QC_FIELDS)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(args.results, "qc_passed.fasta"), "w") as f:
        for rid, s in passed:
            f.write(f">{rid}\n{s}\n")

    n = len(rows)
    flag_counts = Counter(fl for r in rows for fl in r["qc_flags"].split(";") if fl)
    halves = Counter(r.get("cox1_half", "") for r in rows if r.get("cox1_half"))
    lines = [f"Reference COX1: {len(ref)} bp", f"Records checked: {n}",
             f"Passed: {len(passed)} ({100 * len(passed) / n:.1f}%)" if n else "Passed: 0",
             "Flags (a record can have several):"]
    lines += [f"  {k}: {v} ({100 * v / n:.1f}%)" for k, v in flag_counts.most_common()]
    lines += ["Region of COX1 covered:"] + [f"  {k}: {v}" for k, v in sorted(halves.items())]
    summary = "\n".join(lines)
    with open(os.path.join(args.results, "qc_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
