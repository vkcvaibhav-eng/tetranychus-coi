"""Module 8: phylogenetic check of the distance-based results.

For each half of COX1 (5' and 3'), sequences that passed the quality check
and cover most of that half are cut to the region on reference coordinates,
collapsed to unique haplotypes, aligned with MAFFT and analysed with IQ-TREE
(maximum likelihood, ultrafast bootstrap). The tree is midpoint-rooted and
checked for:

- monophyly of each species name (after published synonymies), with the
  bootstrap support of its clade and the number of separate clades
  ("fragments") its haplotypes fall into;
- the same for anchor specimens only (morphologically identified).

Requires the mafft and iqtree2 (or iqtree) programs.

Reads:  results/genus/{merged_records.csv,qc_passed.fasta,qc_records.csv},
        data/anchors/*, reference COX1
Writes: results/phylogeny/<region>/{haplotypes.fasta,alignment.fasta,
        tips.csv,iqtree files}, results/phylogeny/monophyly.csv,
        results/phylogeny/phylogeny_summary.txt

Usage: python scripts/phylogeny.py [--min-cover 0.8] [--bootstrap 1000]
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

from Bio import Phylo, SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_qc import make_aligner  # noqa: E402
from species_analysis import project, species_name  # noqa: E402


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def region_window(qc_rows, half):
    """Reference window covered by most sequences of one COX1 half."""
    rows = [r for r in qc_rows if r["qc_status"] == "PASS" and r["cox1_half"] == half]
    if not rows:
        return None
    starts = sorted(int(r["ref_start"]) - 1 for r in rows)
    ends = sorted(int(r["ref_end"]) for r in rows)
    return starts[len(starts) // 2], ends[len(ends) // 2]


def fragments(clade, is_member):
    """Number of maximal clades whose tips are all members."""
    tips = clade.get_terminals()
    if all(is_member(t.name) for t in tips):
        return 1
    if clade.is_terminal():
        return 0
    return sum(fragments(c, is_member) for c in clade.clades)


def support(node):
    if node.confidence is not None:
        return float(node.confidence)
    if node.name:
        try:
            return float(node.name.split("/")[-1])
        except ValueError:
            return None
    return None


def main():
    p = argparse.ArgumentParser(description="Phylogenetic check of COI species groups")
    p.add_argument("--genus", default=os.path.join("results", "genus"))
    p.add_argument("--anchors", default=os.path.join("data", "anchors"))
    p.add_argument("--reference", default=os.path.join("data", "reference", "coi_reference.fasta"))
    p.add_argument("--out", default=os.path.join("results", "phylogeny"))
    p.add_argument("--min-cover", type=float, default=0.8,
                   help="Minimum fraction of the region window a sequence must cover")
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--threads", default="AUTO")
    p.add_argument("--fast", action="store_true", help="IQ-TREE -fast, no bootstrap (for testing)")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    iqtree = shutil.which("iqtree2") or shutil.which("iqtree")
    if not shutil.which("mafft") or not iqtree:
        sys.exit("mafft and iqtree2 are required")

    syn_path = os.path.join(args.anchors, "synonyms.csv")
    synonyms = {r["name"]: r["accepted_name"] for r in read(syn_path)} if os.path.exists(syn_path) else {}

    def canon(label):
        n = species_name(label)
        return synonyms.get(n, n)

    merged = {r["record_id"]: r for r in read(os.path.join(args.genus, "merged_records.csv"))}
    qc_rows = read(os.path.join(args.genus, "qc_records.csv"))
    by_acc = {r["genbank_accession"]: rid for rid, r in merged.items() if r["genbank_accession"]}
    studies = {s["study_id"]: s for s in read(os.path.join(args.anchors, "anchor_studies.csv"))}
    anchor_name = {}
    for a in read(os.path.join(args.anchors, "anchor_records.csv")):
        rid = by_acc.get(a["accession"])
        if rid and studies[a["study_id"]]["evidence_level"] in "AB":
            anchor_name[rid] = canon(merged[rid]["species_label"])

    ref = str(SeqIO.read(args.reference, "fasta").seq).upper()
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(args.genus, "qc_passed.fasta"), "fasta")}
    aligner = make_aligner()
    projected = {rid: project(aligner, ref, s) for rid, s in seqs.items()}

    mono_rows, lines = [], []
    for half in ("5prime", "3prime"):
        win = region_window(qc_rows, half)
        if not win:
            continue
        w0, w1 = win
        rdir = os.path.join(args.out, half)
        os.makedirs(rdir, exist_ok=True)
        haps, members = {}, defaultdict(list)
        for rid, proj in projected.items():
            seg = proj[w0:w1]
            covered = sum(c != "-" for c in seg) / max(len(seg), 1)
            if covered < args.min_cover:
                continue
            raw = seg.replace("-", "")
            if raw not in haps:
                haps[raw] = f"H{len(haps) + 1:04d}"
            members[haps[raw]].append(rid)
        if len(haps) < 4:
            continue
        hap_fa = os.path.join(rdir, "haplotypes.fasta")
        with open(hap_fa, "w") as f:
            for s, hid in haps.items():
                f.write(f">{hid}\n{s}\n")
        names = {hid: Counter(canon(merged[r]["species_label"]) for r in recs) for hid, recs in members.items()}
        anchors_in = {hid: Counter(anchor_name[r] for r in recs if r in anchor_name) for hid, recs in members.items()}
        with open(os.path.join(rdir, "tips.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["haplotype", "records", "species_names", "anchor_names", "record_ids"])
            for hid, recs in members.items():
                w.writerow([hid, len(recs), "; ".join(f"{n or '(unnamed)'} ({k})" for n, k in names[hid].most_common()),
                            "; ".join(f"{n} ({k})" for n, k in anchors_in[hid].most_common()), " ".join(recs)])

        aln = os.path.join(rdir, "alignment.fasta")
        with open(aln, "w") as f:
            subprocess.run(["mafft", "--auto", "--quiet", "--thread", "-1", hap_fa], stdout=f, check=True)
        prefix = os.path.join(rdir, "iqtree")
        cmd = [iqtree, "-s", aln, "-m", "GTR+F+G4", "-T", str(args.threads), "--prefix", prefix,
               "-redo", "--quiet"]
        cmd += ["-fast"] if args.fast else ["-B", str(args.bootstrap)]
        subprocess.run(cmd, check=True)

        tree = Phylo.read(prefix + ".treefile", "newick")
        tree.root_at_midpoint()
        lines.append(f"{half}: reference window {w0 + 1}-{w1}, {sum(len(v) for v in members.values())} "
                     f"sequences, {len(haps)} haplotypes, model GTR+F+G4"
                     + (", fast search" if args.fast else f", {args.bootstrap} ultrafast bootstraps"))

        for label, tipnames in (("all records", names), ("anchors only", anchors_in)):
            species = sorted({n for c in tipnames.values() for n in c if n})
            for sp in species:
                member_tips = [t for t in tree.get_terminals() if tipnames[t.name].get(sp)]
                if len(member_tips) < 2:
                    continue
                mrca = tree.common_ancestor(member_tips)
                clade_tips = mrca.get_terminals()
                mono = all(tipnames[t.name].get(sp) for t in clade_tips)
                n_frag = fragments(tree.root, lambda n: bool(tipnames[n].get(sp))) if not mono else 1
                others = Counter(o for t in clade_tips for o in names[t.name] if o and o != sp)
                mono_rows.append({"region": half, "basis": label, "species": sp,
                                  "haplotypes": len(member_tips), "monophyletic": "yes" if mono else "no",
                                  "clade_support": "" if not mono or support(mrca) is None else f"{support(mrca):.0f}",
                                  "fragments": n_frag,
                                  "other_names_in_mrca_clade": "; ".join(f"{o} ({k})" for o, k in others.most_common(5))})

    with open(os.path.join(args.out, "monophyly.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mono_rows[0].keys()) if mono_rows else ["region"])
        w.writeheader()
        w.writerows(mono_rows)

    for half in ("5prime", "3prime"):
        for basis in ("anchors only", "all records"):
            rs = [r for r in mono_rows if r["region"] == half and r["basis"] == basis]
            if not rs:
                continue
            lines.append(f"{half}, {basis}: {sum(r['monophyletic'] == 'yes' for r in rs)} of {len(rs)} "
                         f"species monophyletic")
            for r in sorted(rs, key=lambda r: -r["haplotypes"])[:14]:
                extra = f"support {r['clade_support']}" if r["monophyletic"] == "yes" else \
                    f"{r['fragments']} fragments; MRCA also holds {r['other_names_in_mrca_clade'] or '-'}"
                lines.append(f"  {r['species']} ({r['haplotypes']} haplotypes): "
                             f"{'monophyletic' if r['monophyletic'] == 'yes' else 'not monophyletic'}, {extra}")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "phylogeny_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
