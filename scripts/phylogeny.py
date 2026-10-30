"""Module 8: phylogenetic check of the distance-based results.

For each half of COX1 (5' and 3'), sequences that passed the quality check
and cover most of that half are cut to the region on reference coordinates,
collapsed to unique haplotypes, aligned with MAFFT and analysed with IQ-TREE
(maximum likelihood, ultrafast bootstrap). The tree is midpoint-rooted and
checked for:

- for each species name (after published synonymies): the share of its
  records in its exclusive clade and in its permissive group (below), with
  bootstrap support;
- the same for anchor specimens only (morphologically identified).

Two measures are reported for each species. The exclusive clade is the
largest clade whose tips carry only that name (unnamed tips allowed); a
haplotype carried by records under several names conflicts with every one of
them. The permissive group lets such shared haplotypes count towards each of
their names, so it measures how concentrated a species is, not whether it is
exclusive; the number of records of other names inside it is reported.
Because midpoint rooting is arbitrary
without an outgroup, the largest group is also found without regard to the
root: any group of tips cut off by one branch of the tree (a clade, or the
complement of a clade) is accepted.

With --paired, only records that cover both windows are used, and the two
trees are built from the same specimens, so that differences between the
halves are not differences in sampling.

Requires the mafft and iqtree2 (or iqtree) programs. With --reuse, existing
tree files are evaluated again without rebuilding them.

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


def pure_clades(clade, member, neutral):
    """Maximal clades whose tips are all members or neutral, with at least one member."""
    tips = clade.get_terminals()
    if all(member(t.name) or neutral(t.name) for t in tips):
        return [clade] if any(member(t.name) for t in tips) else []
    if clade.is_terminal():
        return []
    return [c for child in clade.clades for c in pure_clades(child, member, neutral)]


def largest_split_group(tree, member, neutral, weight):
    """Largest pure group cut off by one branch: a clade or its complement."""
    all_tips = tree.get_terminals()
    total_names = [t.name for t in all_tips]
    best_n, best_sup = 0, None
    for clade in tree.find_clades():
        inside = {t.name for t in clade.get_terminals()}
        for group, node in ((inside, clade), (set(total_names) - inside, clade)):
            if not group or not all(member(t) or neutral(t) for t in group):
                continue
            n = sum(weight(t) for t in group)
            if n > best_n:
                best_n = n
                best_sup = None if node.is_terminal() else support(node)
    return best_n, best_sup


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
    p.add_argument("--seed", type=int, default=None,
                   help="IQ-TREE random seed (the seed used is always written to the IQ-TREE log)")
    p.add_argument("--fast", action="store_true", help="IQ-TREE -fast, no bootstrap (for testing)")
    p.add_argument("--paired", action="store_true",
                   help="Use only records covering both windows (written to <out>/paired)")
    p.add_argument("--reuse", action="store_true", help="Evaluate existing tree files without rebuilding")
    args = p.parse_args()
    if args.paired:
        args.out = os.path.join(args.out, "paired")
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

    windows = {h: region_window(qc_rows, h) for h in ("5prime", "3prime")}
    eligible = None
    if args.paired:
        eligible = set()
        for rid, proj in projected.items():
            if all(w and sum(c != "-" for c in proj[w[0]:w[1]]) / max(w[1] - w[0], 1) >= args.min_cover
                   for w in windows.values()):
                eligible.add(rid)

    mono_rows, lines = [], []
    if args.paired:
        lines.append(f"Paired mode: {len(eligible)} records cover at least {args.min_cover:.0%} of both windows")
    for half in ("5prime", "3prime"):
        win = windows[half]
        if not win:
            continue
        w0, w1 = win
        rdir = os.path.join(args.out, half)
        os.makedirs(rdir, exist_ok=True)
        haps, members = {}, defaultdict(list)
        for rid, proj in projected.items():
            if eligible is not None and rid not in eligible:
                continue
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
        prefix = os.path.join(rdir, "iqtree")
        if not (args.reuse and os.path.exists(prefix + ".treefile")):
            with open(aln, "w") as f:
                subprocess.run(["mafft", "--auto", "--quiet", "--thread", "-1", hap_fa], stdout=f, check=True)
            cmd = [iqtree, "-s", aln, "-m", "GTR+F+G4", "-T", str(args.threads), "--prefix", prefix,
                   "-redo", "--quiet"]
            cmd += ["-fast"] if args.fast else ["-B", str(args.bootstrap)]
            if args.seed is not None:
                cmd += ["-seed", str(args.seed)]
            subprocess.run(cmd, check=True)
        aligned = list(SeqIO.parse(aln, "fasta"))
        cols = len(aligned[0].seq) if aligned else 0
        internal = set()  # columns with a gap inside a sequence (not at its ends)
        for r in aligned:
            sq = str(r.seq)
            a, b = len(sq) - len(sq.lstrip("-")), len(sq.rstrip("-"))
            internal.update(i for i in range(a, b) if sq[i] == "-")
        gap_cols = len(internal)

        tree = Phylo.read(prefix + ".treefile", "newick")
        tree.root_at_midpoint()
        lines.append(f"{half}: reference window {w0 + 1}-{w1}, {sum(len(v) for v in members.values())} "
                     f"sequences, {len(haps)} haplotypes, model GTR+F+G4"
                     + (", fast search" if args.fast else f", {args.bootstrap} ultrafast bootstraps")
                     + f"; alignment {cols} columns, {gap_cols} with an internal gap")

        for label, tipnames in (("all records", names), ("anchors only", anchors_in)):
            species = sorted({n for c in tipnames.values() for n in c if n})
            for sp in species:
                def strict(t, sp=sp, tipnames=tipnames):
                    # an exclusive member carries this species name and no other name
                    named = {n for n in tipnames[t] if n}
                    return named == {sp}

                def shared(t, sp=sp, tipnames=tipnames):
                    # permissive: any tip carrying the name, even if it also carries others
                    return bool(tipnames[t].get(sp))

                def neutral(t, tipnames=tipnames):
                    # Unnamed tips (all records) or tips without anchors (anchors only)
                    # do not contradict a species.
                    return not any(n for n in tipnames[t])

                def weight(t, sp=sp, tipnames=tipnames):
                    return tipnames[t].get(sp, 0)

                total = sum(c.get(sp, 0) for c in tipnames.values())
                n_haps = sum(1 for c in tipnames.values() if c.get(sp))
                if n_haps < 2:
                    continue
                on_mixed = sum(c.get(sp, 0) for c in tipnames.values()
                               if c.get(sp) and len({n for n in c if n}) > 1)

                def largest(member):
                    clades = pure_clades(tree.root, member, neutral)
                    sized = sorted(((sum(weight(t.name) for t in c.get_terminals()), c) for c in clades),
                                   key=lambda x: -x[0])
                    if not sized:
                        return 0, None, 0, None
                    n_main, main = sized[0]
                    sup = support(main) if not main.is_terminal() else None
                    return n_main, sup, len(sized), main

                ex_n, ex_sup, ex_k, _ = largest(strict)
                pm_n, pm_sup, pm_k, pm_main = largest(shared)
                other_in = sum(sum(v for n, v in tipnames[t.name].items() if n and n != sp)
                               for t in pm_main.get_terminals()) if pm_main is not None else 0
                rf_n, rf_sup = largest_split_group(tree, strict, neutral, weight)
                mono_rows.append({"region": half, "basis": label, "species": sp,
                                  "records": total, "haplotypes": n_haps,
                                  "records_on_shared_haplotypes": on_mixed,
                                  "exclusive_clade_share_pct": f"{100 * ex_n / total:.1f}",
                                  "exclusive_clade_support": "" if ex_sup is None else f"{ex_sup:.0f}",
                                  "exclusive_clades": ex_k,
                                  "monophyletic_exclusive": "yes" if ex_n == total else "no",
                                  "root_free_exclusive_share_pct": f"{100 * rf_n / total:.1f}",
                                  "concentration_share_pct": f"{100 * pm_n / total:.1f}",
                                  "concentration_support": "" if pm_sup is None else f"{pm_sup:.0f}",
                                  "other_named_records_in_concentration_group": other_in})

    with open(os.path.join(args.out, "monophyly.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mono_rows[0].keys()) if mono_rows else ["region"])
        w.writeheader()
        w.writerows(mono_rows)

    for half in ("5prime", "3prime"):
        for basis in ("anchors only", "all records"):
            rs = [r for r in mono_rows if r["region"] == half and r["basis"] == basis]
            if not rs:
                continue
            lines.append(f"{half}, {basis}: {sum(r['monophyletic_exclusive'] == 'yes' for r in rs)} of {len(rs)} "
                         f"species form one exclusive clade (unnamed tips allowed, shared haplotypes not)")
            for r in sorted(rs, key=lambda r: -r["records"])[:14]:
                lines.append(f"  {r['species']} ({r['records']} records, {r['haplotypes']} haplotypes): "
                             f"exclusive clade {r['exclusive_clade_share_pct']}%"
                             + (f" (support {r['exclusive_clade_support']})" if r["exclusive_clade_support"] else "")
                             + f", root-free {r['root_free_exclusive_share_pct']}%; "
                             f"{r['records_on_shared_haplotypes']} records on shared haplotypes; "
                             f"permissive group {r['concentration_share_pct']}%"
                             + (f" (support {r['concentration_support']})" if r["concentration_support"] else "")
                             + f" with {r['other_named_records_in_concentration_group']} records of other names")
    summary = "\n".join(lines)
    with open(os.path.join(args.out, "phylogeny_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
