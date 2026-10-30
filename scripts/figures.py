"""Figures 1-5 of the manuscript, drawn from the result files.

Figure 1  steps of the pipeline
Figure 2  clusters holding the records named T. cinnabarinus
Figure 3  identification of anchor specimens, raw against cleaned database
Figure 4  maximum-likelihood trees of the two halves of COX1 (all records and
          the same records in both halves), with the names carried by each tip
Figure 5  documentation of records by submission period

Each figure is written as PDF (vector) and PNG (600 dpi) to results/figures.

Usage: python scripts/figures.py
"""
import csv
import os
import re
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
from Bio import Phylo  # noqa: E402

OUT = os.path.join("results", "figures")
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical slots 1-4 (adjacent order)
SEQ = ["#ffffff", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.edgecolor": BASE,
                     "axes.linewidth": 0.6, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none",
                     "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans",
                     "mathtext.it": "Liberation Sans:italic", "mathtext.bf": "Liberation Sans:bold"})


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, f"{name}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def italic_species(name):
    """'Tetranychus urticae' -> mathtext italic 'T. urticae'."""
    w = name.split()
    return f"$\\mathit{{T.\\,{w[1]}}}$" if len(w) == 2 and w[0] == "Tetranychus" else name


# Figure 1 ------------------------------------------------------------------
def figure1():
    fig, ax = plt.subplots(figsize=(6.7, 3.3))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.axis("off")

    def box(x, y, w, h, text, fill="#f4f3ef"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    fc=fill, ec=BASE, lw=0.6))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7, color=INK, linespacing=1.25)
        return (x, y, w, h)

    def arrow(a, b, start="r", end="l"):
        (x1, y1, w1, h1), (x2, y2, w2, h2) = a, b
        p = {"r": (x1 + w1 + 0.4, y1 + h1 / 2), "b": (x1 + w1 / 2, y1 - 0.4), "t": (x1 + w1 / 2, y1 + h1 + 0.4)}[start]
        q = {"l": (x2 - 0.4, y2 + h2 / 2), "t": (x2 + w2 / 2, y2 + h2 + 0.4), "b": (x2 + w2 / 2, y2 - 0.4)}[end]
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=7, lw=0.7, color=INK2,
                                     shrinkA=0, shrinkB=0))

    gb = box(1, 38, 15, 7, "GenBank COI\n(NCBI Entrez)")
    bold = box(1, 27, 15, 7, "BOLD COI\n(BCDM files)")
    its = box(1, 8, 15, 7, "GenBank ITS2")
    merge = box(22, 32, 15, 8, "Merge and\nmetadata audit\n(0-7 score)")
    qc = box(43, 32, 15, 8, "Sequence QC\n(length, frame,\nidentity)")
    clus = box(64, 32, 16, 8, "Haplotypes, K2P,\nclusters,\nbarcode gaps")
    anc = box(85, 38, 14, 7, "Anchor check\n(12 sources)")
    ident = box(85, 27.5, 14, 7.5, "Identification\ntest (leave one\nstudy out)")
    phy = box(64, 18, 16, 8, "Trees: 5′, 3′ and\nsame records\nin both halves")
    sens = box(85, 17, 14, 7.5, "Sensitivity:\nthresholds,\nanchors, folds")
    its2 = box(22, 8, 15, 7, "ITS2 distances,\nnearest names")
    pair = box(43, 8, 15, 7, "COI–ITS2 pairs\n(one-to-one links)")
    trend = box(22, 19.5, 15, 7, "Trends by year\nand country")
    arrow(gb, merge)
    arrow(bold, merge)
    arrow(merge, qc)
    arrow(qc, clus)
    arrow(clus, anc)
    arrow(clus, ident)
    arrow(anc, ident, "b", "t")
    arrow(clus, phy, "b", "t")
    arrow(ident, sens, "b", "t")
    arrow(its, its2)
    arrow(its2, pair)
    arrow(qc, pair, "b", "t")
    arrow(merge, trend, "b", "t")
    save(fig, "Figure1")


# Figure 2 ------------------------------------------------------------------
def figure2():
    checks = read("results/species/record_checks.csv")
    anchored = {r["cluster"]: r["anchored_species"] for r in read("results/anchors/cluster_anchors.csv")}
    comp = {r["cluster"]: r["composition"] for r in read("results/species/clusters.csv")}
    c = Counter(r["cluster"] for r in checks if r["species_name"] == "Tetranychus cinnabarinus")

    def label(cid):
        if anchored.get(cid):
            names = [italic_species(n) for n in anchored[cid].split("; ")]
            return f"{cid}: anchored as " + " + ".join(names)
        parts = [p.rsplit(" (", 1)[0] for p in comp[cid].split("; ")]
        named = [p for p in parts if p not in ("(unnamed)", "Tetranychus cinnabarinus")]
        lead = "mostly unnamed" if parts[0] == "(unnamed)" else "no anchor"
        return f"{cid}: {lead}" + (f", with {italic_species(named[0])}" if named else "")

    rows = c.most_common()[::-1]
    fig, ax = plt.subplots(figsize=(5.2, 2.3))
    y = range(len(rows))
    ax.barh(list(y), [n for _, n in rows], height=0.55, color=SERIES[0], edgecolor="none")
    for i, (_, n) in enumerate(rows):
        ax.text(n + 0.4, i, str(n), va="center", fontsize=7, color=INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels([label(cid) for cid, _ in rows], fontsize=7)
    ax.set_xlabel("Records named " + italic_species("Tetranychus cinnabarinus") + " (59 in total)")
    ax.set_xlim(0, max(n for _, n in rows) * 1.15)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, "Figure2")


# Figure 3 ------------------------------------------------------------------
STUDY_LABEL = {"matsuda2013": "Matsuda et al. 2013", "mendonca2011": "de Mendonça et al. 2011",
               "inak2022": "İnak et al. 2022", "srinivasa2012": "Srinivasa et al. 2012",
               "navajas_boursot2003": "Navajas & Boursot 2003", "abdolahadi2026": "Abdolahadi et al. 2025",
               "chen2014": "Chen et al. 2014", "li2010": "Li et al. 2010", "gotoh2009": "Gotoh et al. 2009"}
OUTCOMES = ["correct", "wrong", "ambiguous", "no match"]


def figure3():
    qs = read("results/robustness/queries_main.csv")
    trans = read("results/robustness/transitions.csv")
    order = ["matsuda2013", "mendonca2011", "inak2022", "srinivasa2012", "navajas_boursot2003",
             "abdolahadi2026", "chen2014", "li2010", "gotoh2009"]
    by = defaultdict(list)
    for q in qs:
        by[q["study_id"]].append(q)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.9, 3.4), gridspec_kw={"width_ratios": [1.65, 1]})
    ypos, labels, y = [], [], 0.0
    for sid in order + ["all"]:
        rs = qs if sid == "all" else by[sid]
        for key, tag in (("raw_result", "raw"), ("cleaned_result", "cleaned")):
            c = Counter(r[key] for r in rs)
            left = 0
            for k, o in enumerate(OUTCOMES):
                if c[o]:
                    ax.barh(y, c[o] / len(rs) * 100, left=left, height=0.62, color=SERIES[k],
                            edgecolor="white", linewidth=0.8)
                    left += c[o] / len(rs) * 100
            ypos.append(y)
            name = "All studies" if sid == "all" else STUDY_LABEL[sid]
            labels.append(f"{name} (n = {len(rs)}), {tag}" if tag == "raw" else "cleaned")
            y += 0.75
        y += 0.55
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=6.3)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Anchor queries (%)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=SERIES[k]) for k in range(4)]
    ax.legend(handles, ["Correct", "Wrong", "Ambiguous", "No match"], ncol=4, fontsize=6.5, frameon=False,
              loc="lower left", bbox_to_anchor=(-0.02, 1.0), handlelength=1.0, columnspacing=1.0)
    ax.text(-0.02, 1.09, "a", transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")

    m = [[int(r[o]) for o in OUTCOMES] for r in trans if r["raw"] in OUTCOMES]
    vmax = max(max(r) for r in m)
    for i, row in enumerate(m):
        for j, v in enumerate(row):
            level = 0 if v == 0 else 1 + int((len(SEQ) - 2) * (v ** 0.5) / (vmax ** 0.5))
            bx.add_patch(plt.Rectangle((j, i), 0.96, 0.96, color=SEQ[min(level, len(SEQ) - 1)], lw=0))
            bx.text(j + 0.48, i + 0.5, str(v), ha="center", va="center", fontsize=7.5,
                    color="white" if level >= 5 else INK)
    bx.set_xlim(0, 4)
    bx.set_ylim(4, 0)
    bx.set_xticks([k + 0.48 for k in range(4)])
    bx.set_xticklabels(["Correct", "Wrong", "Ambig.", "No match"], fontsize=6.8)
    bx.set_yticks([k + 0.48 for k in range(4)])
    bx.set_yticklabels(["Correct", "Wrong", "Ambig.", "No match"], fontsize=6.8)
    bx.set_xlabel("Cleaned database")
    bx.set_ylabel("Raw database")
    bx.xaxis.set_label_position("top")
    bx.xaxis.tick_top()
    bx.tick_params(length=0)
    for s in bx.spines.values():
        s.set_visible(False)
    bx.set_aspect("equal")
    bx.text(-0.30, 1.22, "b", transform=bx.transAxes, fontsize=10, fontweight="bold", va="bottom")
    fig.tight_layout(w_pad=2.5)
    save(fig, "Figure3")


# Figure 4 ------------------------------------------------------------------
COLUMNS = [("Tetranychus urticae", SERIES[0]), ("Tetranychus turkestani", SERIES[1]),
           ("Tetranychus truncatus", SERIES[2]), ("Tetranychus kanzawai", SERIES[3]),
           ("other names", MUTED), ("unnamed", BASE), ("anchor", INK)]


def tip_names(path):
    out = {}
    for r in read(path):
        names = {}
        for part in r["species_names"].split("; "):
            m = re.match(r"(.+) \((\d+)\)$", part)
            if m:
                names[m.group(1)] = int(m.group(2))
        out[r["haplotype"]] = (names, bool(r["anchor_names"]))
    return out


def draw_tree(ax, treefile, tipsfile, title, show_support=True, lw=0.35):
    tree = Phylo.read(treefile, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    tips = tree.get_terminals()
    ypos = {t: i for i, t in enumerate(tips)}
    depths = tree.depths()
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)

    def y_of(c):
        if c in ypos:
            return ypos[c]
        ys = [y_of(k) for k in c.clades]
        ypos[c] = (min(ys) + max(ys)) / 2
        return ypos[c]

    y_of(tree.root)
    segs = []
    for c in tree.find_clades(order="preorder"):
        if c.clades:
            segs.append([(depths[c], ypos[c.clades[0]]), (depths[c], ypos[c.clades[-1]])])
            for k in c.clades:
                segs.append([(depths[c], ypos[k]), (depths[k], ypos[k])])
    ax.add_collection(LineCollection(segs, colors=INK2, linewidths=lw))
    xmax = max(depths.values())
    if show_support:
        for c in tree.get_nonterminals():
            sup = c.confidence
            if sup is None and c.name:
                try:
                    sup = float(c.name.split("/")[-1])
                except ValueError:
                    sup = None
            if sup is not None and sup >= 95 and len(c.get_terminals()) >= 3:
                ax.plot(depths[c], ypos[c], "o", ms=1.2, color=INK, mec="none")
    names = tip_names(tipsfile)
    x0, dx = xmax * 1.04, xmax * 0.055
    for i, t in enumerate(tips):
        sp, anchor = names.get(t.name, ({}, False))
        flags = [sp.get(COLUMNS[0][0]), sp.get(COLUMNS[1][0]), sp.get(COLUMNS[2][0]), sp.get(COLUMNS[3][0]),
                 any(k not in ("(unnamed)",) + tuple(n for n, _ in COLUMNS[:4]) for k in sp),
                 "(unnamed)" in sp, anchor]
        for j, on in enumerate(flags):
            if on:
                ax.add_patch(plt.Rectangle((x0 + j * dx, ypos[t] - 0.5), dx * 0.8, 1.0, color=COLUMNS[j][1],
                                           lw=0))
    for j, (n, _) in enumerate(COLUMNS):
        lab = italic_species(n) if n.startswith("Tetranychus") else n
        ax.text(x0 + j * dx + dx * 0.4, -1.5, lab, rotation=90, ha="center", va="bottom", fontsize=6.2, color=INK)
    ax.set_xlim(-xmax * 0.01, x0 + len(COLUMNS) * dx)
    ax.set_ylim(len(tips) + 0.5, -0.5)
    ax.axis("off")
    ax.set_title(title, fontsize=7.5, loc="left", color=INK, pad=50)
    ybar = -0.035 * len(tips) - 1.5
    ax.plot([0, 0.05], [ybar, ybar], color=INK2, lw=0.6, clip_on=False)
    ax.text(0.058, ybar, "0.05 substitutions/site", ha="left", va="center", fontsize=6, color=INK2)


def figure4():
    fig = plt.figure(figsize=(6.9, 9.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.25, 1], hspace=0.30, wspace=0.12)
    panels = [("results/phylogeny/5prime", "a  5′ half, all records (505 haplotypes)"),
              ("results/phylogeny/3prime", "b  3′ half, all records (422 haplotypes)"),
              ("results/phylogeny/paired/5prime", "c  5′ half, the same 98 records"),
              ("results/phylogeny/paired/3prime", "d  3′ half, the same 98 records")]
    for k, (d, title) in enumerate(panels):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        draw_tree(ax, os.path.join(d, "iqtree.treefile"), os.path.join(d, "tips.csv"), title,
                  lw=0.3 if k < 2 else 0.5)
    save(fig, "Figure4")


# Figure 5 ------------------------------------------------------------------
def figure5():
    rows = read("results/trends/by_period.csv")
    periods = [r["period"].replace("before 2006", "Before 2006").replace("2021 onwards", "2021 on")
               .replace("-", "–") for r in rows]
    series = [("Species name", [100 - float(r["unnamed_pct"]) for r in rows]),
              ("Voucher", [float(r["voucher_pct"]) for r in rows]),
              ("Coordinates", [float(r["coordinates_pct"]) for r in rows]),
              ("Publication cited", [float(r["publication_pct"]) for r in rows])]
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    x = list(range(len(periods)))
    for k, (name, ys) in enumerate(series):
        ax.plot(x, ys, color=SERIES[k], lw=1.4, marker="o", ms=3.6, mec="white", mew=0.8, label=name,
                solid_capstyle="round", solid_joinstyle="round")
    ends = sorted(((ys[-1], name) for name, ys in series))
    placed = []
    for v, name in ends:
        yv = v
        while any(abs(yv - p) < 6 for p in placed):
            yv += 6
        placed.append(yv)
        ax.annotate(name, (x[-1], v), xytext=(x[-1] + 0.18, yv), fontsize=6.8, color=INK, va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED) if abs(yv - v) > 0.5 else None)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p}\nn = {r['records']}" for p, r in zip(periods, rows)], fontsize=7)
    ax.set_xlim(-0.2, len(periods) - 0.2 + 1.1)
    ax.set_ylim(0, 105)
    ax.set_ylabel("Records (%)")
    ax.set_xlabel("Submission period")
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.28), ncol=4,
              handlelength=1.6)
    save(fig, "Figure5")


def main():
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    print(f"Figures written to {OUT}")


if __name__ == "__main__":
    main()
