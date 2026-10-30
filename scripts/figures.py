"""Figures 1-5 of the manuscript, drawn from the result files.

Figure 1  steps of the pipeline
Figure 2  clusters holding the records named T. cinnabarinus
Figure 3  identification of anchor specimens, raw against cleaned database
Figure 4  maximum-likelihood trees of the two halves of COX1 (all records and
          the same records in both halves), with the names carried by each tip
Figure 5  documentation of records by submission period
Figures S1-S4  the four trees of Figure 4 at full size, with labelled tips
              (also written as labelled Newick files next to each tree)

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
    ax.set_xlim(0, 101)
    ax.set_ylim(0, 50)
    ax.axis("off")

    def box(x, y, w, h, text, fill="#f4f3ef"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    fc=fill, ec=BASE, lw=0.6))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=6.4 if text.count("\n") >= 3 else 7,
                color=INK, linespacing=1.25)
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
    qc = box(43, 32, 15, 8, "Sequence QC\n(length, identity,\nframe; translation)")
    clus = box(64, 32, 16, 8, "Haplotypes, K2P,\nclusters,\nbarcode gaps")
    anc = box(82.4, 37, 18, 9.5, "Anchor check:\nTetranychus 9 studies,\nOligonychus 10\n(evidence per record)")
    ident = box(82.4, 25.5, 18, 7.5, "Identification test\n(leave one study out;\nfour reference sets)")
    phy = box(64, 18, 16, 8, "Trees: 5′, 3′ and\nsame records\nin both halves")
    sens = box(82.4, 14.5, 18, 7.5, "Sensitivity:\nthresholds, anchors,\nfolds, species pairs")
    its2 = box(22, 8, 15, 7, "ITS2 distances,\nnearest names")
    pair = box(43, 8, 15, 7, "COI–ITS2 links\n(one-to-one,\nby key type)")
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
OUTCOMES = ["correct", "wrong", "ambiguous", "nearest unnamed", "no match"]
OUT_LABEL = ["Correct", "Wrong", "Ambiguous", "Nearest record unnamed", "No match"]
OUT_COLOR = [SERIES[0], SERIES[1], SERIES[2], MUTED, SERIES[3]]
SET_LABEL = {"all": "all records", "named": "named records", "cleaned": "cleaned", "conservative": "conservative"}


def outcome_bar(ax, y, rs, key, height=0.62):
    c = Counter(r[key] for r in rs)
    left = 0
    for k, o in enumerate(OUTCOMES):
        if c[o]:
            ax.barh(y, c[o] / len(rs) * 100, left=left, height=height, color=OUT_COLOR[k],
                    edgecolor="white", linewidth=0.8)
            left += c[o] / len(rs) * 100


def outcome_legend(ax, y=1.0):
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in OUT_COLOR]
    ax.legend(handles, OUT_LABEL, ncol=5, fontsize=6.2, frameon=False, loc="lower left",
              bbox_to_anchor=(-0.02, y), handlelength=1.0, columnspacing=0.9)


def figure3():
    qs = read("results/robustness/queries_main.csv")
    order = ["matsuda2013", "mendonca2011", "inak2022", "srinivasa2012", "navajas_boursot2003",
             "abdolahadi2026", "chen2014", "li2010", "gotoh2009"]
    by = defaultdict(list)
    for q in qs:
        by[q["study_id"]].append(q)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.9, 4.4), gridspec_kw={"width_ratios": [1.7, 1]})
    ypos, labels, y = [], [], 0.0
    for sid in ["all"] + order:
        rs = qs if sid == "all" else by[sid]
        for key in ("all", "named", "cleaned"):
            outcome_bar(ax, y, rs, f"{key}_result")
            ypos.append(y)
            name = f"All studies ({len(order)})" if sid == "all" else STUDY_LABEL[sid]
            labels.append(f"{name}, {len(rs)} quer{'y' if len(rs) == 1 else 'ies'}: {SET_LABEL[key]}"
                          if key == "all" else SET_LABEL[key])
            y += 0.72
        y += 0.5
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=5.6)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Anchor queries (%)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    outcome_legend(ax)
    ax.text(-0.02, 1.07, "a", transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")

    outs = ["correct", "wrong", "ambiguous", "no match"]
    m = [[sum(1 for q in qs if q["named_result"] == a_ and q["cleaned_result"] == b_) for b_ in outs] for a_ in outs]
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
    bx.set_xlabel("Cleaned records")
    bx.set_ylabel("Named records")
    bx.xaxis.set_label_position("top")
    bx.xaxis.tick_top()
    bx.tick_params(length=0)
    for sp in bx.spines.values():
        sp.set_visible(False)
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
    """The two halves of COX1 from the same records; full trees are Figures S1 and S2."""
    fig = plt.figure(figsize=(6.9, 6.4))
    gs = fig.add_gridspec(1, 2, wspace=0.12)
    for k, (d, half) in enumerate((("results/phylogeny/paired/5prime", "5′"), ("results/phylogeny/paired/3prime", "3′"))):
        tips = read(os.path.join(d, "tips.csv"))
        n_rec = sum(int(t["records"]) for t in tips)
        ax = fig.add_subplot(gs[0, k])
        draw_tree(ax, os.path.join(d, "iqtree.treefile"), os.path.join(d, "tips.csv"),
                  f"{'ab'[k]}  {half} half, the same {n_rec} records ({len(tips)} haplotypes)", lw=0.6)
    save(fig, "Figure4")


# Figures S1-S4: the trees of Figure 4 at full size, with labelled tips ------------
SUPP_TREES = [
    ("results/phylogeny/5prime", "FigureS1",
     "Figure S1. Maximum-likelihood tree of the 5′ half of COX1 (reference positions 43–696), all records "
     "covering at least 80% of the window, collapsed to haplotypes (Figure 4a)."),
    ("results/phylogeny/3prime", "FigureS2",
     "Figure S2. Maximum-likelihood tree of the 3′ half of COX1 (reference positions 719–1,130), all records "
     "covering at least 80% of the window, collapsed to haplotypes (Figure 4b)."),
    ("results/phylogeny/paired/5prime", "FigureS3",
     "Figure S3. Maximum-likelihood tree of the 5′ half of COX1 for the 98 records that cover both windows "
     "(Figure 4c)."),
    ("results/phylogeny/paired/3prime", "FigureS4",
     "Figure S4. Maximum-likelihood tree of the 3′ half of COX1 for the 98 records that cover both windows "
     "(Figure 4d)."),
]


def ordered_tree(treefile):
    tree = Phylo.read(treefile, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    return tree


def node_support(c):
    sup = c.confidence
    if sup is None and c.name:
        try:
            sup = float(c.name.split("/")[-1])
        except ValueError:
            sup = None
    return sup


def tip_label(row):
    def names(field):
        parts = []
        for part in (row.get(field) or "").split("; "):
            m = re.match(r"(.+) \((\d+)\)$", part)
            if m:
                n = m.group(1)
                parts.append(f"{italic_species(n) if n.startswith('Tetranychus') else n.strip('()')} ({m.group(2)})")
        return ", ".join(parts)
    text = f"{row['haplotype']}   {names('species_names')}"
    if row.get("anchor_names"):
        text += f"   [anchor: {names('anchor_names')}]"
    return text


def supplementary_tree(d, name, caption):
    tree = ordered_tree(os.path.join(d, "iqtree.treefile"))
    rows = {r["haplotype"]: r for r in read(os.path.join(d, "tips.csv"))}
    tips = tree.get_terminals()
    n = len(tips)
    height = 0.115 * n + 2.5
    fig = plt.figure(figsize=(8.27, height))
    ax = fig.add_axes([0.03, 0.7 / height, 0.95, (height - 2.4) / height])
    ypos = {t: i for i, t in enumerate(tips)}
    depths = tree.depths()

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
    ax.add_collection(LineCollection(segs, colors=INK2, linewidths=0.4))
    xmax = max(depths.values())
    span = xmax / 0.42  # the tree takes 42% of the width
    for c in tree.get_nonterminals():
        sup = node_support(c)
        if sup is not None and sup >= 70 and c is not tree.root:
            ax.text(depths[c] - span * 0.002, ypos[c] - 0.15, f"{sup:.0f}", ha="right", va="bottom",
                    fontsize=3.6, color=INK2)
    x0, dx = xmax + span * 0.012, span * 0.011
    for t in tips:
        sp, anchor = tip_names_row(rows.get(t.name, {}))
        flags = [sp.get(COLUMNS[0][0]), sp.get(COLUMNS[1][0]), sp.get(COLUMNS[2][0]), sp.get(COLUMNS[3][0]),
                 any(k not in ("(unnamed)",) + tuple(nm for nm, _ in COLUMNS[:4]) for k in sp),
                 "(unnamed)" in sp, anchor]
        for j, on in enumerate(flags):
            if on:
                ax.add_patch(plt.Rectangle((x0 + j * dx, ypos[t] - 0.42), dx * 0.8, 0.84, color=COLUMNS[j][1], lw=0))
        ax.plot([depths[t], x0 - span * 0.004], [ypos[t]] * 2, color=GRID, lw=0.3, zorder=0)
        ax.text(x0 + len(COLUMNS) * dx + span * 0.008, ypos[t], tip_label(rows.get(t.name, {"haplotype": t.name})),
                va="center", ha="left", fontsize=5.2, color=INK)
    for j, (nm, _) in enumerate(COLUMNS):
        lab = italic_species(nm) if nm.startswith("Tetranychus") else nm
        ax.text(x0 + j * dx + dx * 0.4, -1.2, lab, rotation=90, ha="center", va="bottom", fontsize=5.2, color=INK)
    ax.set_xlim(-span * 0.005, span)
    ax.set_ylim(n - 0.5, -0.5)
    ax.axis("off")
    ax.plot([0, 0.05], [n + 1.2, n + 1.2], color=INK2, lw=0.6, clip_on=False)
    ax.text(0.055, n + 1.2, "0.05 substitutions/site", ha="left", va="center", fontsize=6, color=INK2, clip_on=False)
    fig.text(0.03, 1 - 0.25 / height, caption, ha="left", va="top", fontsize=7.5, color=INK, wrap=True)
    fig.text(0.03, 1 - 0.55 / height, "Tips are haplotypes, labelled with their identifier, the accepted names of the "
             "records carrying them (number of records) and any anchor specimens; Table S7 lists the accession of "
             "every record. Numbers at nodes are ultrafast bootstrap values of at least 70. Columns as in Figure 4.",
             ha="left", va="top", fontsize=6.2, color=INK2, wrap=True)
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, f"{name}.pdf"))
    plt.close(fig)


def tip_names_row(row):
    names = {}
    for part in (row.get("species_names") or "").split("; "):
        m = re.match(r"(.+) \((\d+)\)$", part)
        if m:
            names[m.group(1)] = int(m.group(2))
    return names, bool(row.get("anchor_names"))


def labelled_newick(d):
    """IQ-TREE tree with tips renamed to haplotype|names|records (Newick-safe)."""
    tree = Phylo.read(os.path.join(d, "iqtree.treefile"), "newick")
    rows = {r["haplotype"]: r for r in read(os.path.join(d, "tips.csv"))}
    for t in tree.get_terminals():
        r = rows.get(t.name)
        if r:
            names = re.sub(r"[^A-Za-z0-9_.+-]", "", r["species_names"].replace("Tetranychus ", "T.")
                           .replace(" (", "_").replace(")", "").replace("; ", "+").replace(" ", "_"))
            t.name = f"{r['haplotype']}|{names}|n{r['records']}"
    Phylo.write(tree, os.path.join(d, "labelled_tree.nwk"), "newick")


def supplementary_figures():
    for d, name, caption in SUPP_TREES:
        supplementary_tree(d, name, caption)
        labelled_newick(d)


# Figure S5 -----------------------------------------------------------------
def figure_s5():
    rows = read("results/trends/by_period.csv")
    periods = [r["period"].replace("before 2006", "Before 2006").replace("2021 onwards", "2021 to\n29 Oct 2026")
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
    pub = series[3][1][-1]
    ax.annotate(f"{pub:.1f}%", (x[-1], pub), xytext=(-6, 9), textcoords="offset points", ha="right",
                fontsize=6.5, color=INK2)
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
    save(fig, "FigureS5")


# Figure 5 ------------------------------------------------------------------
O_STUDY_LABEL = {"matsuda2012": "Matsuda et al. 2012", "arabuli2018": "Arabuli & Gotoh 2018",
                 "mushtaq2023a": "Mushtaq et al. 2023a", "mushtaq2023b": "Mushtaq et al. 2023b",
                 "majeed2026": "Safeena Majeed et al. 2026", "guzman2014": "Guzmán-Valencia et al. 2014",
                 "guzman2017": "Guzmán-Valencia et al. 2017", "sakamoto2019": "Sakamoto et al. 2019",
                 "majeed2022": "Safeena Majeed et al. 2022", "muralidharan2020": "Muralidharan et al. 2020"}


def stacked_rows(ax, rows, key="all_result"):
    """One horizontal stacked bar per (label, queries) row."""
    for y, (label, rs) in enumerate(rows):
        outcome_bar(ax, y, rs, key, height=0.64)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label for label, _ in rows], fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Anchor queries (%)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def figure5():
    """Oligonychus beside Tetranychus, all-records reference set (what a database user meets)."""
    tq = read("results/robustness/queries_main.csv")
    oq = read("results/oligonychus/robustness/queries_main.csv")
    grouped = read("results/oligonychus/robustness/queries_grouped.csv")

    def lab(genus, what, rs):
        return f"$\\mathit{{{genus}}}$, {what}: {len(rs)} queries, {len({r['study_id'] for r in rs})} studies"
    tA = [r for r in tq if r["evidence_level"] == "A"]
    oA = [r for r in oq if r["evidence_level"] == "A"]
    summary = [(lab("Tetranychus", "all studies", tq), tq), (lab("Tetranychus", "level A", tA), tA),
               (lab("Oligonychus", "all studies", oq), oq), (lab("Oligonychus", "level A", oA), oA),
               (lab("Oligonychus", "Mexican studies held out together", grouped), grouped)]
    by = defaultdict(list)
    for r in oq:
        by[r["study_id"]].append(r)
    order = sorted(by, key=lambda sid: (min(r["evidence_level"] for r in by[sid]), -len(by[sid])))
    studies = [(f"{O_STUDY_LABEL[sid]} [{'/'.join(sorted({r['evidence_level'] for r in by[sid]}))}]: "
                f"{len(by[sid])} quer{'y' if len(by[sid]) == 1 else 'ies'}", by[sid]) for sid in order]
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(5.8, 5.2), gridspec_kw={"height_ratios": [5, 10]})
    stacked_rows(ax, summary)
    stacked_rows(bx, studies)
    outcome_legend(ax, 1.02)
    for a_, letter in ((ax, "a"), (bx, "b")):
        a_.text(-0.02, 1.02 if a_ is bx else 1.14, letter, transform=a_.transAxes, fontsize=10,
                fontweight="bold", va="bottom", ha="right")
    fig.tight_layout(h_pad=1.6)
    save(fig, "Figure5")


def main():
    figure1()
    figure2()
    figure3()
    figure4()
    figure_s5()
    figure5()
    supplementary_figures()
    print(f"Figures written to {OUT}")


if __name__ == "__main__":
    main()
