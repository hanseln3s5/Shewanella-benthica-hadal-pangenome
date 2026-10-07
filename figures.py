"""Publication figures for the S. benthica pangenome manuscript.
Usage: python3 figures.py <resdir> <figdir> <input.xlsx>
"""
import sys, json, itertools
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.lines import Line2D

res, out = sys.argv[1], sys.argv[2]
import os; os.makedirs(out, exist_ok=True)
R = json.load(open(f"{res}/results.json"))

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.edgecolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.labelcolor": "#0b0b0b",
    "text.color": "#0b0b0b", "savefig.dpi": 400, "pdf.fonttype": 42, "svg.fonttype": "none",
})
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
UNK = "#c9c8c2"
COMP_COL = {"Core": CAT[0], "Accessory": CAT[1], "Strain-specific": CAT[2]}
STRAINS = ["BOEU19-1", "DSM8812", "DB21MT-2", "KT99"]  # non-hadal first, hadal last
HADAL = {"DB21MT-2", "KT99"}
CATS = ["Core", "Accessory", "Strain-specific"]

def panel(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    fig.savefig(f"{out}/{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{out}/{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

# =====================================================================
# Figure 1 - workflow
# =====================================================================
fig, ax = plt.subplots(figsize=(7.2, 3.0)); ax.set_xlim(-0.5, 101); ax.set_ylim(0, 40); ax.axis("off")
steps = [
    ("1  Genomes", ["4 S. benthica genomes", "ncbi-genome-download", "BOEU19-1, DSM8812", "DB21MT-2*, KT99*"]),
    ("2  Annotation", ["Prokka", "uniform gene calling", "and annotation", "GFF3 + proteins"]),
    ("3  Pangenome", ["Panaroo", "graph-based clustering", "core 4/4, accessory 2-3/4", "strain-specific 1/4"]),
    ("4  Function", ["eggNOG-mapper", "COG, KEGG KO, Pfam", "curated functional", "classes"]),
    ("5  Analysis", ["intersections, Jaccard", "rarefaction, Heaps' law", "Fisher tests, BH-FDR", "trait inventory, anvi'o"]),
]
w, h, gap, y0 = 18.4, 30, 2.0, 8
for i, (title, lines) in enumerate(steps):
    x = 0.5 + i * (w + gap)
    ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.25,rounding_size=1.6", fc="#f4f6fa", ec=CAT[0], lw=0.9))
    ax.add_patch(FancyBboxPatch((x, y0 + h - 7), w, 7, boxstyle="round,pad=0.25,rounding_size=1.6", fc=CAT[0], ec=CAT[0], lw=0.9))
    ax.text(x + w / 2, y0 + h - 3.5, title, ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
    for j, l in enumerate(lines):
        st = "italic" if (i == 0 and j == 0) else "normal"
        ax.text(x + w / 2, y0 + h - 11.5 - j * 4.9, l, ha="center", va="center", fontsize=6.3, color=INK, style=st)
    if i < len(steps) - 1:
        ax.annotate("", xy=(x + w + gap - 0.2, y0 + h / 2), xytext=(x + w + 0.4, y0 + h / 2),
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.0, mutation_scale=9))
ax.text(50, 3.2, "* hadal isolates (> 6,000 m).   Outputs: Figures 2-5 and Supplementary Tables S1-S10",
        ha="center", fontsize=7, color=INK2)
save(fig, "Figure1_workflow")

# =====================================================================
# Figure 2 - pangenome architecture
# =====================================================================
sub = pd.read_csv(f"{res}/subset_rarefaction.csv")
X = pd.read_excel(sys.argv[3], sheet_name=None)
it = X["Intersections"].copy()
it["set"] = it.Combination.str.split(r" \+ ").map(frozenset)
it = it.sort_values("Families", ascending=False).reset_index(drop=True)
fam_per = R["families_per_strain"]

fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1.25, 1.05], hspace=0.62, wspace=0.55)
# a compartments
ax = fig.add_subplot(gs[0, 0])
vals = [1893, 2043, 1858]
b = ax.bar(CATS, vals, color=[COMP_COL[c] for c in CATS], width=0.66, edgecolor="white", linewidth=1)
for r_, v in zip(b, vals):
    ax.text(r_.get_x() + r_.get_width() / 2, v + 40, f"{v:,}\n({100*v/5794:.1f}%)", ha="center", va="bottom", fontsize=6.8, color=INK)
ax.set_ylim(0, 2700); ax.set_ylabel("Gene families"); ax.set_xticks(range(3)); ax.set_xticklabels(["Core", "Accessory", "Strain-\nspecific"])
ax.set_title("Compartments (n = 5,794)", loc="left", fontsize=8.5)
panel(ax, "a", x=-0.32)
# b frequency spectrum
ax = fig.add_subplot(gs[0, 1])
fv = [1858, 1015, 1028, 1893]; cols = [COMP_COL["Strain-specific"], COMP_COL["Accessory"], COMP_COL["Accessory"], COMP_COL["Core"]]
ax.bar([1, 2, 3, 4], fv, color=cols, width=0.66, edgecolor="white", linewidth=1)
for x, v in zip([1, 2, 3, 4], fv): ax.text(x, v + 40, f"{v:,}", ha="center", va="bottom", fontsize=6.8)
ax.set_xticks([1, 2, 3, 4]); ax.set_xlabel("Genomes containing the family"); ax.set_ylabel("Gene families"); ax.set_ylim(0, 2300)
ax.set_title("Frequency spectrum", loc="left", fontsize=8.5)
panel(ax, "b", x=-0.32)
# c Jaccard
ax = fig.add_subplot(gs[0, 2])
J = np.ones((4, 4))
for k, v in R["pairwise"].items():
    a_, b_ = k.split("|"); i, j = STRAINS.index(a_), STRAINS.index(b_); J[i, j] = J[j, i] = v["jaccard"]
blues = matplotlib.colors.LinearSegmentedColormap.from_list("b", ["#f0efec", "#9ec5f4", "#2a78d6", "#104281"])
Jm = np.ma.masked_where(np.eye(4) == 1, J)
im = ax.imshow(Jm, cmap=blues, vmin=0.40, vmax=0.70)
for i in range(4):
    for j in range(4):
        if i != j:
            ax.text(j, i, f"{J[i,j]:.2f}", ha="center", va="center", fontsize=6.8, color="white" if J[i, j] > 0.6 else INK)
        else:
            ax.text(j, i, f"{fam_per[STRAINS[i]]:,}", ha="center", va="center", fontsize=6, color=INK2)
ax.set_xticks(range(4)); ax.set_yticks(range(4)); ax.set_xticklabels(STRAINS, rotation=40, ha="right"); ax.set_yticklabels(STRAINS)
for t in ax.get_xticklabels() + ax.get_yticklabels():
    if t.get_text() in HADAL: t.set_fontweight("bold")
ax.spines[:].set_visible(False); ax.tick_params(length=0)
ax.set_title("Jaccard similarity", loc="left", fontsize=8.5)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)
panel(ax, "c", x=-0.42)
# d UpSet
gsd = gs[1, :].subgridspec(2, 2, width_ratios=[1, 4.6], height_ratios=[1.5, 1], hspace=0.05, wspace=0.03)
axb = fig.add_subplot(gsd[0, 1]); axm = fig.add_subplot(gsd[1, 1], sharex=axb)
x = np.arange(len(it))
def colfor(s):
    n = len(s); return COMP_COL["Core"] if n == 4 else (COMP_COL["Strain-specific"] if n == 1 else COMP_COL["Accessory"])
axb.bar(x, it.Families, color=[colfor(s) for s in it.set], width=0.7, edgecolor="white", linewidth=0.8)
for xi, v in zip(x, it.Families): axb.text(xi, v + 25, f"{v:,}", ha="center", va="bottom", fontsize=6.2)
axb.set_ylabel("Families in\nintersection"); axb.tick_params(axis="x", bottom=False, labelbottom=False); axb.set_ylim(0, 2150)
axb.set_title("Pangenome intersections (exclusive sets)", loc="left", fontsize=8.5)
panel(axb, "d", x=-0.3)
ys = {s: i for i, s in enumerate(STRAINS[::-1])}
for xi, s in zip(x, it.set):
    axm.scatter([xi] * 4, range(4), s=22, color="#e4e3df", zorder=1)
    yy = sorted(ys[k] for k in s)
    axm.plot([xi, xi], [yy[0], yy[-1]], color=INK, lw=1.2, zorder=2)
    axm.scatter([xi] * len(yy), yy, s=26, color=INK, zorder=3)
axm.set_yticks(range(4)); axm.set_yticklabels([f"{s} ({fam_per[s]:,})" for s in STRAINS[::-1]]); axm.tick_params(axis="x", bottom=False, labelbottom=False)
for t in axm.get_yticklabels():
    if t.get_text().split()[0] in HADAL: t.set_fontweight("bold")
axm.spines[["left", "bottom"]].set_visible(False); axm.tick_params(axis="y", length=0)
for i in range(4):
    if i % 2 == 0: axm.axhspan(i - 0.5, i + 0.5, color="#f6f5f2", zorder=0)
axm.set_ylim(-0.5, 3.5)
axb.legend(handles=[Patch(color=COMP_COL[c], label=c) for c in CATS], frameon=False, loc="upper right", ncol=3)
# e rarefaction
ax = fig.add_subplot(gs[2, 0:2])
H = R["heaps_panaroo"]
for k_ in ["pangenome", "core"]:
    col = CAT[0] if k_ == "pangenome" else CAT[1]
    jit = (np.random.RandomState(1).rand(len(sub)) - 0.5) * 0.12
    ax.scatter(sub.n_genomes + jit, sub[k_], s=14, color=col, alpha=0.55, edgecolor="none")
    mm = sub.groupby("n_genomes")[k_].mean()
    ax.plot(mm.index, mm.values, color=col, lw=2, marker="o", ms=4, label=("Pangenome" if k_ == "pangenome" else "Core genome") + " (mean of all subsets)")
NN = np.linspace(1, 4, 50)
ax.plot(NN, H["K"] * NN ** H["gamma"], ls="--", color=INK2, lw=1.2,
        label=f"Heaps' law fit: n = {H['K']:.0f}·N$^{{{H['gamma']:.2f}}}$ (γ = {H['gamma']:.2f} ± {H['gamma_se']:.2f}; R² = {H['R2']:.2f})")
ax.set_xticks([1, 2, 3, 4]); ax.set_xlabel("Number of genomes (N)"); ax.set_ylabel("Gene families"); ax.set_ylim(1500, 6300)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, fontsize=6.4); ax.grid(axis="y", color=GRID, lw=0.5)
ax.set_title("Accumulation curves (all 15 strain subsets)", loc="left", fontsize=8.5)
panel(ax, "e", x=-0.09)
# f new genes per added genome
ax = fig.add_subplot(gs[2, 2])
ng = R["new_genes_mean"]; kk = [2, 3, 4]; vv = [ng[str(k)] for k in kk]
ax.bar(kk, vv, color=CAT[2], width=0.6)
for k, v in zip(kk, vv): ax.text(k, v + 20, f"{v:.0f}", ha="center", va="bottom", fontsize=6.8)
ax.set_xticks(kk); ax.set_xlabel("k-th genome added"); ax.set_ylabel("New families (mean)"); ax.set_ylim(0, 1250)
ax.set_title("New families per genome", loc="left", fontsize=8.5)
panel(ax, "f", x=-0.32)
save(fig, "Figure2_pangenome_architecture")

# =====================================================================
# Figure 4 - strain-specific repertoires and mobilome
# =====================================================================
cs = pd.read_csv(f"{res}/class_by_strain.csv", index_col=0)[STRAINS]
isf = pd.read_csv(f"{res}/IS_by_strain.csv", index_col=0)[STRAINS]
subs = pd.read_csv(f"{res}/subclass_by_strain.csv", index_col=0).reindex(columns=STRAINS).fillna(0)
CLS = list(cs.index)
CLS_COL = dict(zip(CLS, CAT[:8] + [UNK]))
short = {CLS[0]: "Mobilome: IS elements / transposases", CLS[1]: "Mobilome: prophages, integrases, plasmids",
         CLS[2]: "Defence (CRISPR-Cas, R-M, TA)", CLS[3]: "Cell envelope biogenesis (M)",
         CLS[4]: "Signal transduction & transcription (T, K)", CLS[5]: "Ion transport & secretion (P, U)",
         CLS[6]: "Metabolism (C, E, F, G, H, I, Q)", CLS[7]: "Other characterized functions", CLS[8]: "Unknown / poorly characterized (S, none)"}
fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], width_ratios=[1.15, 1], hspace=0.72, wspace=0.75)
ax = fig.add_subplot(gs[0, :])
left = np.zeros(4); yy = np.arange(4)[::-1]
for c in CLS:
    v = cs.loc[c].values
    ax.barh(yy, v, left=left, color=CLS_COL[c], height=0.62, edgecolor="white", linewidth=1.2, label=short[c])
    for y_, l_, v_ in zip(yy, left, v):
        if v_ >= 30: ax.text(l_ + v_ / 2, y_, f"{int(v_)}", ha="center", va="center", fontsize=6.4, color="white" if c != CLS[8] else INK)
    left += v
for y_, t in zip(yy, left): ax.text(t + 8, y_, f"{int(t)}", va="center", fontsize=7.5, fontweight="bold")
ax.set_yticks(yy); ax.set_yticklabels([f"{s}\n({'hadal' if s in HADAL else 'non-hadal'})" for s in STRAINS])
ax.set_xlabel("Strain-specific gene families"); ax.set_xlim(0, 760)
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=6.6, handlelength=1.2, columnspacing=1.2)
ax.set_title("Functional composition of strain-specific gene families", loc="left")
panel(ax, "a", x=-0.1)
# b IS families
ax = fig.add_subplot(gs[1, 0])
isf = isf.loc[isf.sum(axis=1).sort_values(ascending=False).index]
IS_COL = dict(zip(isf.index, (CAT[:8] + [UNK])[:len(isf)]))
if "Unassigned transposase" in IS_COL: IS_COL["Unassigned transposase"] = UNK
bottom = np.zeros(4)
for f in isf.index:
    v = isf.loc[f].values
    ax.bar(range(4), v, bottom=bottom, color=IS_COL[f], width=0.62, edgecolor="white", linewidth=1, label=f.replace(" transposase", ""))
    bottom += v
for i, t in enumerate(bottom): ax.text(i, t + 3, f"{int(t)}", ha="center", va="bottom", fontsize=7.5, fontweight="bold")
ax.set_xticks(range(4)); ax.set_xticklabels(STRAINS); ax.set_ylabel("Strain-specific transposase families"); ax.set_ylim(0, 175); ax.set_xlim(-0.5, 3.5)
for t in ax.get_xticklabels():
    if t.get_text() in HADAL: t.set_fontweight("bold")
ax.legend(frameon=False, fontsize=6.0, loc="upper left", title="IS family (Pfam/eggNOG)", title_fontsize=6.3, handlelength=1.0)
ax.set_title("Transposase families by IS family", loc="left")
panel(ax, "b", x=-0.2)
# c defence & mobile subclasses dot matrix
ax = fig.add_subplot(gs[1, 1])
order = ["Transposase / IS element", "Integrase / recombinase", "Phage structural / lysis", "CRISPR-Cas", "Restriction-modification", "Toxin-antitoxin / abortive infection"]
subs = subs.reindex(order).fillna(0)
for i, r in enumerate(order):
    for j, s in enumerate(STRAINS):
        v = subs.loc[r, s]
        if v > 0:
            ax.scatter(j, i, s=10 + 14 * np.sqrt(v), color=CAT[0] if s in HADAL else MUTED, alpha=0.9, edgecolor="white", linewidth=0.6)
            ax.text(j + 0.22, i, f"{int(v)}", va="center", fontsize=6.3)
        else:
            ax.scatter(j, i, s=6, color="#e4e3df")
ax.set_xticks(range(4)); ax.set_xticklabels(STRAINS, rotation=30, ha="right"); ax.set_yticks(range(len(order)))
ax.set_yticklabels(["Transposases", "Integrases", "Phage structural", "CRISPR-Cas (I-F)", "R-M systems", "TA / Abi"])
for t in ax.get_xticklabels():
    if t.get_text() in HADAL: t.set_fontweight("bold")
ax.invert_yaxis(); ax.set_xlim(-0.5, 3.8); ax.spines[["left", "bottom"]].set_visible(False); ax.tick_params(length=0)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=CAT[0], label="hadal strain"), Line2D([], [], marker="o", ls="", color=MUTED, label="non-hadal strain")],
          frameon=False, fontsize=6.3, loc="lower center", bbox_to_anchor=(0.45, -0.42), ncol=2)
ax.set_title("Mobilome and defence families", loc="left")
panel(ax, "c", x=-0.75)
save(fig, "Figure4_strain_specific_mobilome")

# =====================================================================
# Figure 3 - COG composition & enrichment
# =====================================================================
cog = pd.read_csv(f"{res}/cog_enrichment.csv")
cog = cog[(cog[[f"{c}_n" for c in CATS]].sum(axis=1) >= 5)].reset_index(drop=True)
COGN = {"A": "RNA processing", "B": "Chromatin", "C": "Energy production", "D": "Cell cycle control", "E": "Amino acid metabolism",
        "F": "Nucleotide metabolism", "G": "Carbohydrate metabolism", "H": "Coenzyme metabolism", "I": "Lipid metabolism",
        "J": "Translation & ribosome", "K": "Transcription", "L": "Replication, recombination & repair", "M": "Cell wall/membrane/envelope",
        "N": "Cell motility", "O": "PTM, turnover, chaperones", "P": "Inorganic ion transport", "Q": "Secondary metabolites",
        "S": "Function unknown", "T": "Signal transduction", "U": "Secretion & vesicular transport", "V": "Defence mechanisms", "W": "Extracellular structures", "Z": "Cytoskeleton"}
fig = plt.figure(figsize=(7.2, 5.2))
gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], height_ratios=[1, 0.025], hspace=0.16, wspace=0.08)
ax = fig.add_subplot(gs[0, 0])
M = cog[[f"{c}_pct" for c in CATS]].values
seq = matplotlib.colors.LinearSegmentedColormap.from_list("s", ["#f7f9fc", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
im = ax.imshow(M, cmap=seq, aspect="auto", vmin=0, vmax=24)
for i in range(M.shape[0]):
    for j in range(3):
        ax.text(j, i, f"{M[i,j]:.1f}", ha="center", va="center", fontsize=6.2, color="white" if M[i, j] > 11 else INK)
ax.set_yticks(range(len(cog))); ax.set_yticklabels([f"{c}  {COGN.get(c,'')}" for c in cog.COG], fontsize=6.8)
ax.set_xticks(range(3)); ax.set_xticklabels(["Core", "Accessory", "Strain-\nspecific"], fontsize=7)
ax.xaxis.tick_top(); ax.spines[:].set_visible(False); ax.tick_params(length=0)
cax = fig.add_subplot(gs[1, 0]); cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label("% of COG-annotated families", fontsize=6.5); cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)
ax.text(-0.02, 1.075, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right")
for k, (cmp, title, letter) in enumerate([("Accessory", "Accessory vs core", "b"), ("Strain-specific", "Strain-specific vs core", "c")]):
    ax2 = fig.add_subplot(gs[0, k + 1], sharey=ax)
    fc = cog[f"{cmp}_vs_Core_log2FC"].values; q = cog[f"{cmp}_vs_Core_q"].values
    colors = [("#e34948" if f > 0 else "#2a78d6") if qq < 0.05 else "#d6d5cf" for f, qq in zip(fc, q)]
    ax2.barh(range(len(cog)), fc, color=colors, height=0.66)
    for i, (f, qq) in enumerate(zip(fc, q)):
        if qq < 0.05:
            st = "***" if qq < 0.001 else ("**" if qq < 0.01 else "*")
            ax2.text(f + (0.08 if f > 0 else -0.08), i, st, va="center", ha="left" if f > 0 else "right", fontsize=6.5)
    ax2.axvline(0, color=INK2, lw=0.6); ax2.set_xlim(-3.3, 3.3); ax2.tick_params(axis="y", left=False, labelleft=False)
    ax2.spines["left"].set_visible(False); ax2.set_xlabel("log$_2$ fold change in proportion"); ax2.grid(axis="x", color=GRID, lw=0.5)
    ax2.set_title(title, fontsize=8.5)
    ax2.text(0.0, 1.075, letter, transform=ax2.transAxes, fontsize=11, fontweight="bold")
fig.text(0.62, 0.01, "Fisher's exact test, Benjamini-Hochberg q: * < 0.05, ** < 0.01, *** < 0.001; grey bars, q ≥ 0.05", fontsize=6.5, color=INK2, ha="center")
save(fig, "Figure3_COG_enrichment")

# =====================================================================
# Figure 5 - trait inventory dot matrix
# =====================================================================
tr = pd.read_csv(f"{res}/traits.csv")
cols = ["Core", "Accessory"] + [f"SS_{s}" for s in STRAINS]
labs = ["Core\n(4/4)", "Accessory\n(2-3/4)"] + [f"{s}\n(specific)" for s in STRAINS]
fig, ax = plt.subplots(figsize=(7.2, 7.0))
modules = list(dict.fromkeys(tr.Module))
ypos, y, yt, ytl, mod_y = [], 0, [], [], []
for mname in modules:
    sub_ = tr[tr.Module == mname]
    mod_y.append((mname, y))
    y += 0.9
    for _, r in sub_.iterrows():
        ypos.append((r, y)); yt.append(y); ytl.append(r.Trait); y += 1
    y += 0.35
for r, yy_ in ypos:
    for j, c in enumerate(cols):
        v = r[c]
        if v > 0:
            colr = COMP_COL["Core"] if j == 0 else (COMP_COL["Accessory"] if j == 1 else (CAT[6] if c[3:] in HADAL else MUTED))
            ax.scatter(j, yy_, s=10 + 18 * np.sqrt(v), color=colr, alpha=0.9, edgecolor="white", linewidth=0.5, zorder=3)
            ax.text(j + 0.2, yy_, str(int(v)), va="center", fontsize=6, color=INK2)
        else:
            ax.scatter(j, yy_, s=5, color="#e4e3df", zorder=2)
for mname, yy_ in mod_y:
    ax.text(-3.55, yy_ + 0.25, mname, fontsize=7.6, fontweight="bold", va="center", ha="left", transform=ax.transData, clip_on=False)
ax.set_yticks(yt); ax.set_yticklabels(ytl, fontsize=6.9); ax.invert_yaxis()
ax.set_xticks(range(len(cols))); ax.set_xticklabels(labs, fontsize=7); ax.xaxis.tick_top()
for t in ax.get_xticklabels():
    if any(h in t.get_text() for h in HADAL): t.set_fontweight("bold")
ax.axvline(1.5, color=GRID, lw=0.8); ax.set_xlim(-0.5, len(cols) - 0.35)
ax.spines[:].set_visible(False); ax.tick_params(length=0)
for sz in [1, 10, 50, 150]:
    ax.scatter([], [], s=10 + 18 * np.sqrt(sz), color=MUTED, label=f"{sz}")
ax.legend(title="Gene families", frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.07), ncol=4, fontsize=6.5, title_fontsize=6.8)
save(fig, "Figure5_trait_inventory")
print("ok")
