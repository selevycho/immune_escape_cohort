#!/usr/bin/env python3
"""
Line 1 convergence figure, recomputed on the final patient sets.

Two panels answer two different objections.

Panel A shows that the two cancers rank genes the same way. If breast and
ovarian tumours silenced arbitrary genes, the scatter would be a cloud;
instead it follows the diagonal, and the genes at each extreme are the
same in both.

Panel B answers "what if you had picked a different RSEM cutoff?". The
silenced fraction rises with any threshold - that is arithmetic, not a
finding. What matters is that the gap between the two cancers holds at
every cutoff tested, so the choice of 5 is not doing the work.

Usage:
  python plot_line1_convergence_v2.py <tcga_root> <out_dir>
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = sys.argv[1]
OUT = sys.argv[2]
RSEM_OFF = 5.0
THRESHOLDS = [1, 5, 10, 50]
MIN_MUT = {"brca": 10, "ov": 5}

NONSYN = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
          "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
          "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}

BLUE = "#1F6FA8"
ORANGE = "#D98324"
GREY = "#B8B6AE"

os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#666666", "axes.labelcolor": "#222222",
    "xtick.color": "#444444", "ytick.color": "#444444",
})


def norm(b):
    p = str(b).split("-")
    return "-".join(p[:4])[:15] if len(p) >= 4 else str(b)


data = {}

for cohort in ["brca", "ov"]:
    print(f"loading {cohort.upper()} ...", flush=True)
    m = pd.read_csv(f"{ROOT}/{cohort}/data_mutations.txt", sep="\t",
                    comment="#", low_memory=False)
    e = pd.read_csv(f"{ROOT}/{cohort}/data_mrna_seq_v2_rsem.txt",
                    sep="\t", low_memory=False)
    c = pd.read_csv(f"{ROOT}/{cohort}/data_cna.txt", sep="\t", low_memory=False)

    e = e[e.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol")
    e = e.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"], errors="ignore")
    e.columns = [norm(x) for x in e.columns]
    e = e.loc[:, ~e.columns.duplicated()].clip(lower=0)

    cna_samples = {norm(x) for x in c.columns if str(x).startswith("TCGA")}
    usable = set(e.columns) & cna_samples

    m["sample"] = m.Tumor_Sample_Barcode.map(norm)
    m = m[m.Variant_Classification.isin(NONSYN)]
    m = m[m["sample"].isin(usable)]
    m = m[m.Hugo_Symbol.isin(e.index)]

    idx = pd.MultiIndex.from_arrays([m.Hugo_Symbol.values, m["sample"].values])
    m = m.assign(rsem=e.stack().reindex(idx).values)
    m = m[m.rsem.notna()].copy()

    print(f"  {m['sample'].nunique()} patients, {len(m):,} mutations", flush=True)

    m["silenced"] = m.rsem < RSEM_OFF
    g = m.groupby("Hugo_Symbol").agg(
        n=("silenced", "size"), s=("silenced", "sum")).reset_index()
    g["pct"] = 100 * g.s / g.n
    g = g[g.n >= MIN_MUT[cohort]]

    curve = [100 * (m.rsem < t).mean() for t in THRESHOLDS]

    data[cohort] = {"genes": g.set_index("Hugo_Symbol"),
                    "curve": curve,
                    "patients": m["sample"].nunique(),
                    "mutations": len(m)}

# ---------------- shared genes ----------------
gb, go = data["brca"]["genes"], data["ov"]["genes"]
shared = gb.index.intersection(go.index)
x = gb.loc[shared, "pct"].values
y = go.loc[shared, "pct"].values
rho, p = stats.spearmanr(x, y)
print(f"\nshared genes {len(shared)}, Spearman rho = {rho:.3f}, p = {p:.3g}")

df = pd.DataFrame({"gene": shared, "brca": x, "ov": y})
# Label the genes that carry the most mutations, not the ones nearest the
# cutoff: nlargest on the silencing percentage picks borderline genes and
# leaves TP53 and TTN, the two the audience will look for, unlabelled.
low = df[(df.brca < 5) & (df.ov < 5)].copy()
low["weight"] = gb.loc[low.gene, "n"].values + go.loc[low.gene, "n"].values
low = low.nlargest(8, "weight")
high = df[(df.brca > 95) & (df.ov > 95)].copy()
high["weight"] = gb.loc[high.gene, "n"].values + go.loc[high.gene, "n"].values
high = high.nlargest(8, "weight")
low_names = set(low.gene)
high_names = set(high.gene)
print(f"  shared LOW  : {', '.join(sorted(low_names))}")
print(f"  shared HIGH : {', '.join(sorted(high_names))}")

# =====================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.2))
fig.suptitle("Line 1  |  Transcriptomic silencing converges between cancers",
             fontsize=15, fontweight="bold", x=0.02, ha="left", y=0.98)

# ---------------- panel A ----------------
# Both highlighted groups sit at exactly 0/0 and 100/100, so eight points
# collapse into one mark. A small jitter separates them without moving any
# point far enough to misrepresent its value.
rng = np.random.default_rng(7)


def jitter(v, scale=1.6):
    return v + rng.normal(0, scale, len(v))


# Seven of the eight genes in each group sit at exactly 0/0 or 100/100,
# so the marks collapse into one point. A small jitter separates them.
rng = np.random.default_rng(7)
def jit(v, s=1.7):
    return np.asarray(v) + rng.normal(0, s, len(v))

other = df[~df.gene.isin(low_names | high_names)]
ax1.scatter(jitter(other.brca.values), jitter(other.ov.values),
            s=14, c=GREY, alpha=0.5, edgecolors="none",
            label=f"other genes (n={len(other)})", zorder=2)
ax1.scatter(jitter(low.brca.values), jitter(low.ov.values),
            s=70, c=BLUE, edgecolors="white", linewidths=0.8,
            label=f"top {len(low)} kept expressed", zorder=4)
ax1.scatter(jitter(high.brca.values), jitter(high.ov.values),
            s=70, c=ORANGE, edgecolors="white", linewidths=0.8,
            label=f"top {len(high)} silenced", zorder=4)
ax1.plot([0, 100], [0, 100], "--", c=GREY, lw=1, zorder=1)

# Gene names are omitted from the figure. Both groups sit in tight clusters
# at the two corners, so eight labels each either overlap or trail far from
# their points; the names carry better as text on the slide beside it.

ax1.set_xlabel("Silenced mutations in BREAST (%)", fontsize=12)
ax1.set_ylabel("Silenced mutations in OVARIAN (%)", fontsize=12)
ax1.set_title("A   The same genes, ranked the same way",
              fontsize=13, fontweight="bold", loc="left", pad=12)
n_zero = int(((df.brca == 0) & (df.ov == 0)).sum())
n_full = int(((df.brca == 100) & (df.ov == 100)).sum())
ax1.text(0.03, 0.97,
         f"Spearman rho = {rho:.3f}\np = {p:.1e}\nn = {len(shared)} genes\n"
         f"{n_zero} at 0% in both, {n_full} at 100%",
         transform=ax1.transAxes, va="top", fontsize=11,
         bbox=dict(boxstyle="round,pad=0.5", fc="#F2F1EC", ec="none"))
ax1.legend(loc="lower right", frameon=False, fontsize=10)
ax1.set_xlim(-6, 112)
ax1.set_ylim(-6, 112)
ax1.grid(alpha=0.15, lw=0.5)

# ---------------- panel B ----------------
xs = range(len(THRESHOLDS))
ax2.plot(xs, data["brca"]["curve"], "-o", c=BLUE, lw=2.2, ms=9,
         label=f"BRCA  ({data['brca']['patients']} patients)")
ax2.plot(xs, data["ov"]["curve"], "-s", c=ORANGE, lw=2.2, ms=9,
         label=f"OV  ({data['ov']['patients']} patients)")

for i, (b, o) in enumerate(zip(data["brca"]["curve"], data["ov"]["curve"])):
    ax2.annotate(f"{b:.1f}%", (i, b), xytext=(0, 11),
                 textcoords="offset points", ha="center",
                 fontsize=10, color=BLUE)
    ax2.annotate(f"{o:.1f}%", (i, o), xytext=(0, -18),
                 textcoords="offset points", ha="center",
                 fontsize=10, color=ORANGE)

cut = THRESHOLDS.index(int(RSEM_OFF))
ax2.axvline(cut, ls=":", c=GREY, lw=1.2)
ax2.annotate("chosen\ncutoff", (cut, ax2.get_ylim()[0]),
             xytext=(8, 14), textcoords="offset points",
             fontsize=10, color="#777777")

ax2.set_xticks(list(xs))
ax2.set_xticklabels([str(t) for t in THRESHOLDS])
ax2.set_xlabel('RSEM threshold for "not expressed"', fontsize=12)
ax2.set_ylabel("Silenced fraction of all mutations (%)", fontsize=12)
ax2.set_title("B   The gap holds at every threshold",
              fontsize=13, fontweight="bold", loc="left", pad=12)
ax2.legend(frameon=False, fontsize=11, loc="upper left")
ax2.grid(alpha=0.15, lw=0.5, axis="y")

plt.tight_layout(rect=[0, 0, 1, 0.95])
out = f"{OUT}/line1_convergence_v2.png"
plt.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
print(f"\nwritten to {out}")

print(f"\n{'='*60}\n NUMBERS FOR THE SLIDE\n{'='*60}")
print(f"  shared genes          {len(shared)}")
print(f"  Spearman rho          {rho:.3f}")
print(f"  p                     {p:.2e}")
print(f"\n  kept expressed in both:  {', '.join(sorted(low_names))}")
print(f"  silenced in both:        {', '.join(sorted(high_names))}")
print(f"\n  threshold sensitivity:")
print(f"    {'RSEM <':<10}{'BRCA':>10}{'OV':>10}{'gap':>8}")
for t, b, o in zip(THRESHOLDS, data["brca"]["curve"], data["ov"]["curve"]):
    print(f"    {t:<10}{b:>9.1f}%{o:>9.1f}%{b-o:>7.1f}")
