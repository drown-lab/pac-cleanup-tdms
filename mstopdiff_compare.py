"""
Compare MSTopDiff delta-mass histograms across three sample-cleanup conditions
and annotate common top-down mass shifts.

MSTopDiff parameters used to make the inputs:
  Mass feature filter 0-100 kDa | RT 0-60 min | delta mass 0-150 Da
  RT window 2 min | max charge difference 2 | bin size 0.01 Da

Because each file holds only a few hundred features, the 0.01 Da histogram is
dominated by combinatorial background (~25-35 pairs/bin). So instead of
untargeted peak-picking we:
  (1) rebin to 0.1 Da for a readable overlay, and
  (2) for each candidate modification, sum the signal in a +/-0.05 Da window
      and compute ENRICHMENT over the local background (median of flanking
      bins). Enrichment is condition-internal, so it is directly comparable
      across conditions despite their different total pair counts.

We focus on the intensity x count signal (count is dominated by per-condition
baseline offsets and is uninformative here).

Outputs:
  fig_mstopdiff_ixc.{png,pdf}            - diverging intensity x count overlay
                                           (lower above / higher below axis),
                                           full 0-150 Da + annotated 0-60 zoom
  fig_mstopdiff_mod_enrichment.{png,pdf} - grouped bars of intensity x count
                                           enrichment (log x)
  mstopdiff_mod_table.csv                - per-mod, per-condition stats
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mstopdiff_config as cfg

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 13

XMAX = 90  # delta-mass axis upper limit (Da)
HERE = os.path.dirname(os.path.abspath(__file__))

# Curated common top-down delta masses (monoisotopic, Da): biological PTMs and
# the artifacts/adducts MSTopDiff is designed to surface.
MODS = [
    (0.984, "Deamidation"),
    (14.016, "Methylation"),
    (15.995, "Oxidation"),
    (18.011, "Water / dehydration"),  # Unimod Dehydrated; pyro-Glu(E)
    (21.982, "Na adduct"),
    (27.995, "Formylation"),
    (28.031, "Dimethyl / +C2H4"),
    (31.990, "Dioxidation"),
    (37.956, "K adduct"),
    (42.011, "Acetylation"),
    (42.047, "Trimethylation"),
    (43.006, "Carbamylation"),
    (43.990, "Carboxylation"),
    (44.026, "PEG (+C2H4O)"),     # polyethylene glycol contamination ladder
    (45.988, "Methylthio"),       # Unimod #39; Cys beta-methylthiolation
    (47.985, "Trioxidation"),
    (57.021, "Carbamidomethyl"),
    (58.005, "Carboxymethyl"),    # Unimod #6; iodoacetate alkylation
    (63.980, "Tetra-oxidation"),
    (71.037, "Acrylamide"),
    (75.998, "b-ME adduct"),      # 2-mercaptoethanol mixed disulfide
    (79.966, "Phosphorylation"),
    (79.957, "Sulfation"),
    (86.000, "Malonylation"),     # Unimod #747
    (100.016, "Succinylation"),
    (114.043, "GlyGly (Ub)"),
    (119.004, "Cysteinylation"),
    (128.059, "+Gln"),
    (128.095, "+Lys"),
    (162.053, "Hexose"),  # out of range, will be skipped
]

# Figure annotations: degenerate shifts (within ~0.05 Da) merged into one label.
ANNOT = [
    (0.984, "Isotopic error / +1"),
    (14.016, "Methylation"),
    (15.995, "Oxidation"),
    (21.982, "Na adduct"),
    (28.01, "Formyl / Dimethyl"),
    (30.011, "Me+Ox"),
    (31.990, "Di-oxidation"),
    (37.956, "K adduct"),
    (42.03, "Acetyl / Trimethyl"),
    (43.99, "Carboxyl / PEG"),
    (47.985, "Tri-oxidation"),
    (57.021, "Carbamidomethyl"),
    (63.980, "Tetra-oxidation"),
    (71.037, "Acrylamide"),
    (75.998, "b-ME adduct"),
    (79.96, "Phospho / Sulfation"),
    (100.016, "Succinylation"),
    (114.043, "GlyGly (Ub)"),
    (119.004, "Cysteinylation"),
    (128.08, "+Gln / +Lys"),
]

HALF_WIN = 0.05        # +/- window around a mod mass (Da)


def load(path):
    df = pd.read_csv(path)
    df["ixc"] = (df["intensity_lower_mass x count"].abs()
                 + df["intensity_higher_mass x count"].abs())
    return df


def rebin(df, factor=10):
    """Sum 0.01 Da bins into coarser bins (factor=10 -> 0.1 Da).

    Returns bin starts, count, and the signed intensity x count partners
    (lower >= 0, higher <= 0) for the faithful diverging MSTopDiff plot.
    """
    n = (len(df) // factor) * factor
    rs = lambda s: df[s].to_numpy()[:n].reshape(-1, factor).sum(1)
    binc = df["bin"].to_numpy()[:n].reshape(-1, factor)[:, 0]
    cnt = rs("count")
    ixc_lo = rs("intensity_lower_mass x count")   # >= 0
    ixc_hi = rs("intensity_higher_mass x count")  # <= 0
    return binc, cnt, ixc_lo, ixc_hi


def window_peak(df, mass, col):
    """Tallest bin within +/-HALF_WIN of a modification mass."""
    x = df["bin"].to_numpy()
    y = df[col].to_numpy().astype(float)
    in_win = np.abs(x - mass) <= HALF_WIN
    return float(y[in_win].max()) if in_win.any() else 0.0


def base_peak(df, col, xmax):
    """Tallest bin within 0..xmax Da (the base peak for relative abundance)."""
    x = df["bin"].to_numpy()
    y = df[col].to_numpy().astype(float)
    return float(y[x <= xmax].max())


def stacked_figure(dsets, outname, height_per=3.4):
    """Stacked single-sided relative-intensity histogram, one row per dataset."""
    n = len(dsets)
    fig, axes = plt.subplots(n, 1, figsize=(11, max(3, height_per * n)),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, d in zip(axes, dsets):
        bx, cnt, ixc_lo, ixc_hi = rebin(d["df"])
        peak = ixc_lo[bx <= XMAX].max()
        ax.plot(bx, ixc_lo / peak if peak else ixc_lo, color=d["color"], lw=0.8)
        ax.set_xlim(0, XMAX)
        ax.set_ylim(0, 1.05)
        for m, nm in ANNOT:
            if m <= XMAX:
                ax.axvline(m, color="grey", ls=":", lw=0.5, alpha=0.5)
        ax.set_ylabel("rel. intensity x count")
        ax.text(0.995, 0.92, d["label"], transform=ax.transAxes, ha="right",
                va="top", fontsize=14, color=d["color"], fontweight="bold")
    ymax = axes[0].get_ylim()[1]
    for m, nm in ANNOT:
        if m <= XMAX:
            axes[0].text(m, ymax, nm, rotation=45, va="bottom", ha="left",
                         fontsize=8.5, alpha=0.85)
    axes[-1].set_xlabel("delta mass (Da)")
    fig.suptitle("MSTopDiff intensity x count delta-mass histograms "
                 "(confident masses, 0-60 min, rebinned to 0.1 Da)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"{outname}.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    dsets = cfg.discover()
    if not dsets:
        raise SystemExit("No *_conf_mstopdiff.csv found in flashdeconv/mstopdiff/")
    for d in dsets:
        d["df"] = load(d["csv"])
        d["base"] = base_peak(d["df"], "ixc", XMAX)
        d["total"] = int(d["df"]["count"].sum())
    mods = [(m, n) for m, n in MODS if m <= XMAX]

    # ---------- relative-abundance table (all datasets) ----------
    rows = []
    for mass, name in mods:
        row = {"modification": name, "mod_mass": mass}
        for d in dsets:
            rel = window_peak(d["df"], mass, "ixc") / d["base"] if d["base"] else np.nan
            row[f"{d['label']}|ixc_rel"] = round(rel, 4)
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(HERE, "mstopdiff_mod_table.csv"), index=False)

    # ---------- stacked histograms: publication subset + all datasets ----------
    pub = [d for d in dsets if d["is_pub"]]
    if pub:
        stacked_figure(pub, "fig_mstopdiff_ixc")
    stacked_figure(dsets, "fig_mstopdiff_ixc_all")

    # ---------- relative-abundance bar chart (publication subset) ----------
    if pub:
        fig2, ax = plt.subplots(figsize=(10, 11))
        names = table["modification"] + "  (" + table["mod_mass"].map(
            lambda v: f"{v:g}") + ")"
        yp = np.arange(len(table))
        h = 0.8 / len(pub)
        for j, d in enumerate(pub):
            vals = 100 * table[f"{d['label']}|ixc_rel"]
            ax.barh(yp + (len(pub) - 1 - 2 * j) * h / 2, vals, height=h,
                    color=d["color"], label=d["label"])
        ax.set_yticks(yp); ax.set_yticklabels(names, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("relative abundance (% of base peak, intensity x count)")
        ax.set_title("Common modification abundance by cleanup condition\n"
                     "(tallest intensity x count bin within +/-0.05 Da)")
        ax.legend(fontsize=8, loc="lower right")
        fig2.tight_layout()
        for ext in ("png", "pdf"):
            fig2.savefig(os.path.join(HERE, f"fig_mstopdiff_mod_enrichment.{ext}"),
                         dpi=200, bbox_inches="tight")
        plt.close(fig2)

    # ---------- relative-abundance heatmap (all datasets) ----------
    labels = [d["label"] for d in dsets]
    M = np.array([[100 * table.loc[i, f"{lab}|ixc_rel"] for lab in labels]
                  for i in range(len(table))])
    fig3, ax = plt.subplots(figsize=(1.1 * len(labels) + 4, 0.42 * len(table) + 2))
    im = ax.imshow(M, aspect="auto", cmap="magma_r", vmin=0, vmax=100)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(table)))
    ax.set_yticklabels([f"{n}  ({m:g})" for n, m in
                        zip(table["modification"], table["mod_mass"])], fontsize=8)
    for i in range(len(table)):
        for j in range(len(labels)):
            if M[i, j] >= 5:
                ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center",
                        fontsize=6, color="white" if M[i, j] > 55 else "black")
    fig3.colorbar(im, ax=ax, label="rel. abundance (% base peak)", shrink=0.6)
    ax.set_title("Modification relative abundance across all datasets")
    fig3.tight_layout()
    for ext in ("png", "pdf"):
        fig3.savefig(os.path.join(HERE, f"fig_mstopdiff_mod_heatmap.{ext}"),
                     dpi=200, bbox_inches="tight")
    plt.close(fig3)

    # ---------- console ----------
    print(f"Datasets ({len(dsets)}):",
          ", ".join(f"{d['label']}(n={d['total']})" for d in dsets))
    print(f"Publication subset: {[d['label'] for d in pub]}")
    disp = table[["modification", "mod_mass"]
                 + [f"{lab}|ixc_rel" for lab in labels]].copy()
    disp.columns = ["modification", "mass"] + labels
    with pd.option_context("display.width", 220):
        print(disp.to_string(index=False))


if __name__ == "__main__":
    main()
