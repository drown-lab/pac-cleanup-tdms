"""
shared_proteoforms.py

Compare the kelleher_negLog_pScore of proteoforms shared across a chosen set of
sample-cleanup CONDITIONS in a ProSightPD .tdReport.

Sample background: one HeLa lysate prepared with different cleanup approaches.
The experimental design lives in experimental_design.csv:
    Name, Cleanup, Resuspension, Include
Files with Include == FALSE used a different acquisition method and are dropped.

Conditions compared here (per user request):
    MCW                 (MCW,             "5% ACN, 0.1% FA")
    Cytiva 0.5% TFA     (Cytiva Carboxyl, "0.5% TFA")
    MagReSyn 0.5% TFA   (MagReSyn Hydroxyl,"0.5% TFA")

Decisions:
  - proteoform identity   = ChemicalProteoformId (sequence + mods + mass)
  - per-condition score    = best (max) kelleher_negLog_pScore over that
                             proteoform's FDR-passing hits in that condition's
                             included file(s)
  - "shared"               = identified in EVERY compared condition

Outputs:
  shared_proteoforms_wide.csv     proteoform x condition matrix of max negLog P
  shared_proteoforms_long.csv     tidy (one row per proteoform x condition)
  per_sample_score_summary.csv    per-condition summary statistics
  fig_score_distribution.png      box/strip plot per condition
  fig_score_heatmap.png           clustered heatmap (proteoform x condition)
  fig_score_paired.png            paired lines (same proteoform across conditions)
  fig_venn_identifications.png    Venn of identified proteoform sets per condition
  identification_membership.csv   presence/absence of every proteoform per condition

Run:  conda run -n tdms python analysis_01_shared_proteoforms.py
"""

from __future__ import annotations
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from matplotlib_venn import venn3

HERE = Path(__file__).resolve().parent
DB = HERE / "20250721_ifeltens_BEC_Consensus.tdReport"
DESIGN = HERE / "experimental_design.csv"

# ScoreTypeId 2 == kelleher_negLog_pScore (see ScoreType table / tdReport_schema.md)
NEGLOG_PSCORE_SCORETYPE_ID = 2
# PrSM-level (AggregationLevel 0) q-value threshold (all hits already pass 0.01).
PRSM_QVALUE_MAX = 0.01

# Conditions to compare: (Cleanup, Resuspension, display label).  Order = x-axis.
CONDITIONS_OF_INTEREST = [
    ("MCW",               "5% ACN, 0.1% FA", "MCW"),
    ("Cytiva Carboxyl",   "0.5% TFA",        "Cytiva 0.5% TFA"),
    ("MagReSyn Hydroxyl", "0.5% TFA",        "MagReSyn 0.5% TFA"),
]

FIG_EXTS = ("pdf", "png")  # vector PDF deliverable + PNG preview


def savefig(obj, stem):
    """Save a Matplotlib Figure or seaborn grid to all FIG_EXTS next to this script."""
    for ext in FIG_EXTS:
        obj.savefig(HERE / f"{stem}.{ext}", dpi=150)


def load_selected_files(con: sqlite3.Connection) -> pd.DataFrame:
    """Read the design, drop Include==FALSE, keep only conditions of interest, and
    attach each row's DataFile.Id. Returns columns: DataFileId, Name, Condition."""
    design = pd.read_csv(DESIGN)
    design["Include"] = design["Include"].astype(str).str.strip().str.upper().eq("TRUE")

    cond = pd.DataFrame(CONDITIONS_OF_INTEREST,
                        columns=["Cleanup", "Resuspension", "Condition"])
    cond["CondOrder"] = range(len(cond))

    sel = (design[design["Include"]]
           .merge(cond, on=["Cleanup", "Resuspension"], how="inner"))

    files = pd.read_sql_query("SELECT Id AS DataFileId, Name FROM DataFile", con)
    sel = sel.merge(files, on="Name", how="left")
    missing = sel[sel["DataFileId"].isna()]
    if not missing.empty:
        raise SystemExit(f"Design names not found in DataFile: "
                         f"{missing['Name'].tolist()}")
    return sel.sort_values("CondOrder")[["DataFileId", "Name", "Condition", "CondOrder"]]


def load_per_file_scores(con: sqlite3.Connection, file_ids: list[int]) -> pd.DataFrame:
    """One row per (ChemicalProteoformId, DataFileId): max negLog P-score across
    that proteoform's FDR-passing hits in that file."""
    placeholders = ",".join("?" for _ in file_ids)
    q = f"""
        SELECT h.ChemicalProteoformId AS ChemicalProteoformId,
               h.DataFileId           AS DataFileId,
               MAX(hs.Value)          AS NegLogPScore,
               COUNT(*)               AS NumHits
        FROM Hit h
        JOIN HitScore hs
          ON hs.HitId = h.Id AND hs.ScoreTypeId = ?
        JOIN GlobalQualitativeConfidence g
          ON g.HitId = h.Id AND g.AggregationLevel = 0
        WHERE g.GlobalQvalue <= ?
          AND h.DataFileId IN ({placeholders})
        GROUP BY h.ChemicalProteoformId, h.DataFileId
    """
    params = [NEGLOG_PSCORE_SCORETYPE_ID, PRSM_QVALUE_MAX, *file_ids]
    return pd.read_sql_query(q, con, params=params)


def load_annotation(con: sqlite3.Connection) -> pd.DataFrame:
    """One row per ChemicalProteoformId with protein annotation + chemistry."""
    q = """
        SELECT cp.Id               AS ChemicalProteoformId,
               cp.MonoisotopicMass AS MonoisotopicMass,
               cp.Sequence         AS Sequence,
               cp.ModificationHash AS ModificationHash,
               GROUP_CONCAT(DISTINCT e.AccessionNumber) AS Accession,
               GROUP_CONCAT(DISTINCT e.Description)      AS Description
        FROM ChemicalProteoform cp
        LEFT JOIN BiologicalProteoform bp ON bp.ChemicalProteoformId = cp.Id
        LEFT JOIN Isoform iso             ON iso.Id = bp.IsoformId
        LEFT JOIN Entry e                 ON e.Id = iso.EntryId
        GROUP BY cp.Id
    """
    return pd.read_sql_query(q, con)


def main() -> None:
    if not DB.exists():
        raise SystemExit(f"Database not found: {DB}")
    if not DESIGN.exists():
        raise SystemExit(f"Design file not found: {DESIGN}")

    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    try:
        sel = load_selected_files(con)
        scores = load_per_file_scores(con, sel["DataFileId"].astype(int).tolist())
        annot = load_annotation(con)
    finally:
        con.close()

    cond_order = (sel.drop_duplicates("Condition").sort_values("CondOrder")
                  ["Condition"].tolist())
    n_cond = len(cond_order)

    print("Included files per condition:")
    for c in cond_order:
        names = sel.loc[sel["Condition"] == c, "Name"].tolist()
        print(f"  {c:18s}: {', '.join(names)}")
    print()

    # attach condition, then collapse to one score per (proteoform, condition)
    scores = scores.merge(sel[["DataFileId", "Condition"]], on="DataFileId", how="left")
    per_cond = (scores.groupby(["ChemicalProteoformId", "Condition"], as_index=False)
                .agg(NegLogPScore=("NegLogPScore", "max"),
                     NumHits=("NumHits", "sum")))

    # ---- select proteoforms present in EVERY compared condition -------------
    presence = per_cond.groupby("ChemicalProteoformId")["Condition"].nunique()
    shared_ids = presence[presence == n_cond].index
    shared = per_cond[per_cond["ChemicalProteoformId"].isin(shared_ids)].copy()

    # ---- identification sets per condition (for Venn / membership) -----------
    id_sets = {c: set(per_cond.loc[per_cond["Condition"] == c, "ChemicalProteoformId"])
               for c in cond_order}

    # presence/absence membership table over every proteoform seen in >=1 cond
    membership = (per_cond.assign(present=1)
                  .pivot_table(index="ChemicalProteoformId", columns="Condition",
                               values="present", fill_value=0)
                  .reindex(columns=cond_order).astype(int))
    membership["NumConditions"] = membership.sum(axis=1)
    membership = (annot.set_index("ChemicalProteoformId")
                  [["Accession", "Description", "MonoisotopicMass"]]
                  .join(membership, how="right").reset_index()
                  .sort_values(["NumConditions", "ChemicalProteoformId"],
                               ascending=[False, True]))
    membership.to_csv(HERE / "identification_membership.csv", index=False)

    # ---- wide matrix: proteoform x condition --------------------------------
    wide = (shared.pivot(index="ChemicalProteoformId", columns="Condition",
                         values="NegLogPScore")
            .reindex(columns=cond_order))
    wide = annot.set_index("ChemicalProteoformId").join(wide, how="right").reset_index()
    wide = wide.sort_values("ChemicalProteoformId")

    front = ["ChemicalProteoformId", "Accession", "Description",
             "MonoisotopicMass", "ModificationHash"]
    wide[front + cond_order].to_csv(HERE / "shared_proteoforms_wide.csv", index=False)

    # ---- tidy / long --------------------------------------------------------
    long = shared.merge(
        annot[["ChemicalProteoformId", "Accession", "Description", "MonoisotopicMass"]],
        on="ChemicalProteoformId", how="left")
    long["Condition"] = pd.Categorical(long["Condition"], categories=cond_order, ordered=True)
    long = long.sort_values(["ChemicalProteoformId", "Condition"])
    long.to_csv(HERE / "shared_proteoforms_long.csv", index=False)

    # ---- per-condition summary ----------------------------------------------
    summary = (long.groupby("Condition", observed=True)["NegLogPScore"]
               .agg(n="count", mean="mean", median="median", sd="std",
                    min="min", max="max")
               .reindex(cond_order).round(3))
    summary.to_csv(HERE / "per_sample_score_summary.csv")

    # ---- statistics ---------------------------------------------------------
    mat = wide.set_index("ChemicalProteoformId")[cond_order]
    friedman_stat, friedman_p = stats.friedmanchisquare(*[mat[c].values for c in cond_order])
    # pairwise paired Wilcoxon signed-rank (same proteoforms across conditions)
    pairwise = []
    for i in range(n_cond):
        for j in range(i + 1, n_cond):
            a, b = cond_order[i], cond_order[j]
            w, p = stats.wilcoxon(mat[a].values, mat[b].values)
            pairwise.append((a, b, round(float(w), 1), p))

    # ---- report -------------------------------------------------------------
    print("=" * 70)
    print(f"Conditions compared                   : {n_cond}  ({', '.join(cond_order)})")
    print(f"Proteoforms in >=1 compared condition : {presence.shape[0]}")
    print(f"Proteoforms shared by ALL conditions  : {len(shared_ids)}")
    print("=" * 70)
    print("\nPer-condition max negLog P-score over the shared proteoforms:\n")
    print(summary.to_string())
    print(f"\nFriedman test (n={len(shared_ids)} paired): "
          f"chi2 = {friedman_stat:.2f}, p = {friedman_p:.3e}")
    print("\nPairwise Wilcoxon signed-rank:")
    for a, b, w, p in pairwise:
        print(f"  {a:18s} vs {b:18s}  W = {w:8.1f}  p = {p:.3e}")
    print("\nIdentification overlap (proteoform sets per condition):")
    only = {c: len(id_sets[c] - set().union(*(id_sets[o] for o in cond_order if o != c)))
            for c in cond_order}
    for c in cond_order:
        print(f"  {c:18s}: {len(id_sets[c]):4d} identified   ({only[c]} unique to it)")
    print(f"  {'shared by all 3':18s}: {len(shared_ids):4d}")

    print("\nWrote: shared_proteoforms_wide.csv, shared_proteoforms_long.csv, "
          "per_sample_score_summary.csv,\n       identification_membership.csv, "
          "fig_score_distribution.png, fig_score_heatmap.png,\n       "
          "fig_score_paired.png, fig_venn_identifications.png")

    # ---- figures ------------------------------------------------------------
    sns.set_theme(style="whitegrid")
    palette = sns.color_palette("Set2", n_cond)

    # (1) distribution per condition
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(data=long, x="Condition", y="NegLogPScore", order=cond_order,
                hue="Condition", palette=palette, legend=False, dodge=False,
                fliersize=0, ax=ax)
    sns.stripplot(data=long, x="Condition", y="NegLogPScore", order=cond_order,
                  jitter=0.25, size=4, alpha=0.5, color="0.2", ax=ax)
    ax.set_title(f"kelleher_negLog_pScore across conditions "
                 f"({len(shared_ids)} shared proteoforms)")
    ax.set_xlabel("Condition"); ax.set_ylabel("max negLog P-score per proteoform")
    fig.tight_layout(); savefig(fig, "fig_score_distribution")
    plt.close(fig)

    # (2) clustered heatmap
    g = sns.clustermap(mat, col_cluster=False, cmap="viridis",
                       figsize=(6, 11), yticklabels=False,
                       cbar_kws={"label": "max negLog P-score"})
    g.ax_heatmap.set_xlabel("Condition")
    g.ax_heatmap.set_ylabel(f"{len(shared_ids)} shared proteoforms")
    g.fig.suptitle("Shared proteoforms: negLog P-score by condition", y=1.01)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=30, ha="right")
    savefig(g, "fig_score_heatmap"); plt.close(g.fig)

    # (3) paired lines (one grey line per proteoform) + condition means
    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.arange(n_cond)
    for _, row in mat.iterrows():
        ax.plot(x, row.values, color="0.7", alpha=0.35, lw=0.7, zorder=1)
    ax.plot(x, mat.mean().values, color="#c0392b", lw=2.5, marker="o",
            zorder=3, label="mean")
    ax.plot(x, mat.median().values, color="#2c3e50", lw=2.0, marker="s",
            ls="--", zorder=3, label="median")
    ax.set_xticks(x); ax.set_xticklabels(cond_order, rotation=30, ha="right")
    ax.set_ylabel("max negLog P-score per proteoform")
    ax.set_title(f"Per-proteoform score across conditions (n={len(shared_ids)})")
    ax.legend()
    fig.tight_layout(); savefig(fig, "fig_score_paired")
    plt.close(fig)

    # (4) Venn of identified proteoform sets per condition
    if n_cond == 3:
        fig, ax = plt.subplots(figsize=(8, 7))
        v = venn3([id_sets[c] for c in cond_order], set_labels=cond_order,
                  set_colors=tuple(palette), alpha=0.55, ax=ax)
        for txt in (v.set_labels or []):
            if txt:
                txt.set_fontsize(11); txt.set_fontweight("bold")
        total = len(set().union(*id_sets.values()))
        ax.set_title(f"Identified proteoforms by condition "
                     f"({total} total, {len(shared_ids)} shared by all 3)")
        fig.tight_layout(); savefig(fig, "fig_venn_identifications")
        plt.close(fig)


if __name__ == "__main__":
    main()
