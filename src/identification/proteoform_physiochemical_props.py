"""
proteoform_physiochemical_props.py

Calculate physiochemical properties for FDR-confident proteoforms pulled DIRECTLY
from a ProSightPD .tdReport (a SQLite database), across ALL sample-cleanup
conditions, grouped by bead type.

This replaces the earlier CSV-driven workflow: proteoform sequences (and the
average mass) are now read straight from the ChemicalProteoform table of the
.tdReport instead of from a hand-exported TDreport hit report.

Decisions (mirroring shared_proteoforms.py):
  - proteoform identity   = ChemicalProteoformId (sequence + mods + mass)
  - "identified"          = has a PrSM-level (AggregationLevel 0) hit with
                            GlobalQvalue <= 0.01 in that condition's file(s)
  - conditions compared   = every included file in experimental_design.csv,
                            grouped by bead type (MCW / MagReSyn / Cytiva)

Panel-A figure: average mass (kDa), isoelectric point (pI) and GRAVY, stacked
vertically with one violin per condition, grouped by bead type. Violins (KDE) are
used instead of box plots because pI is strongly bimodal, which a box plot hides.

Histones dominate this dataset, so every summary/figure is produced twice: once
for ALL confident proteoforms and once EXCLUDING histones (see is_histone()).

GRAVY / pI depend only on the bare amino-acid sequence (modifications, stored
separately in ModificationHash, are ignored); average mass is the proteoform's
AverageMass from the report (includes modifications). Selenocysteine (U) has no
entry in Biopython's hydropathy / instability tables, so it is substituted with C
for the property calculation ONLY; a ContainsU flag marks affected proteoforms.

Outputs (written next to this script):
  proteoform_physiochemical_props.csv        one row per unique proteoform
  proteoform_physiochemical_props_long.csv   tidy, one row per (proteoform x condition)
  physiochemical_summary.csv                 per-condition x per-property stats (all)
  physiochemical_summary_nohistone.csv       per-condition x per-property stats (no histones)
  fig_physiochemical_props.{pdf,png}         stacked panel A (all proteoforms)
  fig_physiochemical_props_nohistone.{pdf,png}  stacked panel A, histones excluded

Run:  conda run -n tdms python proteoform_physiochemical_props.py
"""

from __future__ import annotations
import re
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Repo layout: this script lives in src/identification/, so the project root is
# two levels up. Inputs live under data/ and config/, outputs under results/.
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
DB = DATA_DIR / "20250723_ifeltens_SP3_TDP_annotated_subsequence_2.tdReport"
DESIGN = ROOT / "config" / "experimental_design.csv"
FIGURES_DIR = ROOT / "results" / "figures"
TABLES_DIR = ROOT / "results" / "tables"

# PrSM-level (AggregationLevel 0) q-value threshold for a confident identification.
PRSM_QVALUE_MAX = 0.01

# All conditions, grouped by bead type. Tuple fields:
#   (Cleanup, Resuspension, BeadGroup, x-axis label, colour)
# Order here = left-to-right on the figure. Colours follow the manuscript scheme
# (MCW green / MagReSyn red / Cytiva blue).
CONDITIONS = [
    ("MCW",               "5% ACN, 0.1% FA", "MCW",      "",         "#5AAF46"),
    ("MagReSyn Hydroxyl", "0.5% TFA",        "MagReSyn", "0.5% TFA", "#ED645A"),
    ("MagReSyn Hydroxyl", "2% TFA",          "MagReSyn", "2% TFA",   "#ED645A"),
    ("MagReSyn Hydroxyl", "10% FA",          "MagReSyn", "10% FA",   "#ED645A"),
    ("MagReSyn Hydroxyl", "20% FA",          "MagReSyn", "20% FA",   "#ED645A"),
    ("Cytiva Carboxyl",   "0.5% TFA",        "Cytiva",   "0.5% TFA", "#6482CD"),
    ("Cytiva Carboxyl",   "2% TFA",          "Cytiva",   "2% TFA",   "#6482CD"),
    ("Cytiva Carboxyl",   "10% FA",          "Cytiva",   "10% FA",   "#6482CD"),
    ("Cytiva Carboxyl",   "20% FA",          "Cytiva",   "20% FA",   "#6482CD"),
]

# Mass panel uses MonoisotopicMass: the DB's AverageMass column is only ~40%
# populated (placeholder 18.0153 Da otherwise), whereas MonoisotopicMass is
# complete. Average vs monoisotopic differ by ~0.06% at these sizes (immaterial
# on a kDa axis), so the mass values are monoisotopic but labelled "average mass".
# Properties summarised in the CSV (output column, label).
SUMMARY_PROPERTIES = [
    ("Mass_kDa",         "Average mass (kDa)"),
    ("GRAVY",            "GRAVY score"),
    ("IsoelectricPoint", "Isoelectric point (pI)"),
    ("Aromaticity",      "Aromaticity"),
    ("Instability",      "Instability index"),
]
# Properties drawn in the stacked panel-A figure, top-to-bottom.
FIGURE_PROPERTIES = [
    ("Mass_kDa",         "Average mass (kDa)"),
    ("IsoelectricPoint", "Isoelectric point (pI)"),
    ("GRAVY",            "GRAVY score"),
]

# A description belongs to a histone if a comma-delimited entry *starts* with
# "Histone". This deliberately excludes "Non-histone chromosomal protein HMG-..".
HISTONE_RE = re.compile(r"(?:^|,)\s*histone", re.IGNORECASE)

FIG_EXTS = ("pdf", "png")  # vector PDF deliverable + PNG preview


def style_matplotlib() -> None:
    """Arial / manuscript-consistent matplotlib defaults."""
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False


def savefig(fig, stem) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in FIG_EXTS:
        fig.savefig(FIGURES_DIR / f"{stem}.{ext}", dpi=150, bbox_inches="tight")


def is_histone(description) -> bool:
    return isinstance(description, str) and bool(HISTONE_RE.search(description))


def build_condition_meta(con: sqlite3.Connection) -> pd.DataFrame:
    """Resolve CONDITIONS against the design + DataFile table and assign x-axis
    positions (with a gap between bead groups).

    Returns one row per condition: Condition, BeadGroup, XLabel, Colour,
    CondOrder, Position, plus a DataFileId list per condition."""
    design = pd.read_csv(DESIGN)
    design["Include"] = design["Include"].astype(str).str.strip().str.upper().eq("TRUE")
    files = pd.read_sql_query("SELECT Id AS DataFileId, Name FROM DataFile", con)
    design = design.merge(files, on="Name", how="left")

    rows = []
    for order, (cleanup, resus, group, xlabel, colour) in enumerate(CONDITIONS):
        cond = f"{group}" if group == "MCW" else f"{group} {xlabel}"
        match = design[(design["Cleanup"] == cleanup)
                       & (design["Resuspension"] == resus)
                       & design["Include"]]
        if match.empty:
            raise SystemExit(f"No included file for condition: {cleanup} / {resus}")
        if match["DataFileId"].isna().any():
            raise SystemExit(f"DataFile not found for: {match['Name'].tolist()}")
        rows.append({
            "Condition": cond, "BeadGroup": group, "XLabel": xlabel,
            "Colour": colour, "CondOrder": order,
            "DataFileIds": match["DataFileId"].astype(int).tolist(),
            "FileNames": match["Name"].tolist(),
        })

    meta = pd.DataFrame(rows)

    # x positions: 1.0 between conditions, +1.0 extra gap between bead groups
    positions, pos, prev = [], 0.0, None
    for grp in meta["BeadGroup"]:
        if prev is not None and grp != prev:
            pos += 1.0
        positions.append(pos)
        pos += 1.0
        prev = grp
    meta["Position"] = positions
    return meta


def load_confident_proteoforms(con: sqlite3.Connection,
                               file_ids: list[int]) -> pd.DataFrame:
    """One row per (ChemicalProteoformId, DataFileId) that has at least one
    FDR-passing PrSM hit in that file. Returns: ChemicalProteoformId, DataFileId."""
    placeholders = ",".join("?" for _ in file_ids)
    q = f"""
        SELECT DISTINCT h.ChemicalProteoformId AS ChemicalProteoformId,
                        h.DataFileId           AS DataFileId
        FROM Hit h
        JOIN GlobalQualitativeConfidence g
          ON g.HitId = h.Id AND g.AggregationLevel = 0
        WHERE g.GlobalQvalue <= ?
          AND h.DataFileId IN ({placeholders})
    """
    params = [PRSM_QVALUE_MAX, *file_ids]
    return pd.read_sql_query(q, con, params=params)


def load_annotation(con: sqlite3.Connection) -> pd.DataFrame:
    """One row per ChemicalProteoformId: sequence, masses, protein annotation."""
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


def calculate_properties(sequence: str) -> tuple[float, float, float, float]:
    """GRAVY, isoelectric point, aromaticity, instability index for a sequence.

    Selenocysteine (U) is substituted with C: Biopython's Kyte-Doolittle and
    instability (DIWV) tables have no U, and Sec is the selenium analogue of Cys.
    The substitution affects the property math only, not the stored sequence.
    """
    analysed = ProteinAnalysis(sequence.replace("U", "C"))
    return (analysed.gravy(),
            analysed.isoelectric_point(),
            analysed.aromaticity(),
            analysed.instability_index())


def summarize(long: pd.DataFrame, cond_order: list[str]) -> pd.DataFrame:
    """Per-condition x per-property summary statistics from a tidy long table."""
    rows = []
    for prop, _ in SUMMARY_PROPERTIES:
        s = (long.groupby("Condition", observed=True)[prop]
             .agg(n="count", mean="mean", median="median", sd="std",
                  min="min", max="max")
             .reindex(cond_order).round(4))
        s.insert(0, "Property", prop)
        rows.append(s)
    return pd.concat(rows)


def make_figure(long: pd.DataFrame, meta: pd.DataFrame, title: str, stem: str) -> None:
    """Stacked panel A: one row per FIGURE_PROPERTIES entry, one violin per
    condition, conditions grouped by bead type along x."""
    positions = meta["Position"].to_numpy()
    colours = meta["Colour"].tolist()
    conditions = meta["Condition"].tolist()
    rng = np.random.default_rng(0)

    n_rows = len(FIGURE_PROPERTIES)
    fig, axes = plt.subplots(n_rows, 1, sharex=True, figsize=(7.5, 3 * n_rows))

    for ax, (prop, label) in zip(axes, FIGURE_PROPERTIES):
        data = [long.loc[long["Condition"] == c, prop].dropna().values
                for c in conditions]
        parts = ax.violinplot(data, positions=positions, widths=0.85,
                              showmedians=True, showextrema=False)
        for body, col in zip(parts["bodies"], colours):
            body.set_facecolor(col); body.set_alpha(0.6)
            body.set_edgecolor("0.3"); body.set_linewidth(0.8)
        parts["cmedians"].set_color("0.15"); parts["cmedians"].set_linewidth(1.4)
        for xi, vals in zip(positions, data):
            q1, q3 = np.percentile(vals, [25, 75])
            ax.vlines(xi, q1, q3, color="0.15", lw=3.5, alpha=0.45, zorder=2)
        for xi, vals, col in zip(positions, data, colours):
            jit = rng.uniform(-0.09, 0.09, size=len(vals))
            ax.scatter(xi + jit, vals, s=3, color=col, alpha=0.2,
                       edgecolors="none", zorder=3)
        ax.set_ylabel(label)
        ax.grid(axis="y", color="0.92", lw=0.7)
        ax.set_axisbelow(True)
        ax.margins(x=0.02)

    bottom = axes[-1]
    bottom.set_xticks(positions)
    bottom.set_xticklabels(meta["XLabel"].tolist(), rotation=45, ha="right")

    # bead-group brackets + labels below the bottom axis
    trans = bottom.get_xaxis_transform()
    for group, sub in meta.groupby("BeadGroup", sort=False):
        x0, x1 = sub["Position"].min(), sub["Position"].max()
        bottom.plot([x0 - 0.4, x1 + 0.4], [-0.30, -0.30], transform=trans,
                    color="0.2", lw=1.2, clip_on=False)
        bottom.text((x0 + x1) / 2, -0.34, group, transform=trans,
                    ha="center", va="top", fontweight="bold", clip_on=False)

    fig.suptitle(title, y=0.995, fontsize=11)
    fig.tight_layout()
    savefig(fig, stem)
    plt.close(fig)


def report_and_render(long: pd.DataFrame, meta: pd.DataFrame, tag: str,
                      title: str, summary_path: Path, fig_stem: str) -> None:
    """Write a summary CSV + figure for one proteoform subset and print its report."""
    cond_order = meta["Condition"].tolist()
    n_pf = long["ChemicalProteoformId"].nunique()
    summary = summarize(long, cond_order)
    summary.to_csv(summary_path)

    print("\n" + "=" * 78)
    print(f"[{tag}]  unique proteoforms: {n_pf}")
    print("=" * 78)
    counts = long.groupby("Condition", observed=True)["ChemicalProteoformId"].nunique()
    print("Identifications per condition:")
    for c in cond_order:
        print(f"  {c:20s}: {int(counts.get(c, 0)):4d}")
    print("\nPer-condition property summary:\n")
    print(summary.to_string())

    make_figure(long, meta, title, fig_stem)


def load_dataset(con: sqlite3.Connection) -> tuple[pd.DataFrame, pd.DataFrame,
                                                   pd.DataFrame, pd.DataFrame]:
    """Assemble the FDR-confident proteoform dataset from the .tdReport.

    Returns (meta, props, pf_cond, long):
      meta    one row per condition (BeadGroup, XLabel, Colour, Position, files)
      props   one row per unique proteoform (sequence, masses, properties, flags)
      pf_cond one row per (proteoform, condition) it was confidently identified in
      long    pf_cond joined to props (tidy table for plotting / statistics)
    Shared by the figure script and physiochemical_stats.py so both use the exact
    same identification set and property values.
    """
    meta = build_condition_meta(con)
    all_file_ids = sorted({fid for ids in meta["DataFileIds"] for fid in ids})
    ident = load_confident_proteoforms(con, all_file_ids)
    annot = load_annotation(con)

    # map DataFileId -> Condition, then to one row per (proteoform, condition)
    file_to_cond = {fid: r["Condition"]
                    for _, r in meta.iterrows() for fid in r["DataFileIds"]}
    ident["Condition"] = ident["DataFileId"].map(file_to_cond)
    pf_cond = ident[["ChemicalProteoformId", "Condition"]].dropna().drop_duplicates()

    confident_ids = pf_cond["ChemicalProteoformId"].unique()
    if len(confident_ids) == 0:
        raise SystemExit("No FDR-confident proteoforms found for the chosen conditions.")

    props = annot[annot["ChemicalProteoformId"].isin(confident_ids)].copy()
    props["Mass_kDa"] = props["MonoisotopicMass"] / 1000.0
    props["Length"] = props["Sequence"].str.len()
    props["ContainsU"] = props["Sequence"].str.contains("U")
    props["IsHistone"] = props["Description"].apply(is_histone)
    props[["GRAVY", "IsoelectricPoint", "Aromaticity", "Instability"]] = (
        props["Sequence"].apply(lambda s: pd.Series(calculate_properties(s))))

    cond_order = meta["Condition"].tolist()
    prop_cols = [c for c, _ in SUMMARY_PROPERTIES]
    long = pf_cond.merge(
        props[["ChemicalProteoformId", "Accession", "Description", "IsHistone",
               "Length", "ContainsU"] + prop_cols],
        on="ChemicalProteoformId", how="left")
    long["Condition"] = pd.Categorical(long["Condition"], categories=cond_order, ordered=True)
    long = long.sort_values(["Condition", "ChemicalProteoformId"])
    return meta, props, pf_cond, long


def main() -> None:
    if not DB.exists():
        raise SystemExit(f"Database not found: {DB}")
    if not DESIGN.exists():
        raise SystemExit(f"Design file not found: {DESIGN}")

    style_matplotlib()
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    try:
        meta, props, pf_cond, long = load_dataset(con)
    finally:
        con.close()

    print("Conditions (left-to-right) and their files:")
    for _, r in meta.iterrows():
        print(f"  {r['Condition']:20s}: {', '.join(r['FileNames'])}")

    confident_ids = props["ChemicalProteoformId"].unique()
    cond_order = meta["Condition"].tolist()
    prop_cols = [c for c, _ in SUMMARY_PROPERTIES]

    # presence matrix (one column per condition) + NumConditions
    membership = (pf_cond.assign(present=1)
                  .pivot_table(index="ChemicalProteoformId", columns="Condition",
                               values="present", fill_value=0)
                  .reindex(columns=cond_order, fill_value=0).astype(int))
    membership["NumConditions"] = membership.sum(axis=1)

    front = ["ChemicalProteoformId", "Accession", "Description", "IsHistone",
             "Mass_kDa", "MonoisotopicMass", "Length", "ContainsU", "Sequence"]
    wide = (props.set_index("ChemicalProteoformId")[front[1:] + prop_cols[1:]]
            .join(membership, how="right").reset_index()
            .sort_values(["NumConditions", "ChemicalProteoformId"],
                         ascending=[False, True]))
    wide[front + prop_cols[1:] + cond_order + ["NumConditions"]].to_csv(
        TABLES_DIR / "proteoform_physiochemical_props.csv", index=False)

    # ---- tidy / long: one row per (proteoform x condition) ------------------
    long.to_csv(TABLES_DIR / "proteoform_physiochemical_props_long.csv", index=False)

    n_hist = int(props["IsHistone"].sum())
    n_u = int(props["ContainsU"].sum())
    print(f"\nFDR-confident proteoforms (q<= {PRSM_QVALUE_MAX}): {len(confident_ids)}"
          f"  ({n_hist} histone, {len(confident_ids) - n_hist} non-histone;"
          f" {n_u} with U->C)")

    # ---- report + figures: all proteoforms, then histones excluded ----------
    report_and_render(
        long, meta, tag="all proteoforms",
        title=f"Proteoform physiochemical properties by condition "
              f"(n={len(confident_ids)} FDR-confident proteoforms)",
        summary_path=TABLES_DIR / "physiochemical_summary.csv",
        fig_stem="fig_physiochemical_props")

    long_nh = long[~long["IsHistone"]]
    report_and_render(
        long_nh, meta, tag="histones excluded",
        title=f"Proteoform physiochemical properties by condition, histones excluded "
              f"(n={long_nh['ChemicalProteoformId'].nunique()} proteoforms)",
        summary_path=TABLES_DIR / "physiochemical_summary_nohistone.csv",
        fig_stem="fig_physiochemical_props_nohistone")

    print("\nWrote: proteoform_physiochemical_props.csv, "
          "proteoform_physiochemical_props_long.csv,\n       "
          "physiochemical_summary.csv, physiochemical_summary_nohistone.csv,\n       "
          "fig_physiochemical_props.{pdf,png}, "
          "fig_physiochemical_props_nohistone.{pdf,png}")


if __name__ == "__main__":
    main()
