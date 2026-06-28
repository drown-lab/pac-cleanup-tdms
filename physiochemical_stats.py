"""
physiochemical_stats.py

Test whether proteoform physiochemical property distributions (average mass, pI,
GRAVY) vary by cleanup METHOD, using the same FDR-confident dataset as
proteoform_physiochemical_props.py (imported, so the identification set and
property values are identical).

Groupings tested:
  - bead type : MCW vs MagReSyn vs Cytiva (resuspension buffers pooled; each
                proteoform counted once per bead type)
  - condition : all nine cleanup + resuspension conditions

Each grouping is tested on ALL proteoforms and with HISTONES EXCLUDED, for each
of the three properties.

IMPORTANT caveat: a proteoform's property value is a FIXED attribute of its
sequence; methods differ only in WHICH proteoforms they detect, and the detected
sets overlap heavily. So observations across groups are NOT independent and these
tests are descriptive. With several hundred per group, p-values get tiny for even
trivial differences -- read the EFFECT SIZES (eta^2, Cliff's delta), not p.

Tests (scipy only):
  - Kruskal-Wallis omnibus across groups, with eta-squared effect size
  - pairwise Mann-Whitney U (two-sided) with Cliff's delta effect size
  - pairwise two-sample Kolmogorov-Smirnov (distribution shape, sensitive to the
    bimodality of pI)
  Pairwise p-values are Benjamini-Hochberg corrected within each
  grouping x subset x property family.

Outputs:
  physiochemical_stats_omnibus.csv    one row per grouping x subset x property
  physiochemical_stats_pairwise.csv   one row per pairwise comparison

Run: conda run -n tdms python physiochemical_stats.py
"""
from __future__ import annotations
import sqlite3
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from proteoform_physiochemical_props import DB, HERE, load_dataset

# Properties tested (the three drawn in the figure).
PROPERTIES = [
    ("Mass_kDa",         "Average mass (kDa)"),
    ("IsoelectricPoint", "Isoelectric point (pI)"),
    ("GRAVY",            "GRAVY score"),
]

# Cliff's delta magnitude thresholds (Romano et al. 2006).
def delta_magnitude(d: float) -> str:
    a = abs(d)
    if a < 0.147:
        return "negligible"
    if a < 0.330:
        return "small"
    if a < 0.474:
        return "medium"
    return "large"


def eta_squared(H: float, n: int, k: int) -> float:
    """Eta-squared (from H) effect size for Kruskal-Wallis (Tomczak & Tomczak 2014)."""
    return (H - k + 1) / (n - k) if n > k else float("nan")


def cliffs_delta(u: float, nx: int, ny: int) -> float:
    """Cliff's delta from scipy's Mann-Whitney U (U is for the first sample)."""
    return 2.0 * u / (nx * ny) - 1.0


def bh_adjust(pvals: list[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment, returned in the input order."""
    p = np.asarray(pvals, float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def analyze(df: pd.DataFrame, group_col: str, group_order: list[str],
            grouping: str, subset: str,
            omni_rows: list[dict], pair_rows: list[dict]) -> None:
    """Run the omnibus + pairwise battery for every property and append results."""
    for prop, _ in PROPERTIES:
        groups = {g: df.loc[df[group_col] == g, prop].dropna().values
                  for g in group_order}
        groups = {g: v for g, v in groups.items() if len(v) > 0}
        if len(groups) < 2:
            continue
        names = list(groups)
        data = list(groups.values())

        H, p = stats.kruskal(*data)
        n = int(sum(len(v) for v in data))
        k = len(data)
        omni_rows.append(dict(
            Grouping=grouping, Subset=subset, Property=prop, k=k, n=n,
            H=round(float(H), 3), p=float(p),
            eta2=round(eta_squared(H, n, k), 4)))

        recs, p_mwu, p_ks = [], [], []
        for a, b in combinations(names, 2):
            x, y = groups[a], groups[b]
            u, pm = stats.mannwhitneyu(x, y, alternative="two-sided")
            d = cliffs_delta(u, len(x), len(y))
            ks, pk = stats.ks_2samp(x, y)
            recs.append(dict(
                Grouping=grouping, Subset=subset, Property=prop,
                GroupA=a, GroupB=b, nA=len(x), nB=len(y),
                MedianA=round(float(np.median(x)), 3),
                MedianB=round(float(np.median(y)), 3),
                CliffsDelta=round(float(d), 3), DeltaMag=delta_magnitude(d),
                MWU_p=float(pm), KS_stat=round(float(ks), 3), KS_p=float(pk)))
            p_mwu.append(pm)
            p_ks.append(pk)

        for rec, pm_adj, pk_adj in zip(recs, bh_adjust(p_mwu), bh_adjust(p_ks)):
            rec["MWU_p_BH"] = float(pm_adj)
            rec["KS_p_BH"] = float(pk_adj)
            pair_rows.append(rec)


def main() -> None:
    if not DB.exists():
        raise SystemExit(f"Database not found: {DB}")

    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    try:
        meta, props, pf_cond, long = load_dataset(con)
    finally:
        con.close()

    cond_to_bead = dict(zip(meta["Condition"], meta["BeadGroup"]))
    bead_order = list(dict.fromkeys(meta["BeadGroup"]))   # MCW, MagReSyn, Cytiva
    cond_order = meta["Condition"].tolist()

    # bead-type frame: each proteoform counted once per bead type
    bead = long.copy()
    bead["BeadGroup"] = bead["Condition"].map(cond_to_bead)
    bead = bead.drop_duplicates(["ChemicalProteoformId", "BeadGroup"])

    omni_rows: list[dict] = []
    pair_rows: list[dict] = []
    for subset, mask_bead, mask_cond in (
            ("all",         np.ones(len(bead), bool), np.ones(len(long), bool)),
            ("no-histone", ~bead["IsHistone"].to_numpy(), ~long["IsHistone"].to_numpy())):
        analyze(bead[mask_bead], "BeadGroup", bead_order, "bead type", subset,
                omni_rows, pair_rows)
        analyze(long[mask_cond], "Condition", cond_order, "condition", subset,
                omni_rows, pair_rows)

    omni = pd.DataFrame(omni_rows)
    pair = pd.DataFrame(pair_rows)
    omni.to_csv(HERE / "physiochemical_stats_omnibus.csv", index=False)
    pair.to_csv(HERE / "physiochemical_stats_pairwise.csv", index=False)

    # ---- console report -----------------------------------------------------
    pd.set_option("display.width", 200)
    print("Caveat: proteoform property values are fixed per sequence and the per-method")
    print("detection sets overlap, so groups are NOT independent -- these tests are")
    print("descriptive. With n in the hundreds, trust effect sizes over p-values.\n")
    print("eta^2 / Cliff's delta guide: negligible <0.01 / <0.147, "
          "small / medium / large above.\n")

    print("=" * 80)
    print("OMNIBUS  (Kruskal-Wallis across groups)")
    print("=" * 80)
    show = omni.copy()
    show["p"] = show["p"].map(lambda v: f"{v:.2e}")
    print(show.to_string(index=False))

    print("\n" + "=" * 80)
    print("PAIRWISE highlights  (BH-corrected MWU p<0.05 AND |Cliff's delta|>=0.147)")
    print("=" * 80)
    sig = pair[(pair["MWU_p_BH"] < 0.05) & (pair["CliffsDelta"].abs() >= 0.147)]
    if sig.empty:
        print("None: no pair reaches a non-negligible effect size at BH p<0.05.")
    else:
        s = sig[["Grouping", "Subset", "Property", "GroupA", "GroupB",
                 "MedianA", "MedianB", "CliffsDelta", "DeltaMag", "MWU_p_BH"]].copy()
        s["MWU_p_BH"] = s["MWU_p_BH"].map(lambda v: f"{v:.2e}")
        print(s.to_string(index=False))

    print(f"\nLargest |Cliff's delta| overall: {pair['CliffsDelta'].abs().max():.3f}")
    print("\nWrote: physiochemical_stats_omnibus.csv, physiochemical_stats_pairwise.csv")


if __name__ == "__main__":
    main()
