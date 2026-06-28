"""
Evaluate UNANNOTATED delta-mass peaks in the MSTopDiff histograms: peaks that
rise above local background but do not correspond to a known modification
(single, pairwise combination, or the integer-Da / isotope comb).

Approach
--------
1. Detect peaks on the native 0.01 Da intensity x count trace with find_peaks,
   prominence set relative to a rolling-median local background.
2. Classify each peak centroid as explained if within TOL of:
     - a single known modification,
     - a pairwise sum of two known modifications (<=150 Da), or
     - an integer Da value (isotope / unit-mass comb; reported separately).
   Otherwise it is UNANNOTATED.
3. Report per condition: n peaks, summed intensity x count, and the
   intensity-weighted fraction unexplained (robust to threshold choice).
4. Sweep the prominence multiple to show threshold sensitivity.

Outputs:
  mstopdiff_unannotated_summary.csv  - per condition, at the chosen threshold
  mstopdiff_unannotated_peaks.csv    - every unannotated peak (chosen threshold)
  mstopdiff_unannotated_sensitivity.csv - sweep over prominence multiple
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 13

import mstopdiff_config as cfg

FIGURES_DIR = cfg.FIGURES_DIR  # results/figures
TABLES_DIR = cfg.TABLES_DIR    # results/tables
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# Known single modification masses (Da). Same curated list as the figure script.
MOD_MASSES = [
    0.984, 14.016, 15.995, 18.011, 21.982, 27.995, 28.031, 31.990, 37.956,
    42.011, 42.047, 43.006, 43.990, 44.026, 45.988, 47.985, 57.021, 58.005,
    63.980, 71.037, 75.998, 79.966, 79.957, 86.000, 100.016, 114.043, 119.004,
    128.059, 128.095,
]

MAX_MASS = 90.0     # restrict the whole analysis to 0-90 Da delta mass
TOL = 0.02          # Da, annotation match tolerance
PEAK_DIST_DA = 0.04 # min spacing between detected peaks (Da)
PEAK_HALF_DA = 0.02 # +/- window summed as a peak's intensity (Da)
REL_THRESH = 0.10   # a peak must reach >=10% of the base (tallest) peak
REL_SWEEP = (0.03, 0.05, 0.10, 0.15, 0.20)  # sensitivity on the threshold
ISOTOPE = 1.00235   # averagine isotope spacing (Da): off-by-one deconv error
SAT_K = (1,)        # off-by-one satellites only (+/- 1 x ISOTOPE); the
                    # dominant deconvolution error. Extend to (1, 2) to also
                    # capture rarer off-by-two (raises chance coverage ~14->24%).


def load(path):
    df = pd.read_csv(path)
    df["ixc"] = (df["intensity_lower_mass x count"].abs()
                 + df["intensity_higher_mass x count"].abs())
    return df


def explained_singles_and_combos(masses, tol, max_mass=MAX_MASS):
    """Sorted array of explained masses (<= max_mass): singles + pairwise sums."""
    vals = set()
    for m in masses:
        if m <= max_mass + tol:
            vals.add(round(m, 3))
    for i, a in enumerate(masses):
        for b in masses[i:]:
            s = a + b
            if s <= max_mass + tol:
                vals.add(round(s, 3))
    return np.array(sorted(vals))


def satellite_masses(expl_sorted, max_mass=MAX_MASS):
    """Off-by-one(/two) deconvolution satellites of every explained mass."""
    vals = set()
    for e in expl_sorted:
        for k in SAT_K:
            for s in (e + k * ISOTOPE, e - k * ISOTOPE):
                if 0 < s <= max_mass:
                    vals.add(round(s, 3))
    return np.array(sorted(vals))


def _near(mass, sorted_arr, tol):
    if len(sorted_arr) == 0:
        return False
    idx = np.searchsorted(sorted_arr, mass)
    for j in (idx - 1, idx):
        if 0 <= j < len(sorted_arr) and abs(sorted_arr[j] - mass) <= tol:
            return True
    return False


def classify(mass, expl_sorted, sat_sorted, tol):
    """Return: 'known', 'satellite', 'integer', or 'unannotated'.

    Order matters: an exact known match wins over a satellite match, which
    wins over the integer/isotope (zero-shift) comb.
    """
    if _near(mass, expl_sorted, tol):
        return "known"
    if _near(mass, sat_sorted, tol):           # off-by-one satellite of a mod
        return "satellite"
    if abs(mass - round(mass)) <= tol:         # integer / zero-shift isotope comb
        return "integer"
    return "unannotated"


def detect(df, rel_thresh):
    """A peak is a local maximum reaching >= rel_thresh of the base peak."""
    x = df["bin"].to_numpy()
    y = df["ixc"].to_numpy().astype(float)
    step = x[1] - x[0]
    mask = x <= MAX_MASS                  # restrict to 0-MAX_MASS Da
    x, y = x[mask], y[mask]
    base = y.max()                       # base (tallest) peak in range
    idx, _ = find_peaks(y, height=rel_thresh * base,
                        distance=max(1, int(round(PEAK_DIST_DA / step))))
    half = int(round(PEAK_HALF_DA / step))
    peaks = []
    for i in idx:
        lo, hi = max(0, i - half), min(len(y), i + half + 1)
        peaks.append({"mass": float(x[i]), "ixc": float(y[lo:hi].sum()),
                      "rel": float(y[i] / base)})
    return pd.DataFrame(peaks)


def summarize(df, expl_sorted, sat_sorted, rel_thresh):
    pk = detect(df, rel_thresh)
    if pk.empty:
        return None, pk
    pk["class"] = pk["mass"].apply(
        lambda m: classify(m, expl_sorted, sat_sorted, TOL))
    g = pk.groupby("class")
    n = g.size()
    ixc = g["ixc"].sum()
    total_ixc = pk["ixc"].sum()
    una = pk[pk["class"] == "unannotated"]
    row = {
        "n_peaks": len(pk),
        "n_known": int(n.get("known", 0)),
        "n_satellite": int(n.get("satellite", 0)),
        "n_integer": int(n.get("integer", 0)),
        "n_unannotated": int(n.get("unannotated", 0)),
        "ixc_total": total_ixc,
        "ixc_known": float(ixc.get("known", 0.0)),
        "ixc_satellite": float(ixc.get("satellite", 0.0)),
        "ixc_integer": float(ixc.get("integer", 0.0)),
        "ixc_unannotated": float(ixc.get("unannotated", 0.0)),
        "frac_unannot_by_count": round(len(una) / len(pk), 3),
        "frac_unannot_by_ixc": round(una["ixc"].sum() / total_ixc, 3),
    }
    return row, pk


def main():
    dsets = cfg.discover()
    if not dsets:
        raise SystemExit("No *_conf_mstopdiff.csv found in flashdeconv/mstopdiff/")
    data = {d["label"]: load(d["csv"]) for d in dsets}
    expl = explained_singles_and_combos(MOD_MASSES, TOL)
    sat = satellite_masses(expl)
    # chance coverage: fraction of 0-MAX_MASS Da within tol of each window set
    cov_known = len(expl) * 2 * TOL / MAX_MASS
    cov_sat = len(sat) * 2 * TOL / MAX_MASS
    n_in_range = sum(m <= MAX_MASS + TOL for m in MOD_MASSES)
    print(f"Analysis range: 0-{MAX_MASS:g} Da")
    print(f"Explained mass list: {n_in_range} singles -> "
          f"{len(expl)} singles+combos within tol {TOL} Da")
    print(f"Off-by-one satellites (+/-{SAT_K} x {ISOTOPE} Da): {len(sat)} masses")
    print(f"Chance coverage of 0-{MAX_MASS:g} Da: known {cov_known:.1%}, "
          f"+satellites {cov_known + cov_sat:.1%}\n")

    # ---- main summary at chosen threshold ----
    summ_rows, all_unannot = [], []
    for cond, df in data.items():
        row, pk = summarize(df, expl, sat, REL_THRESH)
        if row is None:
            print(f"  WARNING: no peaks >= {REL_THRESH:.0%} of base peak in "
                  f"{cond}; skipped")
            continue
        row = {"condition": cond, **row}
        summ_rows.append(row)
        una = pk[pk["class"] == "unannotated"].copy()
        una["condition"] = cond
        all_unannot.append(una)
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(os.path.join(TABLES_DIR, "mstopdiff_unannotated_summary.csv"),
                index=False)
    unannot = pd.concat(all_unannot, ignore_index=True)[
        ["condition", "mass", "rel", "ixc"]].sort_values(
            ["condition", "ixc"], ascending=[True, False])
    unannot.to_csv(os.path.join(TABLES_DIR, "mstopdiff_unannotated_peaks.csv"),
                   index=False)

    # ---- sensitivity sweep over the relative-abundance threshold ----
    sweep = []
    for rt in REL_SWEEP:
        for cond, df in data.items():
            row, _ = summarize(df, expl, sat, rt)
            if row:
                sweep.append({"rel_thresh": rt, "condition": cond,
                              "n_peaks": row["n_peaks"],
                              "n_unannotated": row["n_unannotated"],
                              "frac_unannot_by_count": row["frac_unannot_by_count"],
                              "frac_unannot_by_ixc": row["frac_unannot_by_ixc"]})
    sweep = pd.DataFrame(sweep)
    sweep.to_csv(os.path.join(TABLES_DIR, "mstopdiff_unannotated_sensitivity.csv"),
                 index=False)

    # ---- composition figure: intensity x count split by class ----
    cats = [("ixc_known", "known single / combo", "#4C72B0"),
            ("ixc_satellite", "off-by-one satellite of mod", "#55A868"),
            ("ixc_integer", "integer / isotope comb", "#999999"),
            ("ixc_unannotated", "unannotated", "#C44E52")]
    conds = list(data)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(max(11, 1.1 * len(conds) + 4), 5.5))
    xs = np.arange(len(conds))
    # absolute
    bottom = np.zeros(len(conds))
    for col, lab, c in cats:
        vals = summ.set_index("condition").loc[conds, col].to_numpy()
        axA.bar(xs, vals, bottom=bottom, color=c, label=lab)
        bottom += vals
    axA.set_xticks(xs); axA.set_xticklabels(conds, fontsize=8, rotation=45, ha="right")
    axA.set_ylabel("summed intensity x count of detected peaks")
    axA.set_title("Absolute")
    # fraction
    bottom = np.zeros(len(conds))
    tot = summ.set_index("condition").loc[conds, "ixc_total"].to_numpy()
    for col, lab, c in cats:
        vals = summ.set_index("condition").loc[conds, col].to_numpy() / tot
        axB.bar(xs, vals, bottom=bottom, color=c, label=lab)
        bottom += vals
    axB.set_xticks(xs); axB.set_xticklabels(conds, fontsize=9)
    axB.set_ylabel("fraction of detected-peak intensity x count")
    axB.set_ylim(0, 1)
    axB.set_title("Composition")
    axB.legend(fontsize=8, loc="lower right")
    fig.suptitle(f"Detected delta-mass peak intensity by class "
                 f"(peaks >= {REL_THRESH:.0%} of base peak)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGURES_DIR, f"fig_mstopdiff_unannotated.{ext}"),
                    dpi=200, bbox_inches="tight")

    # ---- console ----
    pd.options.display.width = 160
    pd.options.display.float_format = lambda v: f"{v:,.3g}"
    print(f"=== Summary: peaks >= {REL_THRESH:.0%} of base peak ===")
    print(summ[["condition", "n_peaks", "n_known", "n_satellite", "n_integer",
                "n_unannotated", "frac_unannot_by_count",
                "frac_unannot_by_ixc"]].to_string(index=False))
    print("\nIntensity x count totals by class:")
    print(summ[["condition", "ixc_known", "ixc_satellite", "ixc_integer",
                "ixc_unannotated", "ixc_total"]].to_string(index=False))
    print("\nIntensity x count fraction by class:")
    fr = summ.set_index("condition")
    for c in ("ixc_known", "ixc_satellite", "ixc_integer", "ixc_unannotated"):
        fr[c.replace("ixc_", "f_")] = fr[c] / fr["ixc_total"]
    print(fr[["f_known", "f_satellite", "f_integer", "f_unannotated"]].round(3)
          .to_string())
    print("\n=== Threshold sensitivity: total n_peaks ===")
    print(sweep.pivot(index="rel_thresh", columns="condition",
                      values="n_peaks").to_string())
    print("\n=== Threshold sensitivity: frac unannotated by ixc ===")
    print(sweep.pivot(index="rel_thresh", columns="condition",
                      values="frac_unannot_by_ixc").to_string())
    print("\nTop unannotated peaks per condition:")
    for cond in conds:
        top = unannot[unannot["condition"] == cond].head(8)
        print(f"\n  {cond}:")
        print(top[["mass", "rel", "ixc"]].to_string(index=False))


if __name__ == "__main__":
    main()
