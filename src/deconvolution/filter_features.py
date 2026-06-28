"""
Filter FLASHDeconv feature TSVs to confident masses for MSTopDiff input.

The raw FLASHDeconv feature tables (the MSTopDiff inputs) contain ~50% decoy
features (IsDecoy=1) and many low-confidence masses. FDR control / q-values were
already computed by FLASHDeconv and live in the spectrum-level *_ms1.tsv files
(per-observation Qvalue, with a FeatureIndex linking each observation back to a
feature). This script maps those q-values up to each feature (best = min q-value
over the feature's observations) and writes filtered copies of the feature TSVs.

Filter (defaults):
  * IsDecoy == 0 and MSLevel == 1
  * feature q-value (min over its _ms1 observations) <= FDR_THRESH
  * ChargeCount >= MIN_CHARGE_STATES

NOTE: the charge filter alone does NOT control FDR -- FLASHDeconv decoys have
similar charge counts to targets, so a ChargeCount>=3 set is ~80-99% decoy-like.
The q-value cut does the real confidence work; ChargeCount is a secondary
constraint (and was explicitly requested).

Output: flashdeconv/filtered/<stem>_conf.tsv  (same columns / tab-separated)
"""
import os
import glob
import numpy as np
import pandas as pd

# Repo layout: this script lives in src/deconvolution/, so the project root is
# two levels up; the FLASHDeconv feature tables live under data/flashdeconv/.
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FD = os.path.join(ROOT, "data", "flashdeconv")
OUT = os.path.join(FD, "filtered")

MIN_CHARGE_STATES = 3
FDR_THRESH = 0.05          # feature q-value cutoff (FLASHDeconv anchor = 5%)


def feature_qvalues(ms1_path):
    """Map FeatureIndex -> best (min) Qvalue from the _ms1 spectrum table."""
    m = pd.read_csv(ms1_path, sep="\t",
                    usecols=["FeatureIndex", "Qvalue"], low_memory=False)
    m = m.dropna(subset=["FeatureIndex"])
    m["FeatureIndex"] = m["FeatureIndex"].astype(int)
    m["Qvalue"] = pd.to_numeric(m["Qvalue"], errors="coerce")
    return m.groupby("FeatureIndex")["Qvalue"].min()


def main():
    os.makedirs(OUT, exist_ok=True)
    files = sorted(f for f in glob.glob(os.path.join(FD, "*.tsv"))
                   if not f.endswith("_ms1.tsv")
                   and not os.path.basename(os.path.dirname(f)) == "filtered")
    rows = []
    for path in files:
        stem = os.path.splitext(os.path.basename(path))[0]
        ms1 = os.path.join(FD, f"{stem}_ms1.tsv")
        if not os.path.exists(ms1):
            print(f"  SKIP {stem}: no _ms1.tsv")
            continue
        df = pd.read_csv(path, sep="\t", low_memory=False)
        qmap = feature_qvalues(ms1)
        df["_qvalue"] = df["FeatureIndex"].map(qmap)
        dec = df["IsDecoy"].astype(int)

        tgt = dec == 0
        base = tgt & (df["MSLevel"] == 1)
        keep = base & (df["_qvalue"] <= FDR_THRESH) & (df["ChargeCount"] >= MIN_CHARGE_STATES)

        out = df[keep].drop(columns="_qvalue")
        out.to_csv(os.path.join(OUT, f"{stem}_conf.tsv"), sep="\t", index=False)

        rows.append({
            "file": stem,
            "total": len(df),
            "targets": int(tgt.sum()),
            "tgt_q<=5%": int((base & (df["_qvalue"] <= 0.05)).sum()),
            "tgt_q<=1%": int((base & (df["_qvalue"] <= 0.01)).sum()),
            "tgt_chg>=3": int((base & (df["ChargeCount"] >= 3)).sum()),
            "KEPT (q<=5% & chg>=3)": int(keep.sum()),
            "no_qval_map": int((tgt & df["_qvalue"].isna()).sum()),
        })
    summ = pd.DataFrame(rows)
    pd.options.display.width = 200
    print(f"Wrote filtered TSV(s) to {OUT}")
    print(f"Filter: IsDecoy==0, MSLevel==1, feature q-value<={FDR_THRESH:.0%} "
          f"(from _ms1), ChargeCount>={MIN_CHARGE_STATES}\n")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
