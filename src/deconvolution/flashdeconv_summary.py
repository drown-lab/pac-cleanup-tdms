#!/usr/bin/env python3
"""Summarize FLASHDeconv FDR output across a folder of spectrum-level TSV files.

For the SP3 top-down method comparison. Each input is a FLASHDeconv
spectrum-level TSV produced with ``-FD:report_FDR`` (e.g. ``*_ms1.tsv``). The
file contains one row per deconvolved mass per MS1 spectrum, with target and
decoy masses interleaved and a per-mass Qscore and q-value.

For every file this script:
  * separates target masses from decoys,
  * counts confident target masses at fixed deconvolution-FDR thresholds
    (the "yield at q <= 0.05" number that anchors the method comparison),
  * summarizes the target Qscore distribution,
and writes a tidy summary table plus distribution / yield figures.

Note on counting: rows are spectrum-level mass observations, so "yield" is the
number of confident masses summed over all MS1 spectra in the file - the same
"positive masses" notion used in the deconvolution-FDR paper. If you want
LC-feature-level counts instead, point this at the feature TSV (the ``-out``
file) rather than the ``-out_spec1`` file.

Outputs: figures to results/figures/, the summary table to results/tables/.

Usage
-----
    # defaults to reading data/flashdeconv
    python src/deconvolution/flashdeconv_summary.py
    python src/deconvolution/flashdeconv_summary.py path/to/tsvs \
        --glob "*_ms1.tsv" --thresholds 0.01 0.05 0.10 --per-file-plots
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Repo layout: this script lives in src/deconvolution/, so the project root is
# two levels up. Inputs default to data/flashdeconv/; figures go to
# results/figures/ and the summary table to results/tables/.
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = ROOT / "data" / "flashdeconv"
FIGURES_DIR = ROOT / "results" / "figures"
TABLES_DIR = ROOT / "results" / "tables"

# Canonical name -> candidate raw column names (matched after normalization).
# Column spelling varies across FLASHDeconv builds, so we resolve by candidates
# rather than hard-coding, and report what was found.
COLUMN_CANDIDATES = {
    "mass": ["monoisotopicmass", "monomass", "monoisotopicmassda", "mass"],
    "qscore": ["qscore"],
    "qvalue": ["qvalue", "qval"],
    "decoy": ["decoy", "decoyflag", "isdecoy", "targetdecoytype", "decoytype"],
}

# String tokens treated as "target" when the decoy column is non-numeric.
TARGET_STRINGS = {"target", "0", "false", "f", "no", "n"}


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def detect_columns(df: pd.DataFrame) -> dict[str, str]:
    """Map canonical roles -> actual column names present in ``df``."""
    norm_map = {_norm(c): c for c in df.columns}
    resolved: dict[str, str] = {}
    for canon, candidates in COLUMN_CANDIDATES.items():
        for cand in candidates:
            if cand in norm_map:
                resolved[canon] = norm_map[cand]
                break
    return resolved


def target_mask(decoy_series: pd.Series, target_value: int = 0) -> pd.Series:
    """Boolean mask of target (non-decoy) rows, robust to numeric/string flags."""
    numeric = pd.to_numeric(decoy_series, errors="coerce")
    if numeric.notna().all():
        return numeric == target_value
    return decoy_series.astype(str).str.strip().str.lower().isin(TARGET_STRINGS)


def load_tsv(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    cols = detect_columns(df)
    missing = [k for k in ("qscore", "qvalue", "decoy") if k not in cols]
    if missing:
        raise ValueError(
            f"{path.name}: could not find column(s) for {missing}. "
            f"Headers present: {list(df.columns)}. "
            "Was FLASHDeconv run with -FD:report_FDR? Adjust COLUMN_CANDIDATES "
            "if your build uses different names."
        )
    return df, cols


def summarize_file(
    path: Path, thresholds: list[float], target_value: int, verbose: bool
) -> dict:
    df, cols = load_tsv(path)
    is_target = target_mask(df[cols["decoy"]], target_value)
    targets = df.loc[is_target]
    decoys = df.loc[~is_target]

    if verbose:
        uniq = sorted(df[cols["decoy"]].dropna().unique().tolist())[:10]
        print(
            f"  {path.name}: columns -> "
            f"mass={cols.get('mass', 'n/a')}, qscore={cols['qscore']}, "
            f"qvalue={cols['qvalue']}, decoy={cols['decoy']} "
            f"(values seen: {uniq})"
        )

    qval = pd.to_numeric(targets[cols["qvalue"]], errors="coerce")
    qscore = pd.to_numeric(targets[cols["qscore"]], errors="coerce")

    row = {
        "file": path.stem,
        "n_target_masses": int(len(targets)),
        "n_decoy_masses": int(len(decoys)),
        "target_qscore_median": float(qscore.median()) if len(qscore) else np.nan,
        "target_qscore_mean": float(qscore.mean()) if len(qscore) else np.nan,
    }
    for t in thresholds:
        row[f"yield_q<={t:g}"] = int((qval <= t).sum())
    return row


def yield_curve(path: Path, target_value: int, grid: np.ndarray) -> np.ndarray:
    """Cumulative confident target-mass count over a grid of FDR thresholds."""
    df, cols = load_tsv(path)
    is_target = target_mask(df[cols["decoy"]], target_value)
    qval = pd.to_numeric(df.loc[is_target, cols["qvalue"]], errors="coerce").to_numpy()
    qval = qval[~np.isnan(qval)]
    return np.array([(qval <= t).sum() for t in grid])


def plot_qscore_distributions(
    files: list[Path], target_value: int, out_path: Path
) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.linspace(0, 1, 41)
    for path in files:
        df, cols = load_tsv(path)
        is_target = target_mask(df[cols["decoy"]], target_value)
        q = pd.to_numeric(df.loc[is_target, cols["qscore"]], errors="coerce").dropna()
        ax.hist(q, bins=bins, histtype="step", linewidth=1.5, label=path.stem)
    ax.set_xlabel("Qscore (target masses)")
    ax.set_ylabel("Mass count")
    ax.set_title("Deconvolution score distribution by file")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_yield_curves(
    files: list[Path], target_value: int, out_path: Path
) -> None:
    grid = np.linspace(0.0, 0.20, 41)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for path in files:
        counts = yield_curve(path, target_value, grid)
        ax.plot(grid * 100, counts, linewidth=1.5, label=path.stem)
    ax.axvline(5, color="0.5", linestyle="--", linewidth=1)
    ax.set_xlabel("Deconvolution FDR threshold (%)")
    ax.set_ylabel("Confident target masses")
    ax.set_title("Confident-mass yield vs FDR threshold")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_target_decoy(path: Path, target_value: int, out_path: Path) -> None:
    """Per-file target vs decoy Qscore step plot (the paper's key diagnostic)."""
    df, cols = load_tsv(path)
    is_target = target_mask(df[cols["decoy"]], target_value)
    bins = np.linspace(0, 1, 41)
    fig, ax = plt.subplots(figsize=(6, 4))
    for mask, label, color in [
        (is_target, "target", "#1f77b4"),
        (~is_target, "decoy", "#d62728"),
    ]:
        q = pd.to_numeric(df.loc[mask, cols["qscore"]], errors="coerce").dropna()
        ax.hist(q, bins=bins, histtype="step", linewidth=1.6, label=label, color=color)
    ax.set_xlabel("Qscore")
    ax.set_ylabel("Mass count")
    ax.set_title(path.stem)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_dir", type=Path, nargs="?", default=DEFAULT_INPUT_DIR,
                   help="Folder containing FLASHDeconv spectrum TSVs "
                        "(default: data/flashdeconv).")
    p.add_argument("--glob", default="*_ms1.tsv", help="Filename pattern (default: %(default)s).")
    p.add_argument("--thresholds", type=float, nargs="+", default=[0.01, 0.05, 0.10], help="FDR thresholds for the yield table.")
    p.add_argument("--target-value", type=int, default=0, help="Decoy-column value denoting a target mass (default: 0).")
    p.add_argument("--per-file-plots", action="store_true", help="Also write a target/decoy step plot per file.")
    p.add_argument("--quiet", action="store_true", help="Suppress per-file column diagnostics.")
    args = p.parse_args(argv)

    files = sorted(args.input_dir.glob(args.glob))
    if not files:
        print(f"No files matching {args.glob!r} in {args.input_dir}", file=sys.stderr)
        return 1

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Found {len(files)} file(s). Summarizing...")
    rows = []
    for path in files:
        try:
            rows.append(
                summarize_file(path, args.thresholds, args.target_value, verbose=not args.quiet)
            )
        except ValueError as e:
            print(f"  SKIPPED {e}", file=sys.stderr)

    if not rows:
        print("No files could be summarized.", file=sys.stderr)
        return 1

    summary = pd.DataFrame(rows).sort_values("file").reset_index(drop=True)
    csv_path = TABLES_DIR / "deconvolution_fdr_summary.csv"
    summary.to_csv(csv_path, index=False)
    print(f"\nWrote {csv_path}")
    print(summary.to_string(index=False))

    ok_files = [f for f in files if f.stem in set(summary["file"])]

    plot_qscore_distributions(ok_files, args.target_value, FIGURES_DIR / "qscore_distributions.png")
    plot_yield_curves(ok_files, args.target_value, FIGURES_DIR / "yield_vs_fdr.png")
    print(f"Wrote {FIGURES_DIR / 'qscore_distributions.png'}")
    print(f"Wrote {FIGURES_DIR / 'yield_vs_fdr.png'}")

    if args.per_file_plots:
        pf_dir = FIGURES_DIR / "flashdeconv_per_file"
        pf_dir.mkdir(exist_ok=True)
        for path in ok_files:
            plot_target_decoy(path, args.target_value, pf_dir / f"{path.stem}_target_decoy.png")
        print(f"Wrote per-file target/decoy plots to {pf_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())