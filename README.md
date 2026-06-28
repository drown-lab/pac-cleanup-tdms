# PAC for TDP

## Overview

This repository contains the custom Python used to support a study from the
Drown lab on bead-assisted **p**rotein **a**ggregation **c**apture (PAC) as a
sample-cleanup approach for **t**op-**d**own **p**roteomics (TDP). The code
accompanies the associated publication and is provided as-is to make the
analysis reproducible.

## Repository layout

```
.
├── config/
│   └── experimental_design.csv     # sample → cleanup/resuspension/include map
├── data/                           # all inputs — git-ignored; see data/README.md
│   ├── raw/  mzML/  flashdeconv/
│   └── *.tdReport
├── src/
│   ├── deconvolution/              # FLASHDeconv mass deconvolution
│   ├── modifications/              # MSTopDiff Δmass modification analysis
│   └── identification/             # ProSightPD .tdReport identifications
└── results/
    ├── figures/                    # publication figures (.pdf + .png), tracked
    └── tables/                     # summary/stat tables (bulky dumps git-ignored)
```

The code is grouped into three independent analysis pipelines under `src/`:

**`src/deconvolution/`** — FLASHDeconv mass deconvolution

| Script | Purpose |
| ------ | ------- |
| `run_FLASHDeconv.ps1`, `run_FLASHDeconv_parallel.ps1` | Run FLASHDeconv over a folder of mzML (serial / one-thread-per-file parallel). |
| `flashdeconv_summary.py` | Summarize FLASHDeconv deconvolution-FDR output (confident-mass yield, Qscore distributions) across spectrum-level `*_ms1.tsv` files. |
| `filter_features.py` | Filter FLASHDeconv feature TSVs to confident masses (drop decoys, q-value ≤ 5%, ≥3 charge states) for re-running MSTopDiff. |

**`src/modifications/`** — MSTopDiff Δmass modification analysis

| Script | Purpose |
| ------ | ------- |
| `mstopdiff_compare.py` | Build the MSTopDiff Δmass modification figures (per-condition intensity×count histograms, relative-abundance bars and heatmap). |
| `mstopdiff_unannotated.py` | Classify detected Δmass peaks as known / off-by-one satellite / isotope / unannotated and quantify the unexplained fraction. |
| `mstopdiff_config.py` | Shared dataset registry (auto-discovers MSTopDiff exports, labels/colours them from `config/experimental_design.csv`); imported by the two scripts above. |

**`src/identification/`** — ProSightPD `.tdReport` identifications

| Script | Purpose |
| ------ | ------- |
| `proteoform_physiochemical_props.py` | Calculate biochemical properties (GRAVY, pI, aromaticity, instability) for FDR-confident proteoforms pulled directly from a ProSightPD `.tdReport`, labeled by cleanup condition. |
| `physiochemical_stats.py` | Test whether the physiochemical property distributions (mass, pI, GRAVY) vary by cleanup method — by bead type and by condition — with effect sizes (shares the data pipeline of `proteoform_physiochemical_props.py`). |
| `shared_proteoforms.py` | Compare identification scores of proteoforms shared across sample-cleanup conditions from a ProSightPD `.tdReport`. |

Each script reads inputs from `data/` and `config/` and writes figures/tables to
`results/`, resolved relative to the repository root, so they can be run from
anywhere. Inputs (raw data and the `.tdReport`) are not included — see
[`data/README.md`](data/README.md).

---

## Requirements

* Python 3.12 or later

Install the dependencies with conda:

```bash
conda env create -f environment.yml      # loose pins, latest compatible
conda activate tdms
```

For an exact reproduction of the published analysis, recreate the pinned
environment instead:

```bash
conda env create -f environment.lock.yml
```

---

## Reproducing the analysis

### Data flow

```
 data/raw/  (Thermo .raw)                 data/<...>.tdReport  (ProSightPD search)
     │  MSConvert (ProteoWizard)                    │
     ▼                                             ▼
 data/mzML/                              ┌──────────────────────────────────────┐
     │  FLASHDeconv                      │ src/identification/                    │
     │  (src/deconvolution/run_*.ps1)    │   shared_proteoforms.py                │
     ▼                                   │   proteoform_physiochemical_props.py   │
 data/flashdeconv/                       │   physiochemical_stats.py              │
   ├─ <stem>.tsv      (features)         └──────────────────┬─────────────────────┘
   └─ <stem>_ms1.tsv  (spectra, q-vals)                     ▼
     │  filter_features.py                              results/
     ▼
 data/flashdeconv/filtered/<stem>_conf.tsv
     │  MSTopDiff (GUI)
     ▼
 data/flashdeconv/mstopdiff/<stem>_conf_mstopdiff.csv
     │  mstopdiff_compare.py / mstopdiff_unannotated.py
     ▼
 results/   (figures + tables)
```

### Steps

Steps marked **(external)** are run outside this repo; place their outputs under
`data/` as shown (see [`data/README.md`](data/README.md)). The Python steps read
`data/`/`config/` and write `results/`.

1. **(external)** Acquire the `.raw` files → `data/raw/`.
2. **(external)** Convert `.raw` → `.mzML` with MSConvert (ProteoWizard) → `data/mzML/`.
3. **(external)** Deconvolve with FLASHDeconv → `data/flashdeconv/`
   (`src/deconvolution/run_FLASHDeconv.ps1` or `run_FLASHDeconv_parallel.ps1`;
   edit the in-script paths for your machine).
4. `flashdeconv_summary.py` — deconvolution-FDR QC (yield, Qscore).
5. `filter_features.py` — write confident features → `data/flashdeconv/filtered/`.
6. **(external)** Run MSTopDiff (GUI) on the filtered features with the parameters
   in the [MSTopDiff section](#mstopdiff-modification-analysis) →
   `data/flashdeconv/mstopdiff/`.
7. `mstopdiff_compare.py` and `mstopdiff_unannotated.py` — Δmass figures/tables.
8. **(external)** ProSightPD search of the `.raw` files → the `.tdReport` in `data/`.
9. `shared_proteoforms.py`, `proteoform_physiochemical_props.py`,
   `physiochemical_stats.py` — identification + physiochemical analysis.

The Python steps (4, 5, 7, 9) can be run individually, or together once `data/`
is populated:

```bash
make all          # runs every Python step (override the interpreter with PY=...)
make identification   # just the .tdReport analyses; see the Makefile for targets
```

`make` invokes `conda run -n tdms python …`; without `make`, run the commands
under each script's **Running** heading below.

### Software versions

The published results were produced with:

| Software | Version |
| -------- | ------- |
| Python | 3.12.13 |
| NumPy | 2.5.0 |
| pandas | 3.0.3 |
| SciPy | 1.18.0 |
| matplotlib | 3.11.0 |
| matplotlib-venn | 1.1.2 |
| seaborn | 0.13.2 |
| Biopython | 1.87 |
| FLASHDeconv (OpenMS) | 3.5.0 |
| MSConvert (ProteoWizard) | 3.0.25326-3a190dc |
| MSTopDiff | 1.1.0 |
| ProSightPD | 4.2 |

Python package versions are pinned exactly in `environment.lock.yml`.

---

## Proteoform physiochemical properties (`proteoform_physiochemical_props.py`)

Pulls FDR-confident proteoform sequences **directly from a ProSightPD
`.tdReport`** (a SQLite database) and calculates several biochemical properties
for each, across **all** sample-cleanup conditions grouped by bead type. This
replaces the earlier workflow of reading sequences from a hand-exported TDreport
hit-report CSV.

### Calculated properties

* **Average mass (kDa):** proteoform mass. The DB's `AverageMass` column is only
  ~40 % populated (placeholder `18.0153` Da otherwise), so the fully-populated
  `MonoisotopicMass` is used instead; the two differ by ~0.06 % at these sizes
  (immaterial on a kDa axis).
* **GRAVY (Grand Average of Hydropathicity):** overall hydrophobicity/hydrophilicity of a protein.
* **Isoelectric Point (pI):** the pH at which the protein carries no net charge.
* **Aromaticity** and **Instability Index** are also computed into the CSV/summary
  (not plotted).

GRAVY / pI / aromaticity / instability depend only on the bare amino-acid
sequence; modifications (stored separately in `ModificationHash`) are ignored.
Selenocysteine (`U`) is substituted with `C` for those calculations only —
Biopython's hydropathy and instability tables have no `U` — and a `ContainsU`
flag marks any affected proteoforms (the original sequence is preserved).

Histones dominate this dataset (roughly half the confident proteoforms), so
every summary and figure is produced **twice**: once for all confident
proteoforms and once **excluding histones**. A proteoform is classed as a
histone (`IsHistone` column) when a comma-delimited entry of its description
starts with "Histone" — which excludes "Non-histone chromosomal protein HMG-…".

### Inputs

Uses the shared inputs — the `.tdReport` in `data/` and
`config/experimental_design.csv` (see [Repository layout](#repository-layout) and
[`data/README.md`](data/README.md)). Identification uses the same
`ChemicalProteoformId` identity and PrSM-level q ≤ 0.01 confidence rule as
`shared_proteoforms.py`.

Conditions and their left-to-right order, colours and bead-type grouping are set
by the `CONDITIONS` list (all nine included files: MCW, then the MagReSyn and
Cytiva resuspension series).

### Outputs

| File | Description |
| ---- | ----------- |
| `proteoform_physiochemical_props.csv` | One row per unique proteoform: properties + `IsHistone` + per-condition presence |
| `proteoform_physiochemical_props_long.csv` | Tidy table, one row per proteoform × condition |
| `physiochemical_summary.csv` | Per-condition × per-property summary statistics (all proteoforms) |
| `physiochemical_summary_nohistone.csv` | Per-condition × per-property summary statistics, histones excluded |
| `fig_physiochemical_props.{pdf,png}` | Panel A: average mass, pI and GRAVY stacked, one violin per condition grouped by bead type (all proteoforms) |
| `fig_physiochemical_props_nohistone.{pdf,png}` | Same stacked panel, histones excluded |

The figures use **violin plots** (kernel-density) rather than box plots because
the pI distribution is strongly bimodal, which a box plot hides. The script also
prints per-condition summary statistics to the console.

### Running

```bash
conda run -n tdms python src/identification/proteoform_physiochemical_props.py
```

---

## Testing properties by cleanup method (`physiochemical_stats.py`)

Tests whether the property distributions (average mass, pI, GRAVY) **vary by
cleanup method**, reusing the exact FDR-confident dataset assembled by
`proteoform_physiochemical_props.py` (imported via its shared `load_dataset()`).

Two groupings are tested — **bead type** (MCW vs MagReSyn vs Cytiva, resuspensions
pooled) and **all nine conditions** — each on all proteoforms and with histones
excluded, for each property.

**Key caveat:** a proteoform's property value is a fixed attribute of its
sequence; methods differ only in *which* proteoforms they detect, and the
detected sets overlap heavily. The groups are therefore not independent and the
tests are **descriptive**. With hundreds of proteoforms per group, p-values are
tiny for even trivial differences, so the **effect sizes** (η² for the omnibus,
Cliff's δ for pairwise) are what matter.

Battery (scipy only): Kruskal–Wallis omnibus + η²; pairwise Mann–Whitney U +
Cliff's δ; pairwise two-sample Kolmogorov–Smirnov (distribution shape).
Pairwise p-values are Benjamini–Hochberg corrected within each
grouping × subset × property family.

### Outputs

| File | Description |
| ---- | ----------- |
| `physiochemical_stats_omnibus.csv` | One row per grouping × subset × property: K–W H, p, η² |
| `physiochemical_stats_pairwise.csv` | One row per pairwise comparison: Cliff's δ, MWU & KS p (raw + BH) |

### Running

```bash
conda run -n tdms python src/identification/physiochemical_stats.py
```

---

## Shared proteoforms across conditions (`shared_proteoforms.py`)

Compares the `kelleher_negLog_pScore` of proteoforms that are identified across a
chosen set of sample-cleanup conditions, using a ProSightPD `.tdReport` (a SQLite
database) as the source of identifications.

The study background is a single HeLa lysate prepared with different cleanup
approaches. The script compares three conditions by default:

| Condition          | Cleanup            | Resuspension      |
| ------------------ | ------------------ | ----------------- |
| MCW                | MCW                | 5% ACN, 0.1% FA   |
| Cytiva 0.5% TFA    | Cytiva Carboxyl    | 0.5% TFA          |
| MagReSyn 0.5% TFA  | MagReSyn Hydroxyl  | 0.5% TFA          |

Key analysis decisions (see the module docstring for details):

* **Proteoform identity** = `ChemicalProteoformId` (sequence + modifications + mass).
* **Per-condition score** = the best (maximum) `kelleher_negLog_pScore` over that
  proteoform's FDR-passing hits (PrSM-level q-value ≤ 0.01) in that condition's files.
* **"Shared"** = identified in *every* compared condition.

### Inputs

Uses the shared inputs — the `.tdReport` in `data/` and
`config/experimental_design.csv` (see [Repository layout](#repository-layout) and
[`data/README.md`](data/README.md)). The design has columns
`Name, Cleanup, Resuspension, Include`; rows with `Include == FALSE` (e.g. files
acquired with a different method) are dropped. To compare different conditions,
edit the `CONDITIONS_OF_INTEREST` list in the script.

### Outputs

| File | Description |
| ---- | ----------- |
| `shared_proteoforms_wide.csv` | Proteoform × condition matrix of max negLog P-scores |
| `shared_proteoforms_long.csv` | Tidy table, one row per proteoform × condition |
| `per_sample_score_summary.csv` | Per-condition summary statistics |
| `identification_membership.csv` | Presence/absence of every proteoform per condition |
| `fig_score_distribution.{pdf,png}` | Box/strip plot of scores per condition |
| `fig_score_heatmap.{pdf,png}` | Clustered heatmap (proteoform × condition) |
| `fig_score_paired.{pdf,png}` | Paired lines: same proteoform across conditions |
| `fig_venn_identifications.{pdf,png}` | Venn diagram of identified proteoform sets |

The script also prints a summary to the console, including a **Friedman test**
across conditions and **pairwise Wilcoxon signed-rank** tests on the shared
proteoforms.

### Running

```bash
conda run -n tdms python src/identification/shared_proteoforms.py
```

---

## MSTopDiff modification analysis

This set of scripts evaluates whether each sample-cleanup approach **introduces
artifactual modifications or suppresses PTMs**, using
[MSTopDiff](https://github.com/PhilippKaulich/MSTopDiff) Δmass histograms of
the FLASHDeconv masses. The approach mirrors Kaulich *et al.* 2024
(*Nat. Methods* **21**, 2397–2407), Figure 5d: per-condition intensity×count
histograms of all pairwise proteoform mass differences, with peaks annotated
against common modifications.

MSTopDiff (v1.1.0) was run on the confident-mass feature tables with these
parameters:

| Parameter | Value |
| --------- | ----- |
| Mass feature filter | 0 – 100 kDa |
| Retention time range | 0 – 60 min |
| Delta mass range | 0 – 150 Da |
| Retention time window | 2 min |
| Maximum charge difference | 2 |
| Bin size | 0.01 Da |

### Confidence filtering (`filter_features.py`)

The raw FLASHDeconv **feature** TSVs (the MSTopDiff inputs) contain ~50 % decoy
features (`IsDecoy=1`) and many low-confidence masses — only ~10 % of features
are confident. The deconvolution FDR / q-values are already computed by
FLASHDeconv and stored in the spectrum-level `*_ms1.tsv` files (per-observation
`Qvalue`, with a `FeatureIndex` linking each observation to a feature). The
charge-state count does **not** control FDR (decoys carry similar charge
envelopes), so confidence must come from the q-value.

`filter_features.py` maps each feature's best (minimum) q-value up from its
`*_ms1.tsv` observations and writes confident copies of the feature tables:

* keep `IsDecoy == 0` and `MSLevel == 1`
* keep feature q-value ≤ `FDR_THRESH` (default 0.05)
* keep `ChargeCount ≥ MIN_CHARGE_STATES` (default 3)

Output: `data/flashdeconv/filtered/<stem>_conf.tsv` (same columns/format). Re-run
MSTopDiff (GUI; no CLI) on these and place the resulting
`<stem>_conf_mstopdiff.csv` into `data/flashdeconv/mstopdiff/`.

### Dataset registry (`mstopdiff_config.py`)

Auto-discovers every `data/flashdeconv/mstopdiff/*_conf_mstopdiff.csv`, then labels,
colours and orders each dataset from `experimental_design.csv` (MCW = green,
Cytiva = blues, MagReSyn = reds, shaded by resuspension), flagging the three
publication datasets (MCW, Cytiva 0.5 % TFA, MagReSyn 0.5 % TFA). Both analysis
scripts import it, so adding a CSV needs no code change. **Note:** the MSTopDiff
CSV's date prefix must match the raw/design date or the dataset is skipped.

### Figures (`mstopdiff_compare.py`)

Reads the MSTopDiff CSV columns (`bin`, `count`, `intensity_lower_mass`,
`intensity_higher_mass`, and the `… x count` variants), rebins to 0.1 Da and
plots single-sided intensity×count histograms (0–90 Da, each normalised to its
own base peak). Outputs:

| File | Description |
| ---- | ----------- |
| `fig_mstopdiff_ixc.{png,pdf}` | Publication figure — stacked histograms for the three datasets |
| `fig_mstopdiff_ixc_all.{png,pdf}` | Stacked histograms for **all** discovered datasets |
| `fig_mstopdiff_mod_heatmap.{png,pdf}` | Modification × dataset heatmap of relative abundance |
| `fig_mstopdiff_mod_enrichment.{png,pdf}` | Relative-abundance bars (publication subset) |
| `mstopdiff_mod_table.csv` | Relative abundance (% of base peak) per modification per dataset |

### Unannotated-peak analysis (`mstopdiff_unannotated.py`)

Detects peaks ≥ `REL_THRESH` (default 10 %) of the base peak and classifies each
by its Δmass as **known** (single modification or pairwise combination),
**off-by-one satellite** of a known mass (±1.0023 Da deconvolution error),
**integer/isotope comb**, or **unannotated**, then reports the intensity-weighted
unexplained fraction (with a threshold sensitivity sweep). Outputs
`fig_mstopdiff_unannotated.{png,pdf}` and `mstopdiff_unannotated_*.csv`.

### Running

```bash
# confident features come from filter_features.py (deconvolution pipeline);
# run MSTopDiff (GUI) on data/flashdeconv/filtered/*_conf.tsv -> data/flashdeconv/mstopdiff/
conda run -n tdms python src/modifications/mstopdiff_compare.py
conda run -n tdms python src/modifications/mstopdiff_unannotated.py
```

---

## License

This project is released under an [MIT license](LICENSE) as-is, without warranty
of any kind.
