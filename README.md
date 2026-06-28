# PAC for TDP

## Overview

This repository contains the custom Python used to support a study from the
Drown lab on bead-assisted **p**rotein **a**ggregation **c**apture (PAC) as a
sample-cleanup approach for **t**op-**d**own **p**roteomics (TDP). The code
accompanies the associated publication and is provided as-is to make the
analysis reproducible.

The repository contains several independent scripts:

| Script | Purpose |
| ------ | ------- |
| `proteoform_physiochemical_props.py` | Calculate biochemical properties (GRAVY, pI, aromaticity, instability) for a list of protein sequences. |
| `shared_proteoforms.py` | Compare identification scores of proteoforms shared across sample-cleanup conditions from a ProSightPD `.tdReport`. |
| `flashdeconv_summary.py` | Summarize FLASHDeconv deconvolution-FDR output (confident-mass yield, Qscore distributions) across spectrum-level `*_ms1.tsv` files. |
| `filter_features.py` | Filter FLASHDeconv feature TSVs to confident masses (drop decoys, q-value ≤ 5%, ≥3 charge states) for re-running MSTopDiff. |
| `mstopdiff_compare.py` | Build the MSTopDiff Δmass modification figures (per-condition intensity×count histograms, relative-abundance bars and heatmap). |
| `mstopdiff_unannotated.py` | Classify detected Δmass peaks as known / off-by-one satellite / isotope / unannotated and quantify the unexplained fraction. |
| `mstopdiff_config.py` | Shared dataset registry (auto-discovers MSTopDiff exports, labels/colours them from `experimental_design.csv`); imported by the two scripts above. |

---

## Requirements

* Python 3.12 or later

Install the dependencies with conda:

```bash
conda env create -f environment.yml
conda activate tdms
```

---

## Script 1 — `proteoform_physiochemical_props.py`

Reads a CSV file of protein sequences, calculates several biochemical properties
for each sequence using Biopython, and writes the results to a new CSV file.

### Calculated properties

* **GRAVY (Grand Average of Hydropathicity):** overall hydrophobicity/hydrophilicity of a protein.
* **Isoelectric Point (pI):** the pH at which the protein carries no net charge.
* **Aromaticity:** relative frequency of aromatic amino acids in the sequence.
* **Instability Index:** predicts whether a protein is stable in vitro.

### Input

A CSV file containing a column named `Sequence`:

| Column   | Description                                       |
| -------- | ------------------------------------------------- |
| Sequence | Protein amino acid sequence in single-letter code |

Example:

```csv
Protein_ID,Sequence
P001,MKWVTFISLLFLFSSAYSRGVFRR
P002,MALWMRLLPLLALLALWGPGPG
```

### Output

The original data plus four additional columns:

| New column        | Description                       |
| ----------------- | --------------------------------- |
| GRAVY             | Hydrophobicity score              |
| Isoelectric Point | Predicted pI value                |
| Aromaticity       | Fraction of aromatic amino acids  |
| Instability       | Predicted instability index       |

Example:

```csv
Protein_ID,Sequence,GRAVY,Isoelectric Point,Aromaticity,Instability
P001,MKWVTFISLLFLFSSAYSRGVFRR,-0.256,8.62,0.120,33.45
```

### Configuring file paths

The script uses hard-coded input/output paths near the top of the file. Edit
them for your environment before running:

```python
file_path = r'C:\Pycharm_Projects\name_of_your_document.csv'
output_file_path = r'C:\Pycharm_Projects\name_of_your_document_finished.csv'
```

### Running

```bash
conda run -n tdms python proteoform_physiochemical_props.py
```

On completion the script prints a confirmation message and writes the
`*_finished.csv` output file.

---

## Script 2 — `shared_proteoforms.py`

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

Both inputs must sit in the same directory as the script:

| File | Description |
| ---- | ----------- |
| `20250721_ifeltens_BEC_Consensus.tdReport` | ProSightPD consensus report (SQLite). **Not included in the repository** — supply your own. |
| `experimental_design.csv` | Maps each raw file to its `Cleanup` / `Resuspension` condition and an `Include` flag. |

`experimental_design.csv` has columns `Name, Cleanup, Resuspension, Include`. Rows
with `Include == FALSE` (e.g. files acquired with a different method) are dropped.

To compare different conditions, edit the `CONDITIONS_OF_INTEREST` list in the
script.

### Outputs

Written next to the script:

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
conda run -n tdms python shared_proteoforms.py
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

Output: `flashdeconv/filtered/<stem>_conf.tsv` (same columns/format). Re-run
MSTopDiff (GUI; no CLI) on these and place the resulting
`<stem>_conf_mstopdiff.csv` into `flashdeconv/mstopdiff/`.

### Dataset registry (`mstopdiff_config.py`)

Auto-discovers every `flashdeconv/mstopdiff/*_conf_mstopdiff.csv`, then labels,
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
conda run -n tdms python filter_features.py        # write confident _conf.tsv
# (re-run MSTopDiff on flashdeconv/filtered/*_conf.tsv -> flashdeconv/mstopdiff/)
conda run -n tdms python mstopdiff_compare.py
conda run -n tdms python mstopdiff_unannotated.py
```

---

## Dependencies

* pandas
* numpy
* scipy
* matplotlib
* seaborn
* matplotlib-venn
* Biopython (`Bio.SeqUtils.ProtParam.ProteinAnalysis`)

---

## License

This project is released under an [MIT license](LICENSE) as-is, without warranty
of any kind.
