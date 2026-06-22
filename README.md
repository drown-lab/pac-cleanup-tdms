# PAC for TDP

## Overview

This repository contains the custom Python used to support a study from the
Drown lab on bead-assisted **p**rotein **a**ggregation **c**apture (PAC) as a
sample-cleanup approach for **t**op-**d**own **p**roteomics (TDP). The code
accompanies the associated publication and is provided as-is to make the
analysis reproducible.

There are two independent scripts:

| Script | Purpose |
| ------ | ------- |
| `Feltenstein_GRAVY_pI_mass_calc.py` | Calculate biochemical properties (GRAVY, pI, aromaticity, instability) for a list of protein sequences. |
| `shared_proteoforms.py` | Compare identification scores of proteoforms shared across sample-cleanup conditions from a ProSightPD `.tdReport`. |

---

## Requirements

* Python 3.12 or later

Install the dependencies with conda:

```bash
conda env create -f environment.yml
conda activate tdms
```

---

## Script 1 — `Feltenstein_GRAVY_pI_mass_calc.py`

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
python Feltenstein_GRAVY_pI_mass_calc.py
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
