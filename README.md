# Protein Property Calculator

## Overview

This Python script reads a CSV file containing protein sequences, calculates several biochemical properties for each sequence using Biopython, and saves the results to a new CSV file.

The calculated properties include:

* **GRAVY (Grand Average of Hydropathicity):** Measures the overall hydrophobicity or hydrophilicity of a protein.
* **Isoelectric Point (pI):** The pH at which the protein carries no net electrical charge.
* **Aromaticity:** The relative frequency of aromatic amino acids in the protein sequence.
* **Instability Index:** Predicts whether a protein is stable under in vitro conditions.

---

## Requirements

### Python Version

* Python 3.8 or later

### Required Libraries

Install the required packages using pip:

```bash
pip install pandas biopython
```

---

## Input File

The script expects a CSV file containing a column named:

| Column Name | Description                                       |
| ----------- | ------------------------------------------------- |
| Sequence    | Protein amino acid sequence in single-letter code |

Example:

```csv
Protein_ID,Sequence
P001,MKWVTFISLLFLFSSAYSRGVFRR
P002,MALWMRLLPLLALLALWGPGPG
```

---

## How It Works

1. Reads the input CSV file.
2. Extracts each protein sequence from the `Sequence` column.
3. Calculates:

   * GRAVY
   * Isoelectric Point
   * Aromaticity
   * Instability Index
4. Adds the calculated values as new columns.
5. Writes the updated data to a new CSV file.

---

## Output

The generated CSV file contains the original data plus four additional columns:

| New Column        | Description                       |
| ----------------- | --------------------------------- |
| GRAVY             | Hydrophobicity score              |
| Isoelectric Point | Predicted pI value                |
| Aromaticity       | Fraction of aromatic amino acids  |
| Instability       | Predicted protein stability index |

Example output:

```csv
Protein_ID,Sequence,GRAVY,Isoelectric Point,Aromaticity,Instability
P001,MKWVTFISLLFLFSSAYSRGVFRR,-0.256,8.62,0.120,33.45
```

---

## File Paths

The script currently uses the following hard-coded file paths:

```python
file_path = r'C:\Pycharm_Projects\20250722_HitReport_NoHistones.csv'
output_file_path = r'C:\Pycharm_Projects\20250722_HitReport_NoHistones_finished.csv'
```

Modify these paths as needed for your environment.

---

## Running the Script

Execute the script from the command line:

```bash
python protein_property_calculator.py
```

Upon completion, the script will display:

```text
Calculations completed and saved to the new Excel file.
```

**Note:** The output file is a CSV file, not an Excel workbook.

---

## Example Workflow

```text
Input CSV
      ↓
Read sequences
      ↓
Calculate protein properties
      ↓
Append new columns
      ↓
Save updated CSV
```

---

## Dependencies

* pandas
* Biopython (`Bio.SeqUtils.ProtParam.ProteinAnalysis`)

---

## License

This project is provided as-is for research and educational purposes.
