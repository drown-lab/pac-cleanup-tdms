# data/

Input data for the analysis scripts. **Everything in this directory is
git-ignored** — the files are large and are distributed through data
repositories, not GitHub. Recreate this layout locally to run the scripts.

| Path | Contents | Source |
| ---- | -------- | ------ |
| `raw/` | Thermo `.raw` acquisitions (one per cleanup condition) | ProteomeXchange/PRIDE accession **<ADD ACCESSION>** |
| `mzML/` | `.raw` converted to `.mzML` (ThermoRawFileParser) | derive from `raw/` |
| `flashdeconv/` | FLASHDeconv feature/spectrum TSVs, plus `filtered/` and `mstopdiff/` subfolders | produced by `src/deconvolution/` + MSTopDiff |
| `20250723_ifeltens_SP3_TDP_annotated_subsequence_2.tdReport` | ProSightPD identification report (SQLite) | ProSightPD search; **<ADD ACCESSION/DOI>** |

Expected sub-structure under `flashdeconv/` (see `src/deconvolution/` and
`src/modifications/` for how each is generated/consumed):

```
flashdeconv/
├── <stem>.tsv                 # FLASHDeconv feature tables
├── <stem>_ms1.tsv             # spectrum-level tables (carry the q-values)
├── filtered/<stem>_conf.tsv   # confident features  (filter_features.py)
└── mstopdiff/<stem>_conf_mstopdiff.csv   # MSTopDiff exports of the above
```

The sample-to-condition mapping lives in `config/experimental_design.csv`.
