# CoPPRA DIA-NN analysis

Downstream CoPPRA analysis from DIA-NN reports: quality control, peptide quantification, PCA, heatmaps, differential analysis, enrichment, and permutation testing. This workflow processes one dataset and writes outputs under `results/` and `plot/`.

The experimental design is fixed to Ctrl, R (ruxolitinib), U (ulixertinib), and RU (combination), each with three biological replicates measured in three technical runs, plus three reference biological samples measured in triplicate: 45 runs in total. Reference samples retain the input label `Carrier`; they are normalized separately and excluded from treatment comparisons.

## Inputs

Place these files in `DIA-NN/Raw_Data/`, or set `COPPRA_INPUT_DIR` to their directory:

| File | Purpose |
| --- | --- |
| `report.parquet` | Original DIA-NN report; `Precursor.Quantity` supplies unnormalized intensities. Required columns are listed in `code/00_prepare_inputs.py`. |
| `clean.parquet` | Companion cleaned report, used to verify the archived cleanup row count. |
| `cleanup_audit.json` | Cleanup provenance containing `source.bytes`, `counts.input_rows`, and `counts.clean_rows`. Values must describe these exact reports. QC independently reapplies the cleanup filters and checks the count. |
| `samples.tsv` | Tab-separated mapping with `Run` and `OriginalExperiment` columns. Experiment labels are `Ctrl1`–`Ctrl3`, `R1`–`R3`, `U1`–`U3`, `RU1`–`RU3`, and `Carrier1`–`Carrier3`. Each run name must end in `_1`, `_2`, or `_3` for its technical replicate. |
| `contaminant_ids.txt` | One contaminant accession per line, matching the report's protein identifiers. |
| `clean.protein_matrix.tsv` | Companion protein matrix, recorded as input provenance; peptide-level inference uses the report. |

The cleaned report and audit must reflect removal of decoys and contaminant-containing groups and a threshold of 0.01 for each of `Q.Value`, `Global.Q.Value`, `Global.PG.Q.Value`, `PG.Q.Value`, and `Lib.Peptidoform.Q.Value`. Input preparation checks report sizes and row counts, then records SHA-256 hashes of the six input files. It does not modify inputs.

The source reports were produced with DIA-NN 2.7.0, Trypsin/P, a maximum of **5 missed cleavages**, minimum peptide length 7, maximum neutral peptide mass 4,600 Da, precursor charges 1–4, fixed carbamidomethylation of cysteine, and variable methionine oxidation and protein N-terminal acetylation (maximum 5 variable modifications). The recorded MS1/MS2 tolerances were 4/7 ppm, with match-between-runs enabled. These are upstream search settings; this repository does not rerun DIA-NN or create a spectral library.

## Environment

The reference analysis used R 4.4.2. Package versions are recorded in [`reference_package_versions.csv`](reference_package_versions.csv). Install the R packages listed in `code/check_environment.R`, including `jsonlite` and `R.utils` for input handling. The environment check requires **limma 3.68.4**. Its source URL and SHA-256 are recorded in [`limma_source.json`](limma_source.json); a source archive can be installed into a local library with:

```sh
mkdir -p DIA-NN/code/R_library
R CMD INSTALL --library=DIA-NN/code/R_library /path/to/limma_3.68.4.tar.gz
```

Use Python 3.11 or newer. For example, from the repository root:

```sh
python3 -m venv DIA-NN/code/.venv-pca
DIA-NN/code/.venv-pca/bin/python -m pip install \
  -r DIA-NN/code/requirements.txt -r DIA-NN/code/requirements-pca.txt
```

The PCA renderer uses Matplotlib without a display server. GO annotations use installed Bioconductor annotation packages; KEGG enrichment queries the live service. A later enrichment run can change with annotation versions or service updates. The workflow saves the mappings and gene sets it uses.

## Run

From the repository root, after installing dependencies and supplying the inputs:

```sh
PYTHON="$PWD/DIA-NN/code/.venv-pca/bin/python" sh DIA-NN/code/run_all.sh
```

The shell runner checks R dependencies, prepares inputs, and runs all analysis stages, including 10,000 joint sample-label permutations and numerical validation. It writes package versions and per-stage session information, warnings, and completion logs. Rerunning a stage replaces its generated outputs.

To prepare inputs or run one stage separately:

```sh
DIA-NN/code/.venv-pca/bin/python DIA-NN/code/00_prepare_inputs.py
Rscript DIA-NN/code/run_analysis.R --help
Rscript DIA-NN/code/run_analysis.R 1.1
Rscript DIA-NN/code/run_analysis.R 2.2
Rscript DIA-NN/code/run_analysis.R validate
```

Stages require their upstream outputs. Stage order is `1.1`, `1.2`, `qc-table`, `2.1`, `2.2`, `heatmaps`, `3.1`, `3.2`, `5.1`, `5.2`, `7.1`, `7.2`, `permutation`, `validate`. R Markdown files are the canonical analysis sources; the runner extracts executable R scripts into `results/logs/`.

| Environment variable | Default / use |
| --- | --- |
| `COPPRA_PACKAGE` | `DIA-NN/`, resolved automatically by the runners. For direct stage scripts, run from `DIA-NN/code/`. |
| `COPPRA_INPUT_DIR` | `<package>/Raw_Data` |
| `COPPRA_RESULTS_DIR` | `<package>/results` |
| `COPPRA_PLOT_DIR` | `<package>/plot` |
| `COPPRA_R_LIBRARY` | `<package>/code/R_library`; supplements installed R libraries. |
| `COPPRA_PCA_PYTHON` | `<package>/code/.venv-pca/bin/python`; set an absolute path for another PCA Python environment. |
| `PYTHON` | `python3`; executable used by `run_all.sh` for input preparation. |
| `COPPRA_PERMUTATIONS` | `10000`; reduce only for a smoke test, not the reference analysis. |

Use absolute paths for overrides so Python and R resolve the same locations. The shell runner requires a POSIX shell; Windows users can invoke the Python and R runners separately and set `COPPRA_PCA_PYTHON` to their Python executable.

## Analysis rules

- Apply identification-quality filters, remove decoys and contaminants, reject nonpositive intensities and ambiguous precursor/run group assignments, then require at least two distinct unmodified peptide sequences per protein **before** the 20% run-detection filter. The protein filter is applied once.
- Collapse repeated modified-sequence/charge measurements by median log2 intensity. Center technical runs using precursors present in all three technical runs of each biological sample. Retain peptides detected in at least 9 of 45 runs and aggregate precursor/technical measurements by medians.
- Require observed values in at least two of three biological replicates in **every** treatment/control group before imputing missing values with each peptide's median across the 12 biological samples.
- Fit limma models; operational response calls use nominal P < 0.05 and absolute fitted log2 fold change > 1. BH-adjusted P values, two direct combination contrasts, and a factorial interaction contrast are also exported.
- The RU-only protein candidate set contains proteins with at least one responsive RU peptide and none under either monotherapy at those thresholds. This set definition does not establish synergy or direct binding. PTM enrichment tests annotated protein categories and does not measure site-specific PTM changes.
- Joint permutations preserve sample-level correlation and three samples per group, refit the empirical-Bayes prior, and compare responsive-peptide counts with the null. Eligibility and imputation remain fixed before permutation.

Key outputs include `QC/tables/QC_stepwise_peptide_protein_counts.csv`, `Limma/tables/all_peptide_differential_results.csv`, `Limma/tables/conjunctional_targeting_proteins.csv`, `Permutation/tables/global_permutation_summary.csv`, and `validation_summary.csv` under the results directory. The `PSEA` and `PEIMAN` folders contain rank-based and over-representation PTM enrichment, respectively; `TopGO` contains GO and KEGG results.

PCA plots include individual biological samples and arithmetic group means. The translucent triangles in the replicate view connect three samples and are not confidence regions. Camera settings and original scores are recorded with the output. Plotting changes do not alter fitted PCA coordinates.

The DIA-NN search, empirical library, and 20,431-entry human FASTA differ from the original MaxQuant analysis; numerical outputs should be attributed to this workflow.

See [`VALIDATION.md`](VALIDATION.md) for the execution checks, numerical comparison with archived results, and the live-annotation difference observed during preparation of this public workflow.
