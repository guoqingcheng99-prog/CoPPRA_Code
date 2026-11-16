# CoPPRA analysis code

R, R Markdown, and Python code for CoPPRA proteomics analysis from DIA-NN reports. The workflow includes input preparation, quality control, peptide quantification, PCA, heatmaps, limma differential analysis, PTM enrichment, GO/KEGG enrichment, joint sample-label permutations, and numerical validation.

## Getting started

1. Follow the [DIA-NN setup guide](DIA-NN/README.md) to supply the required reports and sample metadata and install the R and Python dependencies.
2. Run the complete workflow from the repository root:

```sh
PYTHON="$PWD/DIA-NN/code/.venv-pca/bin/python" sh DIA-NN/code/run_all.sh
```

To list individual stages:

```sh
Rscript DIA-NN/code/run_analysis.R --help
```

The workflow processes one dataset using configurable input, results, and plot directories. Its sample design, search settings, analysis rules, and stage prerequisites are described in the [workflow documentation](DIA-NN/README.md).

## Repository contents

| Location | Contents |
| --- | --- |
| [`DIA-NN/code/`](DIA-NN/code/) | Analysis sources, plotting helpers, configuration, dependency checks, and runners. |
| [`DIA-NN/README.md`](DIA-NN/README.md) | Input requirements, environment setup, commands, and analysis definitions. |
| [`DIA-NN/reference_package_versions.csv`](DIA-NN/reference_package_versions.csv) | R package versions used for the validated workflow. |
| [`DIA-NN/limma_source.json`](DIA-NN/limma_source.json) | Source URL and checksum for the required limma version. |
| [`DIA-NN/VALIDATION.md`](DIA-NN/VALIDATION.md) | Execution checks, comparison with archived results, and live-annotation differences. |

Experimental input data, fitted results, and generated figures are not included. By default, inputs are read from `DIA-NN/Raw_Data/`, results are written to `DIA-NN/results/`, and figures are written to `DIA-NN/plot/`.

The complete workflow was validated through all 14 stages using 100 test permutations. The published default is 10,000 permutations; see the [validation record](DIA-NN/VALIDATION.md) for scope and results.
