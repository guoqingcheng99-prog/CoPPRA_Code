# CoPPRA analysis code

R and R Markdown scripts for CoPPRA proteomics analyses.

## Contents

- [`DIA-NN/`](DIA-NN/README.md): current CoPPRA workflow for DIA-NN reports, including input preparation, QC, PCA, heatmaps, limma, PTM enrichment, GO/KEGG enrichment, joint sample-label permutations, and validation. Start here for the DIA-NN analysis.
- `Base/`: QC, PCA, limma, topGO, PEIMAN PSEA/SEA, and MGSA analyses and plots, plus shared configuration, input validation, and a runner.
- `SEA/`: SEA QC, PCA, limma, biological concept, CT comparison, search validation, and modification analyses and plots.

Experimental input data, fitted results, and generated figures are not included. The DIA-NN workflow uses a single input dataset with configurable paths; see its [setup and running instructions](DIA-NN/README.md).

## Original Base and SEA workflows

The original `Base/` and `SEA/` scripts are retained for provenance. Their input formats and layouts differ from the DIA-NN workflow.

### Inputs and configuration

Experimental input data and generated results are not included. The Base configuration refers to `evidence_new.txt` and `proteinGroups_base.txt`. Inspect the configuration and each script for required inputs, packages, sample labels, and paths before running.

`Base/CoPPRA_config.R` supports `COPPRA_PROJECT_DIR`, `COPPRA_BASE_OUTPUT`, `COPPRA_OUTPUT_ROOT`, and `COPPRA_BASE_PROTEIN_GROUPS`. Its default project directory is a local desktop path that must be adapted for another machine. `Base/run_Base.R` overrides `COPPRA_PROJECT_DIR` to the parent of its own directory and assumes the original `Base/Code/` layout; review this before using the runner with this repository layout.

### Running the analyses

Some scripts also depend on outputs from earlier stages and on the original project directory layout. Numbered filenames indicate the intended stage order; inspect each stage's prerequisites.

The original Base and SEA analyses have not been rerun as part of this update. Their supplied analysis scripts are unchanged.
