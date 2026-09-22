# CoPPRA analysis code

R and R Markdown scripts for CoPPRA proteomics analyses.

## Contents

- `Base/`: QC, PCA, limma, topGO, PEIMAN PSEA/SEA, and MGSA analyses and plots, plus shared configuration, input validation, and a runner.
- `SEA/`: SEA QC, PCA, limma, biological concept, CT comparison, search validation, and modification analyses and plots.

## Inputs and configuration

Experimental input data and generated results are not included. The Base configuration refers to `evidence_new.txt` and `proteinGroups_base.txt`. Inspect the configuration and each script for required inputs, packages, sample labels, and paths before running.

`Base/CoPPRA_config.R` supports `COPPRA_PROJECT_DIR`, `COPPRA_BASE_OUTPUT`, `COPPRA_OUTPUT_ROOT`, and `COPPRA_BASE_PROTEIN_GROUPS`. Its default project directory is a local desktop path that must be adapted for another machine. `Base/run_Base.R` overrides `COPPRA_PROJECT_DIR` to the parent of its own directory and assumes the original `Base/Code/` layout; review this before using the runner with this repository layout.

## Running the analyses

Some scripts also depend on outputs from earlier stages and on the original project directory layout. Numbered filenames indicate the intended stage order; inspect each stage's prerequisites.

The analyses have not been executed or validated as part of preparing this repository for sharing. The supplied analysis scripts are unchanged.
