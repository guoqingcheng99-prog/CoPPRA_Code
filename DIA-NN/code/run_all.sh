#!/bin/sh
set -eu
TASK_CODE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$TASK_CODE"
COPPRA_PACKAGE=${COPPRA_PACKAGE:-$(dirname "$TASK_CODE")}
export COPPRA_PACKAGE
Rscript check_environment.R
# Set PYTHON to a Python >=3.11 executable with requirements.txt installed.
"${PYTHON:-python3}" 00_prepare_inputs.py
Rscript run_analysis.R all
