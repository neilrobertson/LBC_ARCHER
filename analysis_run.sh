#!/bin/bash

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate ARCH_env

cd Scripts

# Run NGF analysis
python full_NGF_run.py
# Run threshold analysis
python full_threshold_run.py
# Run survival analysis
python survival_analysis.py
# Run NGF analysis
python further_analysis.py
