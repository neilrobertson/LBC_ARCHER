#!/bin/bash
# Create and activate conda environment from yml file
echo "Installing and activating conda environment."

eval "$(conda shell.bash hook)"
conda env create --file env/ARCH_env.yml --name ARCH_env

conda activate ARCH_env

# add environment to jupyter notebooks
conda install -c anaconda ipykernel
python -m ipykernel install --user --name=ARCH_env

cd Scripts

# Run NGF analysis
python full_NGF_run.py
# Run threshold analysis
python full_threshold_run.py
# Run survival analysis
python survival_analysis.py
# Run NGF analysis
python further_analysis.py
