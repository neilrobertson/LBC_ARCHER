#!/bin/bash
# Create and activate conda environment from yml file
echo "Installing and activating conda environment."

# Activate conda environment
# conda info | grep -i 'base environment'
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ARCH

cd Scripts

# Exclude participants from study
python excluded_participants.py

# Cohort preprocessing
python cohort_preprocessing.py

# Creating distribution priors
echo "Creating prior_distributions"
python prior_distributions.py

# LiFT filter
echo "LiFT filter"
python LiFT_filter.py

echo "Clonal fit"
python LiFT_clonal_fit.py

echo "survival analysis"
python survival_analysis.py

# Threshold filter
echo "Threshold filter"
python threshold_filter.py

echo "Clonal fit"
python threshold_clonal_fit.py