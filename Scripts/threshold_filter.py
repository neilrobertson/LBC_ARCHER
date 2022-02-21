# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2022-01
# =============================================================================
"""Use 2% VAF threshold to select all fit trajectories.
Export model comparison plots
"""
# =============================================================================
#region Initialisation
# =============================================================================

# =============================================================================
# Import local packages
# =============================================================================
# Append root directory to system's path
import sys
sys.path.append('../ARCH_package')

import filtering
import plot

# =============================================================================
# Import general packages
# =============================================================================
import dill
import os


# =============================================================================
# Check environment being used
# =============================================================================

if os.environ['CONDA_DEFAULT_ENV'] == 'ARCH':
    print('Conda environment ARCH is active')

else:
    print('Conda environment ARCH is not active')
    sys

# =============================================================================
# Global parameters
# =============================================================================

# set bayes factor threshold
VAF_threshold = 0.02

#endregion
# =============================================================================
#region Import Data
# =============================================================================

# Import non-synonymous trajectories as exported with basic.load
with open('../Exports/LBC_non-synonymous.dill', 'rb') as infile:
    lbc = dill.load(infile)

# Deal with error in LBC360919 CEBPA c.542A>C
# drop first time point in participant


for part in lbc:
    if part.id == 'LBC360919':
        for i in range(len(part.trajectories)):
            traj = part.trajectories[i]
            if traj.mutation == 'CEBPA c.542a>C':
                traj.data =  traj.data[traj.data.age > 70]

        break

for part in lbc:
    for traj in part.trajectories:
        traj.filter = False
        traj.optimal_model = 'below_2%'
        if max(traj.data.AF) > VAF_threshold:
            traj.filter = True
            traj.optimal_model = 'fit mutation'

fit_mutation_ids = [traj.mutation for part in lbc for traj in part.trajectories if traj.filter is True]
lbc_cohort = [traj for part in lbc for traj in part.trajectories]

with open('../Exports/LBC_non-synonymous_threshold.dill', 'wb') as outfile:
    dill.dump(lbc, outfile)

#endregion
# =============================================================================
#region Plots
# =============================================================================

# Create path for exporting
path = f'../Results/Threshold/'
if not os.path.exists(path):
    os.makedirs(path)

path = path + "2%_"

# =============================================================================
# Cohort
# =============================================================================

# Plot non-synonymous cohort
fig = plot.plot_cohort_LiFT(lbc_cohort, filter_type='2%')
fig.update_layout(title='Non-synonymous cohort LiFT')
fig.write_image(path + "non-synonymous.png", scale=10)
fig.write_image(path + "non-synonymous.svg")

# Cohort gradients
fig = plot.gradient_threshold_plot(lbc)
fig.write_image(path + 'gradient_vaf.png', scale=10)
fig.write_image(path + 'gradient_vaf.svg')

# # =============================================================================
# # Error measures concordance
# # =============================================================================

# Outlier p_value concordance
fig, outlier_pvalue_mutation_list = (
    plot.outlier_pvalue_comparison(lbc, fit_mutation_ids))
fig.write_image(path + 'vs_outlier_pvalue.png', scale=10)
fig.write_image(path + 'vs_outlier_pvalue.svg')
fig.write_html(path + 'vs_outlier_pvalue.html')

# MDAF concordance
fig, mdaf_mutation_list = plot.mdaf_comparison(lbc, fit_mutation_ids)
fig.write_image(path + 'vs_mdaf.png', scale=10)
fig.write_image(path + 'vs_mdaf.svg')
fig.write_html(path + 'vs_mdaf.html')



# # =============================================================================
# # Gene plot bar
# # =============================================================================

fit_cohort = filtering.model_cohort(lbc)

gene_bar = plot.gene_bar(fit_cohort)[0]
gene_bar.write_image(path + "gene_counts.png", scale=10)
gene_bar.write_image(path + "gene_counts.svg")