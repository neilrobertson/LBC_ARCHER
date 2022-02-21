# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2022-01
# =============================================================================
"""Use LiFT model comparison to select all fit trajectories.
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

import plotly.express as px

from multiprocessing import Pool
from tqdm import tqdm


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

# set number of cores used during mulitprocessing
n_cores = 60

# set bayes factor threshold
bayes_factor_threshold = 4

#endregion
# =============================================================================
#region Import Data
# =============================================================================

# Import non-synonymous trajectories as exported with basic.load
with open('../Exports/LBC_non-synonymous.dill', 'rb') as infile:
    lbc = dill.load(infile)

# Import non-synonymous trajectories as exported with basic.load
with open('../Exports/LBC_synonymous.dill', 'rb') as infile:
    syn = dill.load(infile)

#endregion
# =============================================================================
#region Model comparion
# =============================================================================


# =============================================================================
# Load mutation class objects
# =============================================================================

# Create list of all mutations by cohort with repeats
lbc_mutations = dict()
for part in lbc:
    for traj in part.trajectories:
        if traj.mutation not in  lbc_mutations.keys():
            lbc_mutations[traj.mutation] = []
            lbc_mutations[traj.mutation].append(traj.data)
        else:
            lbc_mutations[traj.mutation].append(traj.data)

lbc_list = []
for key, item in lbc_mutations.items():
    lbc_list.append(
        filtering.model_comparison(id=key, data=item)
    )

# Create list of all mutations by cohort with repeats
syn_mutations = dict()
for part in syn:
    for traj in part.trajectories:
        if traj.mutation not in  syn_mutations.keys():
            syn_mutations[traj.mutation] = []
            syn_mutations[traj.mutation].append(traj.data)
        else:
            syn_mutations[traj.mutation].append(traj.data)

syn_list = []
for key, item in syn_mutations.items():
    syn_list.append(
        filtering.model_comparison(id=key, data=item)
    )

# =============================================================================
# Bayesian model comparison with multiprocessing
# =============================================================================

# Bayesian comarison for non-synonymous mutations
with Pool(n_cores) as p:
    lbc_list = list(tqdm(
            p.imap(filtering.bayesian_model_comparison,
                   lbc_list),
                   total=len(lbc_list)))

# Bayesian comarison for synonymous mutations

with Pool(n_cores) as p:
    syn_list = list(tqdm(
            p.imap(filtering.bayesian_model_comparison,
                   syn_list),
                   total=len(syn_list)))

# =============================================================================
# Create cohort for further fitting
# =============================================================================

# Import non-synonymous trajectories as exported with basic.load
with open('../Exports/LBC_non-synonymous.dill', 'rb') as infile:
    lbc = dill.load(infile)


# Update trajectory type in lbc cohort
fit_mutation_ids = [mutation.id for mutation in lbc_list
                    if mutation.optimal_model =='fit mutation']

# Set filter attribute for trajectories in lbc cohort
for part in lbc:
    for traj in part.trajectories:
        if traj.mutation in fit_mutation_ids:
            traj.filter = True
            traj.optimal_model = 'fit mutation'
        else:
            traj.filter = False
            traj.optimal_model = 'artifact'

# =============================================================================
# Export model comparison objects
# =============================================================================

with open('../Exports/lbc_model_comparison.dill', 'wb') as outfile:
    dill.dump(lbc_list, outfile)

with open('../Exports/syn_model_comparison.dill', 'wb') as outfile:
    dill.dump(syn_list, outfile)

with open('../Exports/LBC_non-synonymous_LiFTed.dill', 'wb') as outfile:
    dill.dump(lbc, outfile)

#endregion
# =============================================================================
#region Plots
# =============================================================================

# Create path for exporting
path = f'../Results/LiFT/'
if not os.path.exists(path):
    os.makedirs(path)

# =============================================================================
# Bayes factor cut-off
# =============================================================================
# Plot effect of bayes factor cut-off on variant selection
fig = plot.bayes_factor_threshold_effect(syn_list)
fig.write_image(path + "bayes_factor_synonymous.png", scale=10)
fig.write_image(path + "bayes_factor_synonymous.svg")

fig = plot.bayes_factor_threshold_effect(lbc_list)
fig.write_image(path + "bayes_factor_non-synonymous.png", scale=10)
fig.write_image(path + "bayes_factor_non-synonymous.svg")

# =============================================================================
# Cohort
# =============================================================================

# Update bayes factor threshold in model comparison
for mutation_list in [lbc_list, syn_list]:
    for mutation in mutation_list:
        mutation.compute_optimal_model(binomial=True,
                                bayes_factor_threshold=bayes_factor_threshold)

# Plot synonymous cohort
fig = plot.plot_cohort_LiFT(syn_list)
fig.update_layout(title='Synonymous cohort LiFT')
fig.write_image(path + "LiFT_synonymous.png", scale=10)
fig.write_image(path + "LiFT_synonymous.svg")

# Plot non-synonymous cohort
fig = plot.plot_cohort_LiFT(lbc_list)
fig.update_layout(title='Non-ynonymous cohort LiFT')
fig.write_image(path + "LiFT_non-synonymous.png", scale=10)
fig.write_image(path + "LiFT_non-synonymous.svg")

# Cohort gradients
fig = plot.gradient_LiFT_plot(lbc, syn)
fig.write_image(path + 'LiFT_gradient_vaf.png', scale=10)
fig.write_image(path + 'LiFT_gradient_vaf.svg')

# =============================================================================
# Error measures concordance
# =============================================================================

# Outlier p_value concordance
fig, outlier_pvalue_mutation_list = (
    plot.outlier_pvalue_comparison(lbc, fit_mutation_ids))
fig.write_image(path + 'LiFT_vs_outlier_pvalue.png', scale=10)
fig.write_image(path + 'LiFT_vs_outlier_pvalue.svg')
fig.write_html(path + 'LiFT_vs_outlier_pvalue.html')

# MDAF concordance
fig, mdaf_mutation_list = plot.mdaf_comparison(lbc, fit_mutation_ids)
fig.write_image(path + 'LiFT_vs_mdaf.png', scale=10)
fig.write_image(path + 'LiFT_vs_mdaf.svg')
fig.write_html(path + 'LiFT_vs_mdaf.html')

# =============================================================================
# Gene plot bar
# =============================================================================

fit_cohort = filtering.model_cohort(lbc)

gene_bar = plot.gene_bar(fit_cohort)[0]
gene_bar.write_image(path + "LiFT_gene_counts.png", scale=10)
gene_bar.write_image(path + "LiFT_gene_counts.svg")