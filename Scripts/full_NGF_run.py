# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================
"""Script to load, filter and fit all trajectories in a cohort using NGF.
"""
# =============================================================================
# Import local packages
# =============================================================================
# Append root directory to system's path
import sys
sys.path.append('../ARCH_package')

import basic
import plot
import filtering
import modelling

# =============================================================================
# Import general packages
# =============================================================================

import os
import dill
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

import multiprocessing as mp
from functools import partial
from tqdm import tqdm


print('NGF filtering and fit start.')

# =============================================================================
# Set exporting paths
# =============================================================================

# Create global path for exporting plots
global_path = '../Results/NGF/'
if not os.path.exists(global_path):
    os.makedirs(global_path)
# Add prefix to plots
global_path = global_path + 'NGF-'

# Create path for exporting fit plots
fit_path = '../Results/NGF/cohort trajectories/'
if not os.path.exists(fit_path):
    os.makedirs(fit_path)
# Add prefix to plots
fit_path = fit_path + 'NGF-'

# Create path for gene trajectories
gene_path = '../Results/NGF/gene trajectories/'
if not os.path.exists(gene_path):
    os.makedirs(gene_path)
gene_path = gene_path + 'NGF-'

# =============================================================================
# Load cohort trajectories
# =============================================================================

print('Loading variant trajectories.')


# Load non-synonymous dataset
df = pd.read_csv(r'../Datasets/LBC_ARCHER.1PCT_VAF.Mar21.non-synonymous.tsv',
                 sep='\t')
lbc = basic.load(df, export_name='LBC_non-synonymous', create_dict=True)

# Load synonymous dataset
df = pd.read_csv(r'../Datasets/LBC_ARCHER.1PCT_VAF.Mar21.synonymous.tsv',
                 sep='\t')
syn = basic.load(df, export_name='LBC_synonymous', create_dict=False)

# =============================================================================
# Use NGF to filter trajectories
# =============================================================================

print('Training NGF and filtering variant trajectories.')


# Fit NGF and filter variants according to this filtering method
cohort = filtering.neutral_filter(lbc, syn)

# =============================================================================
# Export plots associated with the filter.
# =============================================================================

# Plot participant NGF trajectories
plot.participant_filter(lbc[18]).write_image(
    global_path + 'participant_filter_sample.svg')

# NGF stack plot
cohort.neutral_dist.figures[3].write_image(global_path + 'stack.svg')

# Filter on synonymoous variants plot
cohort.neutral_dist.figures[0].write_image(global_path + 'filter.svg',
                                           width=1000)
# Mean linear regression
cohort.neutral_dist.figures[1].write_image(
    global_path + 'mean_regression.svg', width=1000)

# Variance linear regression
cohort.neutral_dist.figures[2].write_image(
    global_path + 'Variance linear regression.svg', width=1000)

# Gene counts bar plot
cohort.gene_bar.write_image(global_path + 'gene_bar.svg')

# Gradient summary plot
cohort.gradient_plot.write_image(global_path + 'gradient_inset.png',
                                 width=600, scale=10)

# =============================================================================
# Export cohort after filtering.
# =============================================================================

with open('../Exports/cohort_neutral.dill', 'wb') as outfile:
    dill.dump(cohort, outfile)


# =============================================================================
# Fit trajectories with model of exponential growth.
# Trajectories are fitted in bulk by participant.
# =============================================================================

print('Fitting model of exponential growth.')

# Set base number of iterations for fitting
fit_iterations = 50

# Select all trajectories with more than 2 datapoints and mean <0.5
model = [traj for traj in cohort.model if len(traj.x) > 2
         and np.mean(traj.y) < 0.5]

# Create a dictionaary of trajectories by participant.
part_ids = set([traj.id for traj in model])
part_dict = {k: [] for k in part_ids}
for traj in model:
    part_dict[traj.id].append(traj)

# Create list of list where each sublist has all trajectories in a participant
fit_list = [part_dict[key] for key in list(part_dict.keys())]

# Fix method for fitting.
fit = partial(modelling.init_fit, method='least_squares',
              n_iterations=fit_iterations)

# Fit trajectories for each participant with multiprocessing.
if __name__ == '__main__':
    with mp.Pool() as p:
        result = list(tqdm(p.imap(fit, fit_list), total=len(fit_list)))

# Post processing of init_fit.
model, fit_dict = modelling.fit_postprocessing(result)

# Plot fitting results
for part_id in list(fit_dict.keys()):
    # Plot distribution of inferred fitness
    fit_fig = modelling.fitting_error_plot(part_id, part_dict, fit_dict)
    fit_fig.write_image(fit_path + f'{part_id} fitness_error_10%.svg')

    # Plot origin vs n_cells
    origin_fig = modelling.cells_origin(fit_dict, part_id)
    origin_fig.write_image(fit_path + f'{part_id} origin vs cells.svg')

    # Plot participant trajectories
    traj_fig = plot.participant_model_plot(model, part_id)
    traj_fig.write_image(fit_path + f'{part_id} trajectories.svg')

# Export fitted trajectories
with open('../Exports/neutral_trajectories.dill', 'wb') as outfile:
    dill.dump(model, outfile)

# Export fit information
with open('../Exports/neutral_fit.dill', 'wb') as outfile:
    dill.dump(fit_dict, outfile)


# =============================================================================
# Fit quality control.
# Exclude trajectories with bad fit and clinical records.
# =============================================================================
# Update participant_list
part_list = list(set([traj.id for traj in model]))
part_dict = {i: [] for i in part_list}
for traj in model:
    part_dict[traj.id].append(traj)

# Compute pearsonr coefficient for every fitted trajectory
for traj in model:
    traj.pearson = pearsonr(list(traj.data_vaf.values()),
                            list(traj.model_vaf.values()))

# Exclude trajectories with negative pearson coefficient
exclude_pearson = []
pearson = []
for traj in model:
    pearson.append(traj.pearson[0])
    if traj.pearson[0] < 0:
        exclude_pearson.append(traj.id)

# Extract error produced during fitting
error = []
for id, part in part_dict.items():
    error.append(part[0].fit.chisqr)

# Exclude participants with large fitting error (>1 std deviation)
exclude = []
for id, part in part_dict.items():
    if part[0].fit.chisqr > np.std(error):
        exclude.append(id)

# Exclude participants with a clinical record
clinical_records = ['LBC360021', 'LBC360725', 'LBC360914']
exclude.append(clinical_records[0])

# Exclude participants
model_filtered = [traj for traj in model if traj.id not in exclude]
model_filtered = [traj for traj in model_filtered if traj.pearson[0] > 0]

# Export fitted trajectories
with open('../Exports/neutral_filtered_trajectories.dill', 'wb') as outfile:
    dill.dump(model_filtered, outfile)

# =============================================================================
# Model analysis.
# =============================================================================

# Plot cohort overview
# # box plot with fitness distribution by fitness
box, gene_dict = plot.gene_box(model_filtered)
box.write_image(global_path + 'gene_box.svg', width=1000)

# Compute heatmap for statistic: 'kruskal-wallis' or 'anova'.
heatmap, statistic_df = plot.gene_statistic(gene_dict)
heatmap.write_image(global_path + 'gene_specific_differences.svg', height=300)

# Export trajectories grouped by gene
for gene in gene_dict.keys():
    fig = modelling.gene_trajectories(model_filtered, gene)
    fig.write_image(gene_path + f'{gene}.svg', width=1000)

# Damage prediction vs fitness plot
fig, df = plot.damage_class(model_filtered)
fig.write_image(global_path + 'damage_box.svg')