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

import modelling
import filtering
import plot
import basic
# =============================================================================
# Import general packages
# =============================================================================
from tqdm import tqdm
from functools import partial
import multiprocessing as mp
from scipy.stats import pearsonr
import numpy as np
import pandas as pd
import argparse
import dill
import os
# =============================================================================
# Parse arguments
# =============================================================================

# Parse number of initialisations per trajectory during fitting
parser = argparse.ArgumentParser(
	description='Choose number of iterations for fitting process.')
parser.add_argument("--iterations", type=int, default=50,
                    help=('number of random initialisations per' +
                          'trajectory during the fitting process.'))

args = parser.parse_args()

# =============================================================================
# Starting message
# =============================================================================

print('NGF filtering and fit start.')
print(f'Number of fitting initialisations: {args.iterations}')

# =============================================================================
# Set exporting paths
# =============================================================================

# Create global path for exporting plots
global_path = '../Results/Threshold/'
if not os.path.exists(global_path):
    os.makedirs(global_path)
# Add prefix to plots
global_path = global_path + 'Threshold-'

# Create path for exporting fit plots
fit_path = '../Results/Threshold/cohort trajectories/'
if not os.path.exists(fit_path):
    os.makedirs(fit_path)
# Add prefix to plots
fit_path = fit_path + 'Threshold-'

# Create path for gene trajectories
gene_path = '../Results/Threshold/gene trajectories/'
if not os.path.exists(gene_path):
    os.makedirs(gene_path)
gene_path = gene_path + 'Threshold-'

# =============================================================================
# Load cohort trajectories
# =============================================================================
print('Loading variant trajectories.')

# Load non-synonymous dataset
df = pd.read_csv(r'../Datasets/LBC_ARCHER.1PCT_VAF.Mar21.non-synonymous.tsv',
                 sep='\t')
lbc = basic.load(df, export_name='LBC_non-synonymous', create_dict=True)

# =============================================================================
# Filter trajectories using a VAF threshold
# =============================================================================

print('Filtering variant trajectories.')

# Filter variants according to a threshold on VAF of 0.02
cohort = filtering.threshold_filter(lbc)

# =============================================================================
# Export plots associated with the filter.
# =============================================================================

# Plot participant NGF trajectories
plot.participant_filter(lbc[18]).write_image(
    global_path + 'participant_filtUse NGF to fer_sample.svg')

# Gene counts bar plot
cohort.gene_bar.write_image(global_path + 'gene_bar.svg')

# Gradient summary plot
cohort.gradient_plot.write_image(global_path + 'gradient_inset.svg',
                                 width=600)

# =============================================================================
# Export cohort after filtering.
# =============================================================================
with open('../Exports/cohort_threshold.dill', 'wb') as outfile:
    dill.dump(cohort, outfile)


# =============================================================================
# Fit trajectories with model of exponential growth.
# Trajectories are fitted in bulk by participant.
# =============================================================================

print('Fitting model of exponential growth.')

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
              n_iterations=args.iterations)

# Fit trajectories for each participant with multiprocessing.
if __name__ == '__main__':
    with mp.Pool() as p:
        result = list(tqdm(p.imap(fit, fit_list), total=len(fit_list)))

# Post processing of init_fit.
model, fit_dict = modelling.fit_postprocessing(result)

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
with open('../Exports/threshold_trajectories.dill', 'wb') as outfile:
    dill.dump(model, outfile)

# Export fit information
with open('../Exports/threshold_fit.dill', 'wb') as outfile:
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
with open('../Exports/threshold_filtered_trajectories.dill', 'wb') as outfile:
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
