# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2022-01
# =============================================================================
"""Fit all trajectories selected as fit using 2% threshold filter.
"""
# =============================================================================
#region Initialisation
# =============================================================================

# =============================================================================
# Import local packages
# =============================================================================
# Append root directory to system's path

import sys

sys.path.append("../ARCH_package")

import clonal_model
import plot
import damage

# =============================================================================
# Import global packages
# =============================================================================

import os
import dill
from tqdm import tqdm

from functools import partial
from multiprocessing import Pool

# =============================================================================
# Global parameters
# =============================================================================

n_cores = 60

resolution_fitness_clonal_structure = 100
resolution_fitness_plots = 201

resolution_N_w = 21

# =============================================================================
# Check environment being used
# =============================================================================
if os.environ['CONDA_DEFAULT_ENV'] == 'ARCH':
     print('Conda environment ARCH is active')

else:
     print('Conda environment ARCH is not active')

# =============================================================================
# Import data
# =============================================================================

# Import non-synonymous trajectories as exported from LiFT
print("Importing cohort")
with open('../Exports/LBC_non-synonymous_threshold.dill', 'rb') as infile:
    lbc = dill.load(infile)
    
# Deal with incomplete mutation
for part in lbc:
    if part.id == 'LBC360919':
        for i in range(len(part.trajectories)):
            traj = part.trajectories[i]
            if traj.mutation != 'CEBPA c.542A>C':
                traj.data = traj.data[1:]


print("Preprocessing cohort")
# Process cohort for clonal model fitting
lbc = clonal_model.cohort_create_clonal_models(lbc)


# Create priors for fitness and N_w with desired resolution
flat_fitness_prior_clonal_structure = clonal_model.create_uniform_prior(0, 1, resolution_fitness_clonal_structure)
flat_fitness_prior_plots = clonal_model.create_uniform_prior(0, 1, resolution_fitness_plots)
flat_N_w_prior = clonal_model.create_uniform_prior(10_000, 200_000, resolution_N_w)


# Set priors for model comparison 
comparison_model = partial(clonal_model.bayesian_clonal_model_comparison,
                           fitness_prior_clonal_structure=flat_fitness_prior_clonal_structure,
                           fitness_prior_plots = flat_fitness_prior_plots,
                           N_w_prior=flat_N_w_prior)


print("Model comparison")
# Fit trajectories for each participant
if __name__ == '__main__':
    with Pool(n_cores) as p:
        lbc = list(tqdm(p.imap(comparison_model, lbc), total=len(lbc)))

# Export fit
with open('../Exports/LBC_non-synonymous_threshold_fitted.dill', 'wb') as outfile:
    dill.dump(lbc, outfile)
# =============================================================================
#region Plots
# =============================================================================
print('Exporting posterior distribution plots')

# Create path for exporting
path = f'../Results/Threshold/fitness posterior/'
if not os.path.exists(path):
    os.makedirs(path)

print("Exporting posterior fitness")
for part in lbc:
     fig = part.optimal_model.posterior_fitness_plot
     fig.write_image(path + f"{part.id}_fitness_distribution.svg")
     fig.write_image(path + f"{part.id}_fitness_distribution.png", scale=5)

     fig = part.optimal_model.contour_plot
     fig.write_image(path + f"{part.id}_contour_plot.svg")
     fig.write_image(path + f"{part.id}_contour_plot.png", scale=5)

print("Exporting deterministic fit plots")
for part in lbc:
    fig = part.deterministic_plot
    fig.write_image(path + f"{part.id}_deterministic_fit.svg")
    fig.write_image(path + f"{part.id}_deterministic_fit.png", scale=5)

cohort = [traj for part in lbc for traj in part.trajectories if traj.fitness >0.02]
gene_set = set([traj.gene for traj in cohort])

# Create path for exporting
path = f'../Results/Threshold/Genes/'
if not os.path.exists(path):
    os.makedirs(path)

for gene in gene_set:
    fig = plot.gene_trajectories(cohort, gene)
    fig.write_image(path + f"{gene}_deterministic_fit.svg")
    fig.write_image(path + f"{gene}_deterministic_fit.png", scale=5)


# =============================================================================
#region Gene statistics
# =============================================================================
# Create path for exporting
path = f'../Results/Threshold/'
if not os.path.exists(path):
    os.makedirs(path)

print('Gene fitness statistics')
fig = plot.gene_fitness_plot(lbc)
fig.write_image(path + 'gene_fitness_summary.png', scale=10)
fig.write_image(path + 'gene_fitness_summary.svg')

gene_statistic_fig, gene_statistic_df = plot.gene_statistic(lbc)
gene_statistic_fig.write_image(path + 'gene_statistic.png', scale=10)
gene_statistic_fig.write_image(path + 'gene_statistic.svg')