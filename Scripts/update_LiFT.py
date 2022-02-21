# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2022-01
# =============================================================================
"""Fit all trajectories selected as fit using LiFT filter.
"""
# =============================================================================
#region Initialisation
# =============================================================================

# =============================================================================
# Import local packages
# =============================================================================
# Append root directory to system's path

from platform import libc_ver
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

resolution_fitness = 201

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
with open('../Exports/LBC_non-synonymous_LiFT_fitted.dill', 'rb') as infile:
    lbc = dill.load(infile)

flat_fitness_prior = clonal_model.create_uniform_prior(0, 1, resolution_fitness)
flat_N_w_prior = clonal_model.create_uniform_prior(10_000, 200_000, resolution_N_w)


for part in lbc:
    clonal_model.determine_optimal_model(part)

flat_fitness_prior = clonal_model.create_uniform_prior(0, 1, resolution_fitness)
flat_N_w_prior = clonal_model.create_uniform_prior(10_000, 200_000, resolution_N_w)


# Set priors for model comparison 
postprocess = partial(clonal_model.optimal_model_post_processing,
                           fitness_prior=flat_fitness_prior,
                           N_w_prior=flat_N_w_prior)

print("Model comparison")
# Fit trajectories for each participant
if __name__ == '__main__':
    with Pool(n_cores) as p:
        lbc = list(tqdm(p.imap(postprocess, lbc), total=len(lbc)))

# Export fit
with open('../Exports/LBC_non-synonymous_LiFT_fitted.dill', 'wb') as outfile:
    dill.dump(lbc, outfile)
# =============================================================================
#region Plots
# =============================================================================
print('Exporting posterior distribution plots')

# Create path for exporting
path = f'../Results/LiFT/fitness posterior/'
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
path = f'../Results/LiFT/Genes/'
if not os.path.exists(path):
    os.makedirs(path)

for gene in gene_set:
    fig = plot.gene_trajectories(cohort, gene)
    fig.write_image(path + f"{gene}_deterministic_fit.svg")
    fig.write_image(path + f"{gene}_deterministic_fit.png", scale=5)


# =============================================================================
#region Damage analysis
# =============================================================================
print('Damage analysis')

# Create path for exporting
path = '../Results/LiFT/'
if not os.path.exists(path):
    os.makedirs(path)


fig = damage.make_damage_plot(lbc)
fig.write_image(path + 'damage_analysis.png', scale=10)
fig.write_image(path + 'damage_analysis.svg')


fig = damage.make_damage_plot_dist(lbc)
fig.write_image(path + 'damage_analysis_dist.png', scale=10)
fig.write_image(path + 'damage_analysis_dist.svg')

# =============================================================================
#region Gene statistics
# =============================================================================
print('Gene fitness statistics')
fig = plot.gene_fitness_plot(lbc)
fig.write_image(path + 'gene_fitness_summary.png', width=1000, scale=10)
fig.write_image(path + 'gene_fitness_summary.svg', width=1000)

gene_statistic_fig, gene_statistic_df = plot.gene_statistic(lbc)
gene_statistic_fig.show()

gene_statistic_fig.write_image(path + 'gene_statistic.png', scale=10)
gene_statistic_fig.write_image(path + 'gene_statistic.svg')