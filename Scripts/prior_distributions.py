# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-11-29
# =============================================================================
"""Create prior distributions for fitness, N_w and binomial artifact
probabilities.
"""

# =============================================================================
# region Import packages
# =============================================================================

# =============================================================================
# Import local packages
# =============================================================================

# Append root directory to system's path
import sys
sys.path.append('../ARCH_package')

# import packages from ARCH module
import filtering

# =============================================================================
# Import python packages
# =============================================================================

from scipy.stats import expon
from scipy.stats import gaussian_kde

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

import os
import dill
from tqdm import tqdm

pio.templates.default = "simple_white"

# endregion
# =============================================================================
# region Priors for fitness and N_w
# =============================================================================

# Fix resolution of priors
fitness_resolution = 50
N_w_resolution = 20
prior_resolution = 50

def create_uniform_prior(min_range, max_range, resolution, dtype='float'):
    range = np.linspace(min_range, max_range,
                        resolution, dtype=dtype)
    prob = [1 / (max_range-min_range)]*len(range)

    return np.vstack((range, prob))


def create_exponential_fitness_prior(min_fitness=0.05,
                                     max_fitness=0.5,
                                     resolution=fitness_resolution):
    """Approximation of the distribution of fitness
    reported in Blundell by an exponential function."""

    # Retrieve lambda from the equation
    # Exponential CDF valued at 0.04 (fitness) = 0.9
    lamb = -np.log(0.1)/0.04

    # Create array of fitness values
    fitness_range = np.linspace(min_fitness, max_fitness, resolution)
    # Retrieve associate pdf values
    fitness_pdf = expon.pdf(fitness_range, scale=1/lamb)
    fitness_prob = fitness_pdf/np.trapz(x=fitness_range, y=fitness_pdf)

    return np.vstack((fitness_range, fitness_prob))


print('Creating N_w and fitness priors')

flat_fitness_prior = create_uniform_prior(0.05, 0.5, fitness_resolution)
flat_N_w_prior = create_uniform_prior(
    10_000, 200_000, N_w_resolution, dtype='int')
exponential_fitness_prior = create_exponential_fitness_prior(fitness_resolution)

# endregion
# =============================================================================
# region Binomial artifact prior probabilies using synonymous mutations
# =============================================================================

# =============================================================================
# Define class objects and functions
# =============================================================================


class mutation ():
    """class object to bulk together data from the same mutation in
    different participants."""

    def __init__(self, id, counts, data):
        self.id = id
        self.counts = counts
        self.data = data


def extract_kde(data, quantiles=[0.01, 0.99], resolution=prior_resolution,
                flat_left_tail=False):
    """Create kde profile from data.
    Parameters:
    -----------
    data: List. List of data for which we want to extract its kde profile.
    quantiles: list. List of low and top quantiles range for kdee evaluation.

    Returns:
    kde_profil: Array. 2-d array containing kde profile.
    -----------

    """
    kernel_dist = gaussian_kde(data)

    # Extract support of distribution
    kernel_support = np.quantile(data, q=quantiles)

    # Sample uniformly from the range of binomial probabilities
    kernel_sample = np.linspace(kernel_support[0],
                                kernel_support[1],
                                resolution)

    # Compute associated prior probabilities using kde
    kernel_values = kernel_dist.evaluate(kernel_sample)

    # Flatten left tail of distribution
    if flat_left_tail is True:
        max_index = np.argmax(kernel_values)
        kernel_values[:max_index] = np.max(kernel_values)

    # normalise kernel to unit integral
    kernel_values = kernel_values / np.trapz(x=kernel_sample, y=kernel_values)

    return np.array([kernel_sample, kernel_values])


# =============================================================================
# Import Data
# =============================================================================

# Import non-synonymous trajectories as exported with basic.load
with open('../Exports/LBC_synonymous.dill', 'rb') as infile:
    syn = dill.load(infile)

# Create list of all mutations by cohort with repeats
total_mutations = dict()
for part in syn:
    for traj in part.trajectories:
        if traj.mutation not in total_mutations.keys():
            total_mutations[traj.mutation] = []
            total_mutations[traj.mutation].append(traj.data)
        else:
            total_mutations[traj.mutation].append(traj.data)

# Create a list of model comparison objects
mutation_list = []
for key, item in total_mutations.items():
    mutation_list.append(
        filtering.model_comparison(id=key, data=item)
    )


# # =============================================================================
# # Extract beta-binomial and binomial parameter's distributions
# # =============================================================================

print('Extracting beta binomial priors')
mutation_list_1 = [
    mutation for mutation in mutation_list if len(mutation.data) == 1]
mutation_list_2 = [
    mutation for mutation in mutation_list if len(mutation.data) > 1]

betabinom_p_prior = []
betabinom_beta_prior = []
for mutation_list in [mutation_list_1, mutation_list_2]:
    betabinom_models = [filtering.betabinom_artifact_fit(mutation.data)
                        for mutation in tqdm(mutation_list)]

    betabinom_p_list = [model.x[0] for model in betabinom_models]
    betabinom_beta_list = [model.x[1] for model in betabinom_models]

    # Extract beta binomial p and beta priors using kde
    betabinom_p_prior.append(extract_kde(betabinom_p_list,
                                         flat_left_tail=False))
    betabinom_beta_prior.append(extract_kde(betabinom_beta_list))


print('Extracting binomial prior')
binom_models = [filtering.binom_artifact_fit(mutation.data)
                for mutation in tqdm(mutation_list_1)]

binom_p_list = [model.x[0] for model in binom_models]

# Extract binomial proportion prior using kde
binom_p_prior = extract_kde(binom_p_list,
                            flat_left_tail=False)

# =============================================================================
# Export prior distributions
# =============================================================================

with open('../Resources/flat_N_w_prior.npy', 'wb') as f:
    np.save(f, flat_N_w_prior)

with open('../Resources/flat_fitness_prior.npy', 'wb') as f:
    np.save(f, flat_fitness_prior)

with open('../Resources/exponential_fitness_prior.npy', 'wb') as f:
    np.save(f, exponential_fitness_prior)

with open('../Resources/betabinom_p_prior.npy', 'wb') as f:
    np.save(f, betabinom_p_prior)

with open('../Resources/betabinom_beta_prior.npy', 'wb') as f:
    np.save(f, betabinom_beta_prior)

with open('../Resources/binom_p_prior.npy', 'wb') as f:
    np.save(f, binom_p_prior)
# endregion

# =============================================================================
# Plot prior distributions
# =============================================================================
# Create path for exporting
path = f'../Results/Priors/'
if not os.path.exists(path):
    os.makedirs(path)

prior_titles = ['beta binomial p prior', 'beta binomial beta prior',
                'binomial p prior']
priors = [betabinom_p_prior[1], betabinom_beta_prior[1], binom_p_prior]
xaxis_labels=['p', 'beta', 'p']

for prior, prior_title, xaxis_label in zip(priors, prior_titles, xaxis_labels):
    fig = px.line(x=prior[0,:], y=prior[1,:])
    fig.update_layout(title=prior_title,
                      xaxis_title=xaxis_label,
                      yaxis_title='Probability density')
    fig.write_image(path + prior_title + ".png", scale=10)
    fig.write_image(path + prior_title + ".svg")

# Plotting single vs multiple mutation occurrence prior comparison
prior_titles = ['Beta binomial p prior', 'Beta binomial beta prior']
priors = [betabinom_p_prior, betabinom_beta_prior]
xaxis_labels = ['p', 'beta']

for prior, title, label in zip([betabinom_p_prior, betabinom_beta_prior], prior_titles, xaxis_labels):
        
    single_traj_prior = prior[0]
    fig =go.Figure()
    fig.add_trace(
        go.Scatter(x=single_traj_prior[0,:],
                y=single_traj_prior[1,:],
                name='Single')
    )

    mult_traj_prior = prior[1]

    fig.add_trace(
        go.Scatter(x=mult_traj_prior[0,:],
                   y=mult_traj_prior[1,:],
                   name='Mutliple'))
    fig.update_layout(title=title,
                    xaxis_title=label,
                    yaxis_title='Probability density',
                    legend=dict(title='Mutation ocurrence',
                        xanchor="right",
                        x=0.95,
                    ))
    fig.show()
    fig.write_image(path + title + "_occurence_comparison.png", scale=10)
    fig.write_image(path + title + "_occurence_comparison.svg")