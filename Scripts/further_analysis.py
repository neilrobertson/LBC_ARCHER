# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================
"""Create and export figures for article.
"""
# =============================================================================
# Imports
# =============================================================================
import dill
import os
import json

import pandas as pd
import numpy as np
from math import isclose
import plotly.graph_objects as go
import plotly.express as px

colors = px.colors.qualitative.D3
# =============================================================================
# Local imports
# =============================================================================

# Append root directory to system's path
import sys
sys.path.append('../ARCH_package')

import plot
import modelling
import basic

# =============================================================================
# Set export path
# =============================================================================

# Load path to Other plots
path = '../Results/Other/'
if not os.path.exists(path):
    os.makedirs(path)

# =============================================================================
# Load datasets
# =============================================================================

# Import non-synonymoous trajectories as exported with basic.load module
with open('../Exports/LBC_non-synonymous.dill', 'rb') as infile:
    lbc = dill.load(infile)

# Import non-synonymoous trajectories as exported with basic.load module
with open('../Exports/LBC_synonymous.dill', 'rb') as infile:
    syn = dill.load(infile)

# load NGF fitted trajectories
with open('../Exports/neutral_filtered_trajectories.dill', 'rb') as infile:
    model = dill.load(infile)

# load Threshold fitted trajectories
with open('../Exports/threshold_filtered_trajectories.dill', 'rb') as infile:
    model_threshold = dill.load(infile)

# load NGF cohort
with open('../Exports/cohort_neutral.dill', 'rb') as infile:
    cohort = dill.load(infile)

# load NGF fit information
with open('../Exports/neutral_fit.dill', 'rb') as infile:
    neutral_fit = dill.load(infile)


# =============================================================================
# 2% gene trajectories and variant summary
# =============================================================================
df = pd.read_csv(r'../Datasets/LBC_ARCHER.2PCT_VAF.Mar21.non-synonymous.tsv',
                 sep='\t')

for gene in ['TET2', 'DNMT3A', 'JAK2']:
    fig = plot.gene_plot(df, gene)
    fig.write_image(path + f'2% {gene}.svg')

# Load variant color dictionary from resources
with open('../Resources/var_dict.json') as json_file:
    var_dict = json.load(json_file)

gene_dict = {element: 0
             for element in set(df.PreferredSymbol.unique())}
var_count_dict = {element: 0
                  for element in set(df.Variant_Classification.unique())}
# create list of all keys

for part in df.participant_id.unique():
    data = df[df['participant_id'] == part]
    for key in data.key.unique():
        data_key = data[data['key'] == key]
        gene = data_key.PreferredSymbol.unique()[0]
        variant = data_key.Variant_Classification.unique()[0]
        # update dict
        gene_dict[gene] += 1
        var_count_dict[variant] += 1

# sort dictionary in descending order
gene_dict = dict(sorted(gene_dict.items(),
                        key=lambda item: item[1], reverse=True))
var_count_dict = dict(sorted(var_count_dict.items(),
                      key=lambda item: item[1]))

# Bar plot
fig = go.Figure()
for key, item in gene_dict.items():
    fig.add_trace(
        go.Bar(x=[key], y=[item],
               name=f'{key}',
               marker_color='Grey',
               showlegend=False))
fig.update_layout(
            template="simple_white",
            yaxis_title='Count',
            xaxis_tickangle=-45)
fig.write_image(path + 'gene_counts_2%.svg', width=800, height=400)

var_name_dict = dict()
for key in var_count_dict.keys():
    new_key = key.replace('_', ' ')
    new_key = new_key.replace('Ins', 'Insertion')

    # create dictionary of names
    var_name_dict[key] = new_key

# Bar plot
fig = go.Figure()
for key, item in var_count_dict.items():
    fig.add_trace(
        go.Bar(y=[var_name_dict[key]], x=[item],
               marker_color=var_dict[key],
               showlegend=False,
               orientation='h'))
fig.update_layout(
            template="simple_white",
            xaxis_title='Count',
            )

fig.write_image(path + 'variant_counts_2%.svg', width=600)


# =============================================================================
# Variants distributions.
# =============================================================================
class participant:
    def __init__(self, id=None, max_fitness=0, max_VAF=0, max_aux=0):
        self.id = id
        self.max_fitness = max_fitness
        self.max_VAF = max_VAF
        self.max_aux = max_aux


part_list = []
for traj in model:
    part_list.append(traj.id)
part_list = list(set(part_list))

part_class = []
for part in part_list:
    part_class.append(participant(id=part))

# Extract maximum VAF, fitness, auxiliary for each participant
for traj in model:
    part = [x for x in part_class if x.id == traj.id][0]
    # update participant values
    mean_vaf = np.mean(traj.y)
    part.max_VAF = max(mean_vaf, part.max_VAF)
    part.max_fitness = max(traj.fitness, part.max_fitness)
    part.max_aux = max(traj.fitness*mean_vaf, part.max_aux)

vaf = [part.max_VAF for part in part_class]
fig = go.Figure(
        go.Box(y=vaf, boxmean=True, boxpoints='all'))
fig.update_layout(template='simple_white',
                  title='Highest mean VAF')
fig.update_yaxes(title='VAF')
fig.write_image(path + 'highest_VAF.svg', width=500)

fitness = [part.max_fitness for part in part_class]
print(f'Median maximum fitness per participant: {round(np.median(fitness),3)}')

fig = go.Figure(
        go.Box(y=fitness, boxmean=True, boxpoints='all'))
fig.update_layout(template='simple_white',
                  title='Highest Fitness')
fig.update_yaxes(title='Fitness')
fig.write_image(path + 'highest_fitness.svg', width=500)


aux = [part.max_aux for part in part_class]
print(f'Median maximum fitness*VAF per participant: {round(np.median(aux),3)}')

fig = go.Figure(
        go.Box(y=aux, boxmean=True, boxpoints='all'))
fig.update_layout(template='simple_white',
                  title='Highest fitnees x VAF')
fig.update_yaxes(title='fitness x VAF')
fig.write_image(path + 'highest_fitness_VAF.svg', width=500)

vaf_list = []
age_list = []
for part in syn:
    max_vaf = 0
    age = []
    for traj in part.trajectories:
        if np.mean(traj.data.AF) > max_vaf:
            # Update max_vaf
            max_vaf = np.mean(traj.data.AF)
            # append age for each trajectory
            max_age = np.mean(traj.data['age'])
    if max_vaf != 0:
        # Update vaf and age lists
        vaf_list.append(max_vaf)
        age_list.append(max_age)

print(f'Mean max_vaf: {np.mean(vaf_list)}')
print(f'Mean age: {np.mean(age_list)}')
fig = go.Figure(
        go.Box(y=vaf_list, boxmean=True, boxpoints='all'))
fig.update_layout(template='simple_white')
fig.update_yaxes(title='VAF')
fig.write_image(path + 'synonymous_size.svg', width=400)

# Set max_trajectory attribute for each participant
for part in lbc:
    part.max_trajectory = 0
    for traj in part.trajectories:
        mean_vaf = np.mean(traj.data.AF)
        # Exclude cases where mean(VAF)>0.5 as these are due to LOH
        if mean_vaf < 0.5:
            part.max_trajectory = max(part.max_trajectory, mean_vaf)


cohort_max_vaf = [part.max_trajectory for part in lbc
                  if part.max_trajectory > 0]

fig = go.Figure(
        go.Box(y=cohort_max_vaf, boxmean=True, boxpoints='all'))
fig.update_layout(template='simple_white',
                  title='Maximum mean VAF full cohort')
fig.update_yaxes(title='VAF')
fig.write_image(path+'cohort_max_vaf.svg', width=500)
print(f'Median of maximum mean VAF per participant: {round(np.median(cohort_max_vaf), 3)}')


# =============================================================================
# Grouped pathway analysis
# =============================================================================

# Create gene_to_category dictionary

category_to_gene = dict()
category_to_gene['histone_regulation'] = ['EZH2', 'ASXL1', 'KMT2A', 'KDM6A']
category_to_gene['splicing'] = ['SF3B1', 'U2AF1', 'SRSF2', 'U2AF2', 'ZRSR2']
category_to_gene['dna_damage'] = ['TP53', 'CDKN2A']
category_to_gene['mitogenic'] = ['KRAS', 'NF1', 'JAK2', 'JAK3']
category_to_gene['cohesin'] = ['RAD21', 'STAG2']
category_to_gene['tf_development'] = ['GATA2', 'RUNX1']
category_to_gene['dna_methylation'] = ['TET2', 'DNMT3A']

gene_to_category = dict()
for k, v in category_to_gene.items():
    for gene in v:
        gene_to_category[gene] = k

# Create fitness distribution by category
category_fitness = {element: [] for element in category_to_gene.keys()}
unclassified_genes = ['NOTCH1', 'NPM1', 'LUC7L2', 'DDX41', 'BCORL1', 'CUX1',
                      'PPM1D']
for traj in model:
    if traj.gene not in unclassified_genes:
        category_fitness[gene_to_category[traj.gene]].append(traj.fitness)

# sort dictionary by mean order
category_fitness = dict(sorted(category_fitness.items(),
                        key=lambda item: np.mean(item[1]),
                        reverse=True))

# Plot distribution of fitness by category
fig = go.Figure()
for i, key in enumerate(category_fitness):
    fig.add_trace(
            go.Box(y=category_fitness[key],
                   name=key, boxpoints='all', showlegend=False))
fig.update_xaxes(linewidth=2, tickangle=-45)
fig.update_layout(template='simple_white')
fig.update_yaxes(linewidth=2,
                 type='log', tickvals=[0.05, 0.1, 0.2, 0.4])
fig.write_image(path + 'gene_category.svg')

# Compute Kruskal Wallis test
fig, dict_2 = plot.gene_statistic(category_fitness, statistic='kruskal-wallis')
fig.write_image(path + 'gene_category_kruskal.svg', width=200)

# =============================================================================
# Range of N_cells
# =============================================================================


def parameter_surface(part_id, rtol=0.002):

    # Extract fitting data from participant
    data = neutral_fit[part_id][0]

    # Select subset of data to plot using relative tolerance
    new_data = data[np.isclose(data['error'], np.min(data.error), rtol=rtol)]

    # extract data presenting the minimum fitting error
    optimal = new_data[new_data['error'] == min(new_data['error'])]

    # extract dictionary optimal fitness for each variant
    optimal_dict = dict()
    for i, row in optimal.iterrows():
        optimal_dict[row['Gene name']] = row['Fitness']

    opt_dist = []
    for i, row in new_data.iterrows():
        opt_dist.append(np.abs(row['Fitness']
                        - optimal_dict[row['Gene name']]))
    new_data['fitness_distance'] = opt_dist

    fig = px.scatter(new_data, x='Origin', y='Cells',
                     hover_data=['Fitness', 'error'],
                     symbol='Gene name',
                     color='fitness_distance',
                     color_continuous_scale='Cividis')
    fig.update_xaxes(title='Ages (years)',
                     linewidth=2)
    fig.update_yaxes(title='HSPCs (counts)',
                     linewidth=2)
    fig.update_layout(template='simple_white')

    return fig


id = 'LBC0001A'
fig = parameter_surface(id, rtol=0.001)
fig.write_image(path + f'hypersurface of paramter solutions {id}.svg')

id = 'LBC0242M'
fig = parameter_surface(id)
fig.write_image(path + f'hypersurface of paramter solutions {id}.svg')
