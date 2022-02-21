# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================
"""Module containing plotting functions.
"""
# =============================================================================
# Import packages
# =============================================================================

import json
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations


import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

from plotly.validators.scatter.marker import SymbolValidator
import random
import json
# =============================================================================
# Global parameters
# =============================================================================
colors = px.colors.qualitative.D3
pio.templates.default = "simple_white"

# Load variant color dictionary
with open('../Resources/var_dict.json') as json_file:
    var_dict = json.load(json_file)

# Load gene color dictionary
with open('../Resources/gene_color_dict.json') as json_file:
    gene_dict = json.load(json_file)

# =============================================================================
# Plotting functions for participants
# =============================================================================


def profile(self, germline=True):
    """Plot all variant trajectories in a participant class object.
    Parameters:
    - germilne: Boolean. Include trajectories with
                         attribute germline = True.
    Return: plotly graphical object.
    """

    # Initialize figure
    fig = go.Figure()

    # If germline is True  plot all all trajectories
    if germline is True:
        for traj in self.trajectories:
            fig.add_trace(go.Scatter(x=traj.data.age, y=traj.data.AF,
                                     mode='lines+markers',
                                     name=traj.mutation))

    # If germline is False plot trajectories with germline attribute==False
    else:
        for traj in self.trajectories:
            if traj.germline is False:
                fig.add_trace(go.Scatter(x=traj.data.age, y=traj.data.AF,
                                         mode='lines+markers',
                                         name=traj.mutation))
    # Update figure layout
    fig.update_layout(title=f'Trajectories of participant {self.id}',
                      xaxis_title='Age (in years)',
                      yaxis_title='VAF')
    return fig


def plot_id(cohort, participant_id, germline=False):
    """Given a participant's id of a cohort, plot all its variant trajectories.
    Parameters:
    - cohort: list of participant class objects. Cohort where we search
                                                 for the participant.
    - participant_id: string. Participant's id.
    - germilne: Boolean. Include trajectories with
                         attribute germline = True.
    Return: plotly graphical object.
    """
    # Plot trajectories by participant_id
    for part in cohort:
        if part.id == participant_id:
            fig = part.profile(germline=germline)
    return fig


def synonymous_profile(part_lbc, syn):

    fig = go.Figure()
    # Plot non-synonymous trajectory for legend
    traj = part_lbc.trajectories[0]
    fig.add_trace(
            go.Scatter(x=traj.data.age,
                       y=traj.data.AF,
                       marker_color=colors[0],
                       opacity=0.5,
                       name='Non-Synonymous variant',
                       legendgroup='Non-synonymous variant'))
    # Plot all non-synonymous trajectories
    for traj in part_lbc.trajectories:
        fig.add_trace(
            go.Scatter(x=traj.data.age,
                       y=traj.data.AF,
                       marker_color=colors[0],
                       opacity=0.3,
                       showlegend=False,
                       legendgroup='Non-synonymous variant'))
    # Find synonymous trajectories of participant
    for part in syn:
        if part.id == part_lbc.id:
            # Plot synonymous trajectory for legend
            traj = part.trajectories[0]
            fig.add_trace(
                go.Scatter(x=traj.data.age,
                           y=traj.data.AF,
                           marker_color='Orange',
                           name='Synonymous variant',
                           legendgroup='Synonymous variant'))
            # Plot all synonymous trajectories
            for traj in part.trajectories:
                fig.add_trace(
                    go.Scatter(x=traj.data.age,
                               y=traj.data.AF,
                               marker_color='Orange',
                               showlegend=False,
                               legendgroup='Synonymous variant'))
    fig.update_layout(title='Synonymous mutations',
                      xaxis_title='Years',
                      yaxis_title='VAF',
                      template='plotly_white',
                      legend=dict(
                          y=0.95,
                          x=0.1,
                      ))

    return fig


# =============================================================================
# Plotting functions for longitudinal trajectories
# =============================================================================


def mutation(cohort, mutation):
    # plot all trajectories with a mutations
    fig = go.Figure()

    for part in cohort:
        for i, word in enumerate(part.cohort):
            if mutation in word.split():
                traj = part.trajectories[i]
                fig.add_trace(go.Scatter(x=traj.data.age,
                                         y=traj.data.AF,
                                         mode='lines+markers',
                                         name=traj.mutation,
                                         hovertemplate=f"{part.id}"
                                         ))
    # Edit the layout
    fig.update_layout(title=f'Trajectories containing mutation {mutation}',
                      xaxis_title='Time (years since first age)',
                      yaxis_title='VAF')

    return fig


def gene_plot(df, gene):
    # filter by gene
    data_gene = df[df['PreferredSymbol'] == gene]

    # Create figure
    fig = go.Figure()
    # create list of all keys
    participants = data_gene.participant_id.unique()
    for part in participants:
        data_part = data_gene[data_gene['participant_id'] == part]
        keys = data_part.key.unique()
        for key in keys:
            data = data_part[data_part['key'] == key]
            v_type = data['Variant_Classification'].unique()[0]
            if gene == 'TET2' and max(data.AF) > 0.4:
                continue
            fig.add_trace(
                go.Scatter(x=data['wave'], y=data['AF'],
                           marker_size=10,
                           marker_color=var_dict[v_type], showlegend=False))
    fig.update_layout(template='simple_white')
    fig.update_yaxes(title='VAF',
                     linewidth=2,
                     dtick=0.1)
    fig.update_layout(title=gene,
                      xaxis=dict(linewidth=2,
                                 tickmode='array',
                                 tickvals=[1, 2, 3, 4, 5],
                                 ticktext=['1 <br>~70 years<br>~79 years',
                                           '2 <br>~73 years<br>~82 years',
                                           '3 <br>~76 years<br>~85 years',
                                           '4 <br>~79 years<br>~88 years',
                                           '5 <br>~82 years<br>~91 years']))

    return fig


def participant_model_plot(model_list, id):
    """ Returns scatter plot of data and model predictions
    of fit trajectories in a participant"""

    part = []
    for traj in model_list:
        if traj.id == id:
            part.append(traj)
    if len(part) == 0:
        return go.Figure()

    # Extract min and max time
    min_time = []
    max_time = []
    for traj in part:
        min_time.append(min(list(traj.data_vaf.keys())))
        max_time.append(max(list(traj.data_vaf.keys())))
    min_time = min(min_time)
    max_time = min(max_time)

    fig = go.Figure()
    for i, traj in enumerate(part):
        x = list(traj.data_vaf.keys())
        y = list(traj.data_vaf.values())
        fig.add_trace(
            go.Scatter(x=x, y=y,
                       mode='markers',
                       marker_color=var_dict[traj.variant_class],
                       name=traj.mutation))
        # x_prediction = list(traj.vaf_plot.keys())
        # y_prediction = list(traj.vaf_plot.values())
        x_prediction = [time
                        for time in list(traj.vaf_plot.keys())
                        if min_time - 3 < time < max_time + 3]
        y_prediction = [traj.vaf_plot[time]
                        for time in list(traj.vaf_plot.keys())
                        if min_time - 3 < time < max_time + 3]
        fig.add_trace(
            go.Scatter(x=x_prediction,
                       y=y_prediction, mode='lines',
                       marker_color=var_dict[traj.variant_class],
                       text=(f'id: {traj.id}<br>'
                             f'fitness: {round(traj.fitness,3)}<br>'
                             f'origin: {round(traj.origin,3)}<br>'
                             f'r2: {round(traj.r2,3)}'),
                       name=traj.mutation))
    fig.update_layout(
        title=(f'Trajectory fit of participant {part[0].id} <br>'
               f'aic: {int(traj.fit.aic)}'),
        xaxis_title='Age (in years)',
        yaxis_title='VAF')

    # fig.update_xaxes(range=[min_time - 3, max_time + 3])

    return fig


# =============================================================================
# Plotting functions for cohort statistics
# =============================================================================


def top_bar(cohort, n_genes=10, all=False):
    gene_list = []
    for part in cohort:
        for traj in part.trajectories:
            gene_list.append(traj.mutation.split()[0])
    gene_dict = {element: 0 for element in set(gene_list)}
    for part in cohort:
        for traj in part.trajectories:
            gene_dict[traj.mutation.split()[0]] = gene_dict[traj.mutation.split()[
                                          0]] + 1
    gene_dict = dict(sorted(gene_dict.items(),
                            key=lambda item: item[1], reverse=True))

    if all is False:
        # Filter top mutations
        top_genes = list(gene_dict.keys())[0:n_genes]
        gene_dict = {gene: gene_dict[gene] for gene in top_genes}

    # Bar plot
    fig = go.Figure([go.Bar(x=list(gene_dict.keys()),
                            y=list(gene_dict.values()))])
    return fig


def gradients(cohort, mutations):
    # violin plots of gradients by mutation
    data = pd.DataFrame(columns=['gradient', 'participant', 'mutation'])

    for part in cohort:
        for traj in part.trajectories:
            if traj.mutation.split()[0] in mutations:
                data = data.append({'gradient': traj.gradient,
                                    'participant': part.id,
                                    'mutation': traj.mutation.split()[0]},
                                   ignore_index=True)
    # violin plot of data
    fig = px.box(data, y="gradient", x='mutation',
                 color='mutation', points='all',
                 hover_name="participant", hover_data=["mutation"])

    fig.update_layout(title='Trajectory gradients by mutation')
    return fig


def gene_box(cohort, order='median', scale='log'):
    """Box plot with counts of filtered mutations by gene.
       percentage computes fitness as the increase with respect to
       the self-renewing replication rate lambda=1.3.
       Color allows you to use a dictionary of colors by gene.
       Returns a figure."""

    # Load gene color dictionary
    with open('../Resources/gene_color_dict.json') as json_file:
        color_dict = json.load(json_file)
    # Create a dictionary with all filtered genes
    gene_list = []
    for traj in cohort:
        gene_list.append(traj.gene)
    gene_dict = {element: [] for element in set(gene_list)}

    # update the counts for each gene
    y_label = 'Fitness'
    for traj in cohort:
        fitness = traj.fitness
        gene_dict[traj.gene].append(fitness)

    # sort dictionary in descending order
    if order == 'mean':
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: np.mean(item[1]),
                                reverse=True))

    if order == 'median':
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: np.median(item[1]),
                                reverse=True))
    if order == 'max':
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: np.max(item[1]),
                                reverse=True))
    # Bar plot
    fig = go.Figure()
    # color_dict = dict()
    # if isinstance(color, dict):
    #     color_dict = color
    for i, key in enumerate(gene_dict):
        fig.add_trace(
            go.Box(y=gene_dict[key],
                   marker_color=color_dict[key],
                   name=key, boxpoints='all', showlegend=False))

    fig.update_layout(title='Gene distribution of filtered mutations',
                      yaxis_title=y_label,
                      template="simple_white")

    fig.update_xaxes(linewidth=2)
    fig.update_yaxes(linewidth=2)

    if scale is 'log':
        fig.update_yaxes(type=scale, tickvals=[0.05, 0.1, 0.2, 0.4])
    
    fig.update_layout(xaxis_tickangle=-45)

    return fig, gene_dict


def gene_statistic(lbc, statistic='kruskal-wallis', fitness_threshold = 0.02, filter=True):
    """ compute a statistical test to find significant differences in the
    distribution of fitness by gene.
    statistic parameter accepts: 'kruskal' or 'anova'.
    Returns:
    * heatmap with significant statistical differences.
    * dataframe."""

    cohort = [traj for part in lbc for traj in part.trajectories if traj.fitness>fitness_threshold]

    # Create gene_dict
    gene_fitness_list = [(traj.gene, traj.fitness) for traj in cohort]

    gene_fitness_dict =dict()
    for x, y in gene_fitness_list:
        gene_fitness_dict.setdefault(x, []).append(y)

    gene_fitness_dict = dict(sorted(gene_fitness_dict.items(), key=lambda item: np.median(item[1]), reverse=True))
    # Check if statistic is allowed
    if statistic not in ['kruskal-wallis', 'anova']:
        return 'Statistic not recognised.'

    # extract all possible gene combinations
    gene_list = []
    for gene in gene_fitness_dict.keys():
        if len(gene_fitness_dict[gene]) > 1:
            gene_list.append(gene)

    # Create dataframe to store statistics
    test_df = pd.DataFrame(index=gene_list, columns=gene_list)
    for gene1, gene2 in combinations(gene_list, 2):
        # compute statistic for each possible comination of genes
        if statistic == 'kruskal-wallis':
            stat, pvalue = stats.kruskal(gene_fitness_dict[gene1], gene_fitness_dict[gene2])
        if statistic == 'anova':
            stat, pvalue = stats.f_oneway(gene_fitness_dict[gene1], gene_fitness_dict[gene2])
        # if statistic is significant store value in dataframe
        if pvalue < 0.05:
            test_df.loc[gene1, gene2] = stat

    # Clean dataset from nan
    if filter is True:
        test_df = test_df.dropna(how='all', axis=1)
        test_df = test_df.dropna(how='all', axis=0)

    test_df = test_df.reindex(index=test_df.index[::-1])
    y = test_df.index
    x = test_df.columns
    fig = go.Figure(data=go.Heatmap(
                    z=np.array(test_df),
                    x=x,
                    y=y,
                    colorscale='Cividis',
                    colorbar=dict(title=f'{statistic} score')))
    fig.update_xaxes(side="top", mirror=True)
    fig.update_yaxes(side='top', mirror=True)
    fig.update_layout(template='simple_white')
    fig.update_layout(margin=dict(pad=4))

    return fig, test_df


# =============================================================================
# Plotting functions for protein damage prediction
# =============================================================================


def damage_class(model):
    """ Box plot of variant predicted protein damage class ~ Fitness.
    Parameters:
    model: List. List of all fitted trajectories
    Returns:
    * Box plot
    * Modified model with damage_class attribute."""

    # Set default color
    default_color = 'Grey'

    # Load dataset with prediction damage
    xls = pd.ExcelFile('../Datasets/variants_damage.xlsx')
    df = pd.DataFrame(columns=['p_key', 'damaging'])

    damaging_dict = {1: 'likely damaging',
                     0: 'possibly damaging',
                     -1: 'likely benign'}
    # Access each sheet of xls file
    for name in xls.sheet_names:
        df_temp = pd.read_excel(xls, name)
        # replace integer for damage_class
        for key, value in damaging_dict.items():
            df_temp['Likely damaging'] = (df_temp['Likely damaging'].
                                          replace(key, value))

        # Extract p_key and damage class
        p_keys = name + ' p.' + df_temp.iloc[:, 0]
        damaging = df_temp['Likely damaging']
        for i, j in zip(p_keys, damaging):
            df = df.append({'p_key': i, 'damaging': j}, ignore_index=True)

    # Exclude trajectories in model without a p_key
    damage_model = [traj for traj in model if isinstance(traj.p_key, str)]

    # Assign damage to trajectories
    for traj in damage_model:
        if traj.p_key in set(df['p_key']):
            traj.damage_class = df.loc[df.p_key == traj.p_key,
                                       'damaging'].values[0]
        elif 'fs' in traj.p_key or '*' in traj.p_key:
            traj.damage_class = 'fs - ter'
        else:
            traj.damage_class = None

    # filter trajectories in model without damage_class
    damage_model = [traj for traj in damage_model
                    if traj.damage_class is not None]

    # Create dataframe for plotting
    df = pd.DataFrame(columns=['damage_class', 'fitness'])
    for traj in damage_model:
        df = df.append({'damage_class': traj.damage_class,
                        'fitness': traj.fitness*1}, ignore_index=True)

    # order df by damage class in a specific order
    df.damage_class = pd.Categorical(df.damage_class,
                                     categories=['likely benign',
                                                 'possibly damaging',
                                                 'likely damaging',
                                                 'fs - ter'],
                                     ordered=True)

    df.sort_values('damage_class', inplace=True)

    # Start subplot figure with shared xaxis
    fig = make_subplots(rows=2, cols=1,
                        row_heights=[0.5, 0.5],
                        shared_xaxes=True,
                        vertical_spacing=0.05)

    fig.add_trace(
        go.Histogram(x=df['damage_class'],
                     marker_color=default_color,
                     showlegend=False),
        row=1, col=1)

    fig.add_trace(
        go.Box(x=df['damage_class'],
               y=df['fitness'],
               marker_color=default_color,
               boxpoints='all',
               showlegend=False),
        row=2, col=1)

    fig.update_layout(
        template='simple_white',
        boxmode='group'
    )
    fig.update_yaxes(title='fitness', type='log',
                     dtick=[0.05, 0.1, 0.2, 0.4], row=2, col=1)
    fig.update_yaxes(title='trajectory counts', row=1, col=1)
    fig.update_yaxes(title='fitness',
                     linewidth=2,
                     type='log',
                     tickvals=[0.05, 0.1, 0.2, 0.4],
                     row=2, col=1)
    fig.update_yaxes(linewidth=2, row=1, col=1)
    fig.update_xaxes(linewidth=2, row=2, col=1)
    return fig, df

#
# def damage(model, fs_score=-1, termination_score=-1,
#            unknown_score=-1):
#     """ Scatter plot of variant predicted protein damage ~ Fitness
#         Returns:
#         * Scatter plot
#         * Modified model with damage_class attribute."""
#
#     # Load xls file
#     xls = pd.ExcelFile('Datasets/new_kristina_variants.xlsx')
#     df = pd.DataFrame(columns=['p_key', 'damaging'])
#     # For each sheet, compute the mean of gnomad and clinvar
#     for name in xls.sheet_names:
#         df_temp = pd.read_excel(xls, name)
#         df_temp.columns = [col.lower() for col in df_temp.columns]
#
#         # Detect the presence of gnomad mean and clinvar mean columns
#         predictors = {'gnomad mean', 'clinvar mean'}
#         col = set(df_temp.columns)
#         inter = predictors.intersection(col)
#
#         # calculate the mean were possible
#         df_temp['score'] = df_temp[inter].mean(axis=1)
#
#         # extract p_keys and score for each variant
#         p_keys = name + ' p.' + df_temp.iloc[:, 0]
#         score = df_temp['score']
#         for i, j in zip(p_keys, score):
#             df = df.append({'p_key': i, 'score': j}, ignore_index=True)
#
#     # Exclude trajectories without a p_key
#     damage_model = [traj for traj in model if isinstance(traj.p_key, str)]
#
#     # Assign damage score to trajectories
#     for traj in damage_model:
#         if traj.p_key in set(df['p_key']):
#             traj.damage_score = df.loc[df.p_key == traj.p_key,
#                                        'score'].values[0]
#         elif 'fs' in traj.p_key:
#             traj.damage_score = fs_score
#         elif '*' in traj.p_key:
#             traj.damage_score = termination_score
#         else:
#             traj.damage_score = unknown_score
#
#     # Scatter plot of fitness vs score
#     fitness = []
#     damage = []
#     scatter_fig = go.Figure()
#     for traj in damage_model:
#         if 0 < traj.damage_score <= 1:
#             scatter_fig.add_trace(
#                 go.Scatter(x=[traj.fitness], y=[traj.damage_score],
#                            marker_color=colors[0],
#                            showlegend=False))
#             fitness.append(traj.fitness)
#             damage.append(traj.damage_score)
#
#     tau, p_value = stats.kendalltau(fitness, damage)
#
#     scatter_fig.update_layout(
#         xaxis_title='Fitness',
#         yaxis_title='damage prediction',
#         title=(f'Fitness ~ damage prediction: '
#                f'tau: {round(tau,2)} p-value: {round(p_value,3)}'))
#
#     # Split trajectories by damage calss 1, 0, -1
#     damaging_var = []
#     possibly_var = []
#     benign_var = []
#     for traj in damage_model:
#         if traj.damage_score == 1:
#             damaging_var.append(traj.nelder.params['fitness'])
#             traj.damage_class = 'Likely damaging'
#         elif traj.damage_score == 0:
#             possibly_var.append(traj.nelder.params['fitness'])
#             traj.damage_class = 'Possibly damaging'
#         elif traj.damage_score == -1:
#             benign_var.append(traj.nelder.params['fitness'])
#             traj.damage_class = 'Likely benign'
#         else:
#             traj.damage_class = 'Unknown'
#
#     # Box plot of fitness by damage class
#     box = go.Figure()
#     box.add_trace(
#         go.Box(y=damaging_var, name='Likely damaging',
#                boxpoints='all', marker_color='rgb(102,102,102)'))
#     box.add_trace(
#         go.Box(y=possibly_var, name='Possibly damaging',
#                boxpoints='all', marker_color='rgb(102,102,102)'))
#     box.add_trace(
#         go.Box(y=benign_var, name='Likely benign',
#                boxpoints='all', marker_color='rgb(102,102,102)'))
#
#     box.update_yaxes(title='Fitness')
#     box.update_layout(showlegend=False)
#
#     return scatter_fig, box, damage_model
#

# =============================================================================
# Plotting functions for filtering
# =============================================================================


def find_fit(cohort):
    """Find the first fit trajectory of the cohort"""
    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                return traj


def find_neutral(cohort):
    """Find the first fit trajectory of the cohort"""
    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is False:
                return traj


def gene_bar(cohort, split=True, relative=False, color=colors[0]):
    """Bar plot with counts of filtered mutations by gene."""

    if split is True:
        """Bar plot with counts of filtered mutations by gene."""

        # Create a dictionary with all filtered genes
        gene_list = []
        for traj in cohort:
            gene_list.append(traj.mutation.split()[0])

        gene_dict = {element: 0 for element in set(gene_list)}
        gene_dict_increasing = {element: 0 for element in set(gene_list)}
        gene_dict_decreasing = {element: 0 for element in set(gene_list)}

        # update the counts for each gene
        for traj in cohort:
            gene_dict[traj.mutation.split()[0]] = \
                gene_dict[traj.mutation.split()[0]] + 1

            # Separate counts for increasing and decreasing trajectories
            gradient = traj.y[-1]-traj.y[0] > 0
            if gradient is True:
                gene_dict_increasing[traj.mutation.split()[0]] = \
                    gene_dict_increasing[traj.mutation.split()[0]] + 1
            else:
                gene_dict_decreasing[traj.mutation.split()[0]] = \
                    gene_dict_decreasing[traj.mutation.split()[0]] + 1

        # sort dictionary in descending order
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: item[1], reverse=True))

        # extract ordered list of genes and values
        genes = [key for key in gene_dict.keys()]
        values_increasing = [gene_dict_increasing[gene] for gene in genes]
        if relative is False:
            values_decreasing = [gene_dict_decreasing[gene] for gene in genes]
        else:
            values_decreasing = [-gene_dict_decreasing[gene] for gene in genes]

        total_inc = sum(values_increasing)
        total_dec = np.abs(sum(values_decreasing))

        # Bar plot
        fig = go.Figure()
        fig.add_trace(
            go.Bar(x=genes, y=values_increasing,
                   name=f'{total_inc} Increasing',
                   marker_color=color))
        fig.add_trace(
            go.Bar(x=genes, y=values_decreasing,
                   name=f'{total_dec} Decreasing',
                   marker_color=color,
                   opacity=0.3))

        fig.update_layout(
            title='Gene distribution of filtered mutations',
            barmode='relative',
            template="simple_white",
            yaxis_title='Trajectory counts',
            xaxis_tickangle=-45)

    else:                   # unsplit barplot

        # Create a dictionary with all filtered genes
        gene_list = []
        for traj in cohort:
            gene_list.append(traj.mutation.split()[0])
        gene_dict = {element: 0 for element in set(gene_list)}

        # update the counts for each gene
        for traj in cohort:
            gene_dict[traj.mutation.split()[0]] = \
                gene_dict[traj.mutation.split()[0]] + 1
        # sort dictionary in descending order
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: item[1], reverse=True))
        genes = [key for key in gene_dict.keys()]
        values = [value for value in gene_dict.values()]

        # Bar plot
        fig = go.Figure()
        fig.add_trace(go.Bar(x=genes, y=values))
        fig.update_layout(title='Gene distribution of filtered mutations',
                          template="simple_white",
                          xaxis_tickangle=-45)
    # Labels position
    fig.update_layout(legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95,
        ))
    fig.update_layout(margin=dict(pad=4))

    return fig, gene_dict


def gradient_threshold_plot(cohort):
    """ Plot the selection of filtered variants.
    threshold: Default = 0.02.
    x_axis: VAF.
    y_axis: gradient.
    color: filter attribute.
    """


    # Find fit and netural trajectories to create group
    fit_traj = find_fit(cohort)
    artifact_traj = find_neutral(cohort)

    def legend_group(filter):
        """Legend group assigns a legend group based on filter attribute"""
        if filter is True:
            return 'fit mutation'
        else:
            return 'artifact'

    def color_group(filter):
        """ returns a different color based on filter attribute"""
        return colors[color_dict[legend_group(filter)]]

    
    def opacity_group(filter):
        """ returns a different color based on filter attribute"""
        return opacity_dict[legend_group(filter)]

    # Create figure with inset
    fig = go.Figure()
    # Create trajectories for legend groups
    # Fit trajectories
    fig.add_trace(
        go.Scatter(x=[fit_traj.data.AF.max()],
                    y=[fit_traj.gradient],
                    name='Fit mutation',
                    mode='markers',
                    marker_color=color_group(fit_traj.filter),
                    marker_size=5,
                    legendgroup=legend_group(fit_traj.filter),
                    showlegend=True))


    # Neutral trajectories
    fig.add_trace(
        go.Scatter(x=[artifact_traj.data.AF.max()],
                    y=[artifact_traj.gradient],
                    name='<2%',
                    mode='markers',
                    marker_color=color_group(artifact_traj.filter),
                    marker_size=5,
                    showlegend=True))

   # Add vertical delimitation of threshold
    fig.add_vrect(x0=0.01, x1=0.02, line_width=0, fillcolor="red", opacity=0.1)

    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is False:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=5,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)


    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=5,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

    fig.update_layout(
        title=None,
        template='simple_white',
        xaxis=dict(title='VAF', type='log', tickvals=[0.01, 0.02, 0.05, 0.1, 0.2]),
        yaxis=dict(title='Gradient',
                   range=[-0.02, 0.02]
                  )
    )

    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)


    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0)
    )

    return fig


def CHIP_plot(cohort, threshold=0.02):
    """ Plot the selection of filtered variants.
    threshold: Default = 0.02.
    x_axis: VAF.
    y_axis: gradient.
    color: filter attribute.
    """

    fit_traj = find_fit(cohort)
    neutral_traj = find_neutral(cohort)

    def legend_group(filter):
        if filter is True:
            return 'Fit'
        else:
            return 'Neutral'

    def color_group(filter):
        if filter is True:
            return colors[0]
        else:
            return colors[4]

    fig = go.Figure()
    # Create trajectories for legend groups
    # Fit trajectories
    if fit_traj is not None:
        fig.add_trace(
            go.Scatter(x=[fit_traj.data.AF.max()],
                       y=[fit_traj.gradient],
                       name='CHIP mutations',
                       mode='markers',
                       marker=dict(color=color_group(fit_traj.filter),
                                   size=5),
                       legendgroup=legend_group(fit_traj.filter),
                       showlegend=True))

    # Neutral trajectories
    if neutral_traj is not None:
        fig.add_trace(
            go.Scatter(x=[neutral_traj.data.AF.max()],
                       y=[neutral_traj.gradient],
                       name='non-CHIP mutations',
                       mode='markers',
                       marker=dict(color=color_group(neutral_traj.filter),
                                   size=5),
                       legendgroup=legend_group(neutral_traj.filter),
                       showlegend=True))

    for part in cohort:
        for traj in part.trajectories:
            fig.add_trace(
                go.Scatter(x=[traj.data.AF.max()],
                           y=[traj.gradient],
                           mode='markers',
                           marker=dict(color=color_group(traj.filter),
                                       size=5),
                           legendgroup=legend_group(traj.filter),
                           showlegend=False))
    # Add vertical delimitation of threshold
    fig.add_trace(
        go.Scatter(x=[0.02, 0.02], y=[-0.05, 0.1],
                   name=f'{threshold} VAF threshold',
                   mode='lines',
                   line=dict(color=colors[2], width=3, dash='dash')))
    fig.update_layout(title=('Filtering trajectories '
                             'achieving VAF > 0.02 over timespan'))
    fig.update_xaxes(title='VAF', range=[0, 0.35])
    fig.update_yaxes(title='Gradient', range=[-0.03, 0.08])

    fig.update_layout(template='plotly_white')

    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.95,
        xanchor="right",
        x=0.95,
    ))

    return fig


def CHIP_plot_inset(cohort, threshold=0.02,
                    xrange_2=[0, 0.03], yrange_2=[-0.01, 0.01]):
    """ Plot the selection of filtered variants.
    threshold: Default = 0.02.
    x_axis: VAF.
    y_axis: gradient.
    color: filter attribute.
    """

    # Find fit and netural trajectories to create group
    fit_traj = find_fit(cohort)
    neutral_traj = find_neutral(cohort)

    def legend_group(filter):
        """Legend group assigns a legend group based on filter attribute"""
        if filter is True:
            return 'Fit'
        else:
            return 'Neutral'

    def color_group(filter):
        """ returns a different color based on filter attribute"""
        if filter is True:
            return colors[0]
        else:
            return colors[4]

    # Create figure with inset
    fig = make_subplots(insets=[{'cell': (1, 1), 'l': 0.7, 'b': 0.3}])

    # Create trajectories for legend groups
    # Fit trajectories
    if fit_traj is not None:
        fig.add_trace(
            go.Scatter(x=[fit_traj.data.AF.max()],
                       y=[fit_traj.gradient],
                       name='CHIP mutations',
                       mode='markers',
                       marker=dict(color=color_group(fit_traj.filter),
                                   size=5),
                       legendgroup=legend_group(fit_traj.filter),
                       showlegend=True))

    # Neutral trajectories
    if neutral_traj is not None:
        fig.add_trace(
            go.Scatter(x=[neutral_traj.data.AF.max()],
                       y=[neutral_traj.gradient],
                       name='non-CHIP mutations',
                       mode='markers',
                       marker=dict(color=color_group(neutral_traj.filter),
                                   size=5),
                       legendgroup=legend_group(neutral_traj.filter),
                       showlegend=True))

    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is False:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker=dict(color=color_group(traj.filter),
                                            size=5),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

                # add inset scatter plot
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                xaxis='x2', yaxis='y2',
                                mode='markers',
                                marker=dict(color=color_group(traj.filter),
                                            size=3),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker=dict(color=color_group(traj.filter),
                                            size=5),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

                # add inset scatter plot
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                xaxis='x2', yaxis='y2',
                                mode='markers',
                                marker=dict(color=color_group(traj.filter),
                                            size=3),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

    # Add vertical delimitation of threshold
    fig.add_trace(
        go.Scatter(x=[0.02, 0.02], y=[-0.05, 0.1],
                   name=f'{threshold} VAF threshold',
                   mode='lines',
                   line=dict(color=colors[2], width=3, dash='dash')))

    # Add rectangle in zoomed area
    fig.add_shape(type="rect",
                  x0=xrange_2[0], y0=yrange_2[0],
                  x1=xrange_2[1], y1=yrange_2[1],
                  line=dict(color="Black", width=1))

    fig.update_layout(
        title=None,
        template='plotly_white',
        xaxis=dict(
            title='VAF',
            range=[-0.01, 0.35]),
        yaxis=dict(
            title='Gradient',
            range=[-0.017, 0.08]),
        yaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=yrange_2,
            domain=[0.55, 1],
            showline=None,
            linewidth=1,
            linecolor='black',
            mirror=True),
        xaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=xrange_2,
            domain=[0.6, 1],
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True)
    )

    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0
    ))

    return fig



def neutral_plot(cohort, neutral, mean_regr, var_regr, n_std=2):
    # Compute the standard deviation of absolute regularized total gradients
    fig = go.Figure()

    # Add traces for labeling
    fig.add_trace(go.Scatter(x=[0], y=[0],
                             line=dict(color=colors[0]),
                             mode='markers',
                             legendgroup="fit",
                             name="Fit variants"
                             ))

    fig.add_trace(go.Scatter(x=[0], y=[0],
                             line=dict(color=colors[3]),
                             mode='markers',
                             legendgroup="non-fit",
                             name="Neutral variants"
                             ))

    fig.add_trace(go.Scatter(x=[0],
                             y=[0],
                             line=dict(color='Orange'),
                             mode='markers',
                             legendgroup="synonymous",
                             name="Synonymous variants"
                             ))

    # Add traces
    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                fig.add_trace(
                    go.Scatter(x=[traj.data.AF.iloc[0]],
                               y=[traj.gradient],
                               line=dict(color=colors[0]),
                               legendgroup="fit",
                               showlegend=False,
                               hovertemplate=f'Mutation: {traj.mutation}'
                               ))
            else:
                fig.add_trace(
                    go.Scatter(x=[traj.data.AF.iloc[0]],
                               y=[traj.gradient],
                               line=dict(color=colors[3]),
                               legendgroup="non-fit",
                               showlegend=False,
                               hovertemplate=f'Mutation: {traj.mutation}'
                               ))
    for part in neutral:
        for traj in part.trajectories:
            fig.add_trace(
                go.Scatter(x=[traj.data.AF.iloc[0]],
                           y=[traj.gradient],
                           line=dict(color=colors[4]),
                           marker=dict(size=3),
                           legendgroup="synonymous",
                           showlegend=False,
                           hovertemplate=f'Mutation: {traj.mutation}'
                           ))

    x = np.linspace(0, 0.5, 1000)
    y_std = mean_regr.predict(np.c_[x]) \
        + n_std*np.sqrt(var_regr.predict(np.c_[x]))

    y_mean = mean_regr.predict(np.c_[x])

    fig.add_trace(
        go.Scatter(x=x, y=y_std,
                   name='Neutral growth filter',
                   marker=dict(color=colors[1])))
    fig.add_trace(
        go.Scatter(x=x, y=y_mean,
                   name='Neutral growth mean',
                   marker=dict(color=colors[2])))

    fig.update_layout(title=('Filtering according to '
                             'neutral growth distribution'))
    fig.update_xaxes(title='VAF', range=[0, 0.35])
    fig.update_yaxes(title='Gradient', range=[-0.03, 0.08])

    fig.update_layout(template='plotly_white')
    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0
    ))

    return fig


def neutral_plot_inset(cohort, neutral, mean_regr, var_regr,
                       n_std=2, xrange_2=[0, 0.03], yrange_2=[-0.01, 0.01]):

    # Find fit and netural trajectories to create group
    fit_traj = find_fit(cohort)
    neutral_traj = find_neutral(cohort)

    def legend_group(filter):
        """Legend group assigns a legend group based on filter attribute"""
        if filter is True:
            return 'Fit'
        else:
            return 'Neutral'

    def color_group(filter):
        """ returns a different color based on filter attribute"""
        if filter is True:
            return colors[0]
        else:
            return colors[3]

    # Create figure with inset
    fig = make_subplots(insets=[{'cell': (1, 1), 'l': 0.7, 'b': 0.3}])

    # Create trajectories for legend groups
    # Fit trajectories
    if fit_traj is not None:
        fig.add_trace(
            go.Scatter(x=[fit_traj.data.AF.iloc[0]],
                       y=[fit_traj.gradient],
                       name='Fit',
                       mode='markers',
                       marker=dict(color=color_group(fit_traj.filter),
                                   size=5),
                       legendgroup=legend_group(fit_traj.filter),
                       showlegend=True))

    # Neutral trajectories label
    if neutral_traj is not None:
        fig.add_trace(
            go.Scatter(x=[neutral_traj.data.AF.iloc[0]],
                       y=[neutral_traj.gradient],
                       name='Neutral',
                       mode='markers',
                       marker=dict(color=color_group(neutral_traj.filter),
                                   size=5),
                       legendgroup='synonymous',
                       showlegend=True))

    # Synonymous trajectories label:
    synonymous = neutral[0].trajectories[0]
    fig.add_trace(
        go.Scatter(x=[synonymous.data.AF.iloc[0]],
                   y=[synonymous.gradient],
                   name='Synonymous',
                   mode='markers',
                   marker=dict(color='Orange',
                               size=5),
                   legendgroup='synonymous',
                   showlegend=True))

    for part in cohort:
        for traj in part.trajectories:
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            mode='markers',
                            marker=dict(color=color_group(traj.filter),
                                        size=5),
                            legendgroup=legend_group(traj.filter),
                            showlegend=False)

            # add inset scatter plot
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            xaxis='x2', yaxis='y2',
                            mode='markers',
                            marker=dict(color=color_group(traj.filter),
                                        size=3),
                            legendgroup=legend_group(traj.filter),
                            showlegend=False)
    for part in neutral:
        for traj in part.trajectories:
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            mode='markers',
                            marker=dict(color='Orange',
                                        size=3),
                            legendgroup='synonymous',
                            showlegend=False)

            # add inset scatter plot
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            xaxis='x2', yaxis='y2',
                            mode='markers',
                            marker=dict(color='Orange',
                                        size=3),
                            legendgroup='synonymous',
                            showlegend=False)

    # Add mean and filter lines
    x = np.linspace(0, 0.5, 1000)
    y_std = mean_regr.predict(np.c_[x]) \
        + n_std*np.sqrt(var_regr.predict(np.c_[x]))
    y_mean = mean_regr.predict(np.c_[x])
    fig.add_trace(
        go.Scatter(x=x, y=y_std,
                   name='Neutral growth filter',
                   xaxis='x2', yaxis='y2',
                   legendgroup='filter',
                   marker=dict(color=colors[1])))
    fig.add_trace(
        go.Scatter(x=x, y=y_mean,
                   name='Neutral growth mean',
                   xaxis='x2', yaxis='y2',
                   legendgroup='mean',
                   marker=dict(color=colors[2])))
    # Add mean and filter lines to inset
    fig.add_trace(
        go.Scatter(x=x, y=y_std,
                   name='Neutral growth filter',
                   legendgroup='filter',
                   showlegend=False,
                   marker=dict(color=colors[1])))
    fig.add_trace(
        go.Scatter(x=x, y=y_mean,
                   name='Neutral growth mean',
                   legendgroup='mean',
                   showlegend=False,
                   marker=dict(color=colors[2])))

    # Add rectangle in zoomed area
    fig.add_shape(type="rect",
                  x0=xrange_2[0], y0=yrange_2[0],
                  x1=xrange_2[1], y1=yrange_2[1],
                  line=dict(color="Black", width=1))

    # Update inset layout
    fig.update_layout(
        title=None,
        template='plotly_white',
        xaxis=dict(
            title='VAF',
            range=[0, 0.35]),
        yaxis=dict(
            title='Gradient',
            range=[-0.017, 0.08]),
        yaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=yrange_2,
            domain=[0.55, 1],
            showline=None,
            linewidth=1,
            linecolor='black',
            mirror=True),
        xaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=xrange_2,
            domain=[0.6, 1],
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True)
    )

    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0
    ))

    return fig


def stack_plot(melt_syn, bin_center,
               mean, mean_regr, mean_quant,
               variance, var_regr, var_quant):
    # Set some overall variables
    mean_color = colors[2]
    var_color = colors[1]
    bin_color = 'Grey'

    # NGF stack plot
    # Start subplot figure with shared xaxis
    NGF_stack = make_subplots(rows=3, cols=1,
                              row_heights=[0.4, 0.4, 0.2],
                              shared_xaxes=True,
                              vertical_spacing=0.05)

    # Row 1 variance plot
    # add trace with data
    NGF_stack.add_trace(
        go.Scatter(
            x=bin_center,
            y=variance,
            mode='markers',
            marker_color=bin_color,
            showlegend=False,
            error_y=dict(
                type='data',
                symmetric=False,
                array=var_quant[0],
                arrayminus=var_quant[1])
            ),
        row=1, col=1)
    # add trace with linear fit
    x = [0, max(bin_center)]
    NGF_stack.add_trace(
        go.Scatter(x=x,
                   y=var_regr.predict(np.c_[x]),
                   mode='lines',
                   marker_color=var_color,
                   showlegend=False),
        row=1, col=1)

    # Row 2 mean plot
    # add trace with data
    NGF_stack.add_trace(
        go.Scatter(
            x=bin_center,
            y=mean,
            mode='markers',
            marker_color=bin_color,
            showlegend=False,
            error_y=dict(
                type='data',
                symmetric=False,
                array=mean_quant[0],
                arrayminus=mean_quant[1])
            ),
        row=2, col=1)

    # add trace with linear fit
    x = [0, max(bin_center)]
    NGF_stack.add_trace(
        go.Scatter(
            x=x,
            y=mean_regr.predict(np.c_[x]),
            mode='lines',
            marker_color=mean_color,
            showlegend=False),
        row=2, col=1)

    # Compute bins for barplot
    bar_center = []
    bar_width = []
    bar_size = []
    for bin in melt_syn.bin.unique():
        bar_size.append(
            len(melt_syn[melt_syn['bin'] == bin]))
        bar_center.append((bin.left+bin.right)/2)
        bar_width.append(bin.right-bin.left)

    # Row 3 Plot bin histogram
    NGF_stack.add_trace(
        go.Bar(
            x=bar_center,
            y=bar_size,
            marker_color='Orange',
            width=bar_width,
            showlegend=False),
        row=3, col=1)

    NGF_stack.update_yaxes(row=1, col=1,
                           title_text='Gradient Variance',
                           range=[0, max(variance)])
    NGF_stack.update_yaxes(title_text='Gradient', row=2, col=1)
    NGF_stack.update_yaxes(title_text='Bin counts',
                           tickmode='linear', dtick=400, row=3, col=1)
    NGF_stack.update_xaxes(title_text="VAF", row=3, col=1)
    NGF_stack.update_layout(template='plotly_white')
    NGF_stack.update_xaxes(ticks='outside', row=3, col=1)
    NGF_stack.update_yaxes(ticks='outside')
    NGF_stack.update_layout(
        template='plotly_white',
        legend=dict(orientation="h",
                    yanchor="bottom", y=1.02,
                    xanchor="right", x=1))
    return NGF_stack


def check_plot(cohort, mean_regr, var_regr, mean_color, var_color, n_std):
    """Plot showing the filter acting on a cohort.
    Parameters:
    cohort: List of participant class objects.
    mean_regr: Scikit LR model. Linear regression model of DIF gradients
               as a function of VAF.
    var_regr: Scikit LR model. Linear regression model of DIF gradient variance
               as a function of VAF.
    n_std: int. Standard deviations from the expected mean of DIF fluctuations
           used to filter variants.
    Returns:
    fig: Plotly figure."""

    # Create a dataframe with the regularized gradient
    # between first and last timepoint
    total_gradient = pd.DataFrame(columns=['AF', 'gradient'])
    for part in cohort:
        for traj in part.trajectories:
            if len(traj.data) > 1:
                total_gradient = (
                    total_gradient.append({'AF': traj.data.iloc[0]['AF'],
                                           'gradient': traj.gradient},
                                          ignore_index=True))

    # Compute the percentage of succesfully filtered trajectories:
    x = total_gradient['AF']
    total_gradient['filter'] = total_gradient['gradient'] \
        - (mean_regr.predict(np.c_[x])
           + n_std*np.sqrt(var_regr.predict(np.c_[x])))
    percentage = 100 * len(total_gradient[total_gradient['filter'] < 0]) \
        / (len(total_gradient)-1)

    # Plot filter
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=total_gradient['AF'],
                   y=total_gradient['gradient'],
                   name='Synonymous mutations',
                   mode='markers'))

    # Confidence line
    x = np.linspace(0, max(total_gradient['AF']), 1000)
    y_std = mean_regr.predict(np.c_[x]) \
        + n_std*np.sqrt(var_regr.predict(np.c_[x]))
    y_mean = mean_regr.predict(np.c_[x])

    fig.add_trace(
        go.Scatter(x=x, y=y_mean,
                   name='Mean linear regression',
                   marker_color=mean_color))

    fig.add_trace(
        go.Scatter(x=x, y=y_std,
                   name='Neutral growth filter',
                   marker_color=var_color))

    fig.update_layout(title=(f'Neutral growth filter <br>Filters '
                             f'{round(percentage,1)}% of all '
                             f'synonymous mutations'),
                      template='plotly_white')
    fig.update_xaxes(title='VAF')
    fig.update_yaxes(title='Regularized gradient')

    return fig


def participant_filter(part):
    """ Plot a participant's profile using different colors for
    fit and non-fit mutations"""

    with open('../Resources/var_dict.json') as json_file:
        var_dict = json.load(json_file)

    fig = go.Figure()

    # Empty trajectories for labeling fit and non-fit trajectories
    fig.add_trace(
        go.Scatter(x=[79], y=[0],
                   name='fit mutations',
                   mode='lines',
                   marker=dict(color='Black'),
                   legendgroup='non-fit',
                   showlegend=True))

    fig.add_trace(
        go.Scatter(x=[79], y=[0],
                   name='non-fit mutations',
                   mode='lines',
                   marker=dict(color='lightgrey'),
                   legendgroup='non-fit',
                   showlegend=True))

    for traj in part.trajectories:
        # Plot trajectories  using different colors according to traj.filter
        if traj.filter is False:
            fig.add_trace(
                go.Scatter(x=traj.data.age, y=traj.data.AF,
                           name='non-fit',
                           mode='lines',
                           marker=dict(color='lightgrey'),
                           legendgroup='non-fit',
                           showlegend=False))

    for traj in part.trajectories:
        if traj.filter is True:
            fig.add_trace(
                go.Scatter(x=traj.data.age, y=traj.data.AF,
                           name=traj.mutation,
                           mode='lines',
                           marker=dict(color=var_dict[traj.variant_class]),
                           legendgroup='fit',
                           showlegend=True))

    # Update layout
    fig.update_layout(title=f'Participant {part.id}',
                      xaxis_title='Age',
                      yaxis_title='VAF',
                      legend=dict(
                        y=0.95,
                        x=0.1))

    return fig




#endregion
# =============================================================================
#region LiFT
# =============================================================================

color_dict = {'artifact': 1, 'below_2%': 1, 'fit mutation':0}
opacity_dict = {'artifact': 0.3, 'below_2%': 0.3, 'fit mutation':1}

def plot_cohort_LiFT(cohort, filter_type = 'LiFT'):
    fig = go.Figure()

    if filter_type == "2%":
        model_types = ['below_2%', 'fit mutation']
    else:
        model_types = ['artifact', 'fit mutation']

    for optimal_model in model_types:
        fig.add_trace(
            go.Scatter(x=[70], y=[0],
                        name=optimal_model,
                        marker_color=colors[color_dict[optimal_model]],
                        opacity=opacity_dict[optimal_model]))

    for optimal_model in model_types:
        marker_color = colors[color_dict[optimal_model]]
        opacity = opacity_dict[optimal_model]
        for mutation in cohort:
            if mutation.optimal_model == optimal_model:
                if type(mutation.data) is not list:
                    fig.add_trace(
                        go.Scatter(
                            x=mutation.data.age,
                            y=mutation.data.AF,
                            showlegend=False,
                            marker_color=marker_color,
                            opacity=opacity
                        )
                    )   
                else:           
                    for traj in mutation.data:
                        fig.add_trace(
                                go.Scatter(
                                    x=traj.age,
                                    y=traj.AF,
                                    showlegend=False,
                                    marker_color=marker_color,
                                    opacity=opacity
                                )
                            )
    fig.update_layout(xaxis_title='Age',
                        yaxis_title='VAF')

    fig.update_layout(margin=dict(pad=4))
    return fig


def bayes_factor_threshold_effect(mutation_list):
    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    k_range = np.arange(1,20)
    k_mutation_counts = []
    max_counts = 0
    max_percentage = 0
    for k in k_range:
        counts = len([mutation for mutation in mutation_list
                    if (mutation.bd_prob/mutation.betabinom_artifact_prob > k
                        )])
        k_mutation_counts.append(counts)
    k_mutation_percentage = [100*count/len(mutation_list) for count in k_mutation_counts]

    max_counts = max(max_counts, max(k_mutation_counts))
    max_percentage = max(max_percentage, max(k_mutation_percentage))

    fig.add_trace(
        go.Scatter(x=k_range, y=k_mutation_percentage, marker_color=colors[0], name='Beta-binomial'), secondary_y=False)

    fig.add_trace(
        go.Scatter(x=k_range, y=k_mutation_counts, marker_color=colors[0], name='Beta-binomial', visible=False, showlegend=False), secondary_y=True)
    k_mutation_counts = []

    for k in k_range:
        counts = len([mutation for mutation in mutation_list
                    if (mutation.bd_prob/mutation.binom_artifact_prob > k
                        )])
        k_mutation_counts.append(counts)
    k_mutation_percentage = [100*count/len(mutation_list) for count in k_mutation_counts]
    max_counts = max(max_counts, max(k_mutation_counts))
    max_percentage = max(max_percentage, max(k_mutation_percentage))

    for mutation in mutation_list:
        mutation.compute_optimal_model(binomial=False, bayes_factor_threshold=4)


    fig.add_trace(
        go.Scatter(x=k_range, y=k_mutation_percentage, marker_color=colors[1], name='Binomial'), secondary_y=False)
    fig.add_trace(
        go.Scatter(x=k_range, y=k_mutation_counts, 
                   marker_color=colors[1],
                   visible= False, 
                   showlegend=False,
                   name='Binomial'),
            secondary_y=True)
    
    fig.update_layout(
        title="Bayes Factor filter",
        xaxis_title='Bayes Factor',
        margin=dict(pad=4),
        legend=dict(
            xanchor="right",
            x=0.9))

    fig.update_yaxes(title='% of total mutations', range = [-0.05, max_percentage], secondary_y=False)
    fig.update_yaxes(title='Mutation counts', range = [-0.05, max_counts], secondary_y=True)

    return fig


class pvalue_mutation ():
    def __init__(self, id, outlier_pvalue=None, mean_AF=None, optimal_model=None):
        self.id = id
        self.outlier_pvalue = outlier_pvalue
        self.mean_AF = mean_AF
        self.optiaml_model = optimal_model 

def outlier_pvalue_comparison(cohort, fit_mutation_ids):
    """Plot comparing LiFT results with outlier_pvalue."""

    # Create list of pvalue_mutation class objects
    # Remark: Not every mutation has a pvalue.
    pvalue_mutation_list = []
    ids = []
    for part in cohort:
        for traj in part.trajectories:
            if traj.AF_outlier_pvalue is not None:
                if traj.mutation not in  ids:
                    ids.append(traj.mutation)
                    pvalue_mutation_list.append(
                        pvalue_mutation(id=traj.mutation,
                                outlier_pvalue=[],
                                mean_AF=[])
                    )
                for mutation in pvalue_mutation_list:
                    if mutation.id == traj.mutation:
                        mutation.outlier_pvalue.append(traj.AF_outlier_pvalue)
                        mutation.mean_AF.append(np.mean(traj.data.AF))

    for mutation in pvalue_mutation_list:
        mutation.outlier_pvalue = np.mean(mutation.outlier_pvalue)
        mutation.mean_AF = np.mean(mutation.mean_AF)
        
        if mutation.id in fit_mutation_ids:
            mutation.optimal_model = 'fit mutation'
        else:
            mutation.optimal_model = 'artifact' 


    # Create figure
    fig = go.Figure()

    for optimal_model in ['artifact', 'fit mutation']:
        fig.add_trace(
            go.Scatter(x=[0], y=[0.01],
                       name=optimal_model,
                       marker_color=colors[color_dict[optimal_model]],
                       opacity=opacity_dict[optimal_model]
                    )
        )

    for optimal_model in ['artifact', 'fit mutation']:
        marker_color = colors[color_dict[optimal_model]]
        opacity = opacity_dict[optimal_model]
        for mutation in pvalue_mutation_list:
            if mutation.optimal_model == optimal_model:
                fig.add_trace(
                    go.Scatter(y=[mutation.mean_AF],
                            x=[mutation.outlier_pvalue],
                            showlegend=False,
                            marker_color = marker_color,
                            hovertext = mutation.id,
                            opacity=opacity))
    fig.update_layout(title='LiFT vs outlier pvalue',
                      xaxis_title='outlier p value',
                      yaxis_title='mean VAF')
    fig.update_layout(margin=dict(pad=4))

    return fig, pvalue_mutation_list

class mdaf_mutation ():
    def __init__(self, id, mdaf=None, mean_AF=None, optimal_model=None):
        self.id = id
        self.mdaf = mdaf
        self.mean_AF = mean_AF
        self.optiaml_model = optimal_model 

def mdaf_comparison(cohort, fit_mutation_ids):
    """Plot comparing LiFT results with outlier_pvalue."""
    # Create list of all mutations by cohort with repeats
    mdaf_mutation_list = []
    ids = []
    for part in cohort:
        for traj in part.trajectories:
            if traj.mdaf is not None:
                if traj.mutation not in  ids:
                    ids.append(traj.mutation)
                    mdaf_mutation_list.append(
                        mdaf_mutation(id=traj.mutation,
                                mdaf=[],
                                mean_AF=[])
                    )
                for mutation in mdaf_mutation_list:
                    if mutation.id == traj.mutation:
                        mutation.mdaf.append(traj.mdaf)
                        mutation.mean_AF.append(np.mean(traj.data.AF))

    for mutation in mdaf_mutation_list:
        mutation.mdaf = np.mean(mutation.mdaf)
        mutation.mean_AF = np.mean(mutation.mean_AF)
        
        if mutation.id in fit_mutation_ids:
            mutation.optimal_model = 'fit mutation'
        else:
            mutation.optimal_model = 'artifact' 


    # Plot figure
    fig = go.Figure()

    for optimal_model in ['artifact', 'fit mutation']:
        fig.add_trace(
            go.Scatter(x=[0.02], y=[0.01],
                        name=optimal_model,
                        marker_color=colors[color_dict[optimal_model]],
                        opacity=opacity_dict[optimal_model]
                    )
        )

    for optimal_model in ['artifact', 'fit mutation']:
        marker_color = colors[color_dict[optimal_model]]
        opacity = opacity_dict[optimal_model]
        for mutation in mdaf_mutation_list:
            if mutation.optimal_model == optimal_model:
                fig.add_trace(
                    go.Scatter(y=[mutation.mean_AF],
                            x=[mutation.mdaf],
                            showlegend=False,
                            marker_color = marker_color,
                            hovertext = mutation.id,
                            opacity=opacity))

    fig.update_layout(title='LiFT vs MDAF',
                      xaxis_title='MDAF',
                      yaxis_title='mean VAF')
    fig.update_layout(margin=dict(pad=4))

    return fig, mdaf_mutation_list



def LiFT_plot_inset(cohort,
                    xrange_2=[0.01, 0.03],yrange_2=[-0.01, 0.01]):
    """ Plot the selection of filtered variants.
    threshold: Default = 0.02.
    x_axis: VAF.
    y_axis: gradient.
    color: filter attribute.
    """


    # Find fit and netural trajectories to create group
    fit_traj = find_fit(cohort)
    artifact_traj = find_neutral(cohort)

    def legend_group(filter):
        """Legend group assigns a legend group based on filter attribute"""
        if filter is True:
            return 'fit mutation'
        else:
            return 'artifact'

    def color_group(filter):
        """ returns a different color based on filter attribute"""
        return colors[color_dict[legend_group(filter)]]

    
    def opacity_group(filter):
        """ returns a different color based on filter attribute"""
        return opacity_dict[legend_group(filter)]


    # Create figure with inset
    fig = make_subplots(insets=[{'cell': (1, 1), 'l': 0.7, 'b': 0.3}])

    # Create trajectories for legend groups
    # Fit trajectories
    fig.add_trace(
        go.Scatter(x=[fit_traj.data.AF.max()],
                    y=[fit_traj.gradient],
                    name='LiFT-positive',
                    mode='markers',
                    marker_color=color_group(fit_traj.filter),
                    marker_size=5,
                    marker_opacity=opacity_group(fit_traj.filter),
                    legendgroup=legend_group(fit_traj.filter),
                    showlegend=True))

    # Neutral trajectories
    fig.add_trace(
        go.Scatter(x=[artifact_traj.data.AF.max()],
                    y=[artifact_traj.gradient],
                    name='LiFT-negative',
                    mode='markers',
                    marker_color=color_group(artifact_traj.filter),
                    marker_size=5,
                    marker_opacity=opacity_group(artifact_traj.filter),
                    showlegend=True))

    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is False:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=5,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

                # add inset scatter plot
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                xaxis='x2', yaxis='y2',
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=3,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=5,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

                # add inset scatter plot
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                xaxis='x2', yaxis='y2',
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=3,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

    # Add vertical delimitation of threshold
    fig.add_trace(
        go.Scatter(x=[0.02, 0.02], y=[-0.05, 0.1],
                   name=f'0.02 VAF threshold',
                   mode='lines',
                   line=dict(color=colors[2], width=3, dash='dash')))

    # Add rectangle in zoomed area
    fig.add_shape(type="rect",
                  x0=xrange_2[0], y0=yrange_2[0],
                  x1=xrange_2[1], y1=yrange_2[1],
                  line=dict(color="Black", width=1))

    fig.update_layout(
        title=None,
        template='plotly_white',
        xaxis=dict(
            title='VAF',
            range=[0.01, 0.35]),
        yaxis=dict(type='log',
            title='Gradient',
            range=[-100, 0.08]),
        yaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=yrange_2,
            domain=[0.55, 1],
            showline=None,
            linewidth=1,
            linecolor='black',
            mirror=True),
        xaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=xrange_2,
            domain=[0.6, 1],
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True)
    )

    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0
    ))

    return fig




#endregion
# =============================================================================
#region Clonal_model
# =============================================================================

# Load variant color dictionary
with open("../Resources/var_dict.json") as json_file:
    var_dict = json.load(json_file)

default_color="Grey"

def gene_trajectories(model_list, gene, jitter=1, showlegend=False):
    """Plot all variant trajectories present in a given gene.
    Parameters:
    model_list: List. Collection of all fitted trajectories.
    gene: Float. Gene name.
    jitter: float. Amount of jitter used in fitness boxplot.
    Return:
    fig: Plotly figure. Side-by-side plot of variant trajectories and fitness
                        boxplot.
    """
    symbols = SymbolValidator().values

    fig = make_subplots(rows=1, cols=2,
                        column_widths=[0.7, 0.3],
                        subplot_titles=(f'{gene} trajectories',
                                        'Fitness distribution'))
    fitness_color = []
    fitness = []

    traj_list = [traj for traj in model_list if traj.gene == gene]
    for i, traj in enumerate(traj_list):
        x = traj.data.age
        y = traj.data.AF
        fig.add_trace(
            go.Scatter(x=x, y=y,
                       mode='markers',
                       marker_size=10,
                       marker_symbol=symbols[8*i + 1],
                       marker_color=var_dict[traj.variant_class],
                       name=traj.variant_class,
                       showlegend=showlegend),
            row=1, col=1)
        x_prediction = traj.deterministic_fit[0]
        y_prediction = traj.deterministic_fit[1]
        fig.add_trace(
            go.Scatter(x=x_prediction,
                       y=y_prediction, mode='lines',
                       marker_color=var_dict[traj.variant_class],
                       showlegend=showlegend,
                       name=traj.mutation),
            row=1, col=1)
        fitness.append(traj.fitness)
        fitness_color.append(var_dict[traj.variant_class])

    # box plot with fitness estimates
    # Frist create the box
    fig.add_trace(
        go.Box(y=fitness, boxpoints=False,
               name=gene,
               marker_color=default_color,
               showlegend=False),
        row=1, col=2)

    # Seccond create the dots adding jitter
    for i, item in enumerate(fitness):
        fig.add_trace(
            go.Box(y=[item],
                   marker_size=10,
                   pointpos=-2, boxpoints='all',
                   name=gene,
                   jitter=random.uniform(0, jitter),
                   marker_symbol=symbols[8*i + 1],
                   marker_color=fitness_color[i],
                   line=dict(color='rgba(0,0,0,0)'),
                   fillcolor='rgba(0,0,0,0)',
                   showlegend=False),
            row=1, col=2)

    fig.update_layout(legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0
        ))

    fig.update_xaxes(title_text='Age (in years)', linewidth=2, row=1, col=1)
    fig.update_yaxes(title_text='VAF', linewidth=2, row=1, col=1)
    fig.update_xaxes(showticklabels=False, linewidth=2, row=1, col=2)
    fig.update_yaxes(title_text='Fitness', linewidth=2, row=1, col=2)
    fig.update_layout(template='simple_white')


    return fig


def gene_fitness_plot(lbc, fitness_threshold=0.02):

    cohort = [traj for part in lbc for traj in part.trajectories if traj.fitness>fitness_threshold]
    
    # Create a list with the gene names
    gene_names=list()
    for traj in cohort:
        if traj.gene not in gene_names:
            gene_names.append(traj.gene)

    #Count the number of appearances of each gene in trajectories with over 0.02 of fitness
    gene_count=pd.DataFrame(columns=['gene','counts'])
    for gn in gene_names:
        count=0
        for traj in cohort:
            if traj.gene==gn:
                if traj.fitness>0.02:
                    count+=1
        gene_count=gene_count.append({'gene':gn,'counts':count},ignore_index=True)

    #Create a dataframe with the most relevant information of each variant

    fitness_data = pd.DataFrame(columns=['gene','fitness','clonal_size','cs','percentiles_minus','percentiles_plus','x_mpd','y_mpd','unique_name'])

    for traj in cohort:
        if traj.fitness>0.02:
            #Count the number of appearances of that gene
            g_c=list(gene_count.counts[gene_count.gene==traj.gene])[0]

            #Look if the variant belongs to a clone with one or more mutations
            if traj.clonal_population==1:
                numb_clones='Clone with 1 mutation'
                nc=1
            if traj.clonal_population>1:
                numb_clones='Clone with 2 or more mutations'
                nc=2
            
            fitness_data=fitness_data.append({'gene':traj.gene,'gene and count':traj.gene+' ('+str(g_c)+')','fitness':traj.fitness,'clonal_size':numb_clones,'cs':nc,'percentiles_minus':traj.fitness_quantiles[0],'percentiles_plus':traj.fitness_quantiles[1],
            'x_mpd':traj.fitness_posterior[0],'y_mpd':traj.fitness_posterior[1],
            'unique_name':traj.part_id + ' - ' + traj.mutation},ignore_index=True)

    #Sort the dataframe by criteria of fitness (and by percentile size in case of a tie in fitness)
    fitness_data=fitness_data.sort_values(by=['fitness','percentiles_plus'],ascending=False)
    fitness_data=fitness_data.reset_index()

    #Create a dataframe with the information about the median fitness for each gene
    median_fitness=pd.DataFrame(columns=['gene','mf','gene_and_count'])
    for gn in list(set(fitness_data.gene)):
        g_c=list(gene_count.counts[gene_count.gene==gn])[0]
        median_fitness=median_fitness.append({'gene':gn,'mf':np.median(fitness_data.fitness[fitness_data.gene==gn]),'gene_and_count':gn+' ('+str(g_c)+')'},ignore_index=True)

    #Arrange the dataframe from highest to lowest median fitness
    median_fitness=median_fitness.sort_values(by='mf',ascending=False)

    #The plot

    fig=go.Figure()

    #Parameters
    sep=3 #Separation between variants of the same gene (lower sep -> more separation)
    bxp_w=1.5 #Width of the boxplot lines
    palette= 2 #Select the color palette (1-> personalised for each gene; 2-> plotly palette)

    count=0
    gene_count=0
    x_markers=[0] #List with the positions for the gene labels

    #Look at the genes in order of median fitness
    for gn in median_fitness.gene:
        #Create a dataframe the rows with a mutation in that gene in the main datafrane (already ordered by fitness)
        d=fitness_data[fitness_data.gene==gn]

        #Select the color for that gene depending on the palette choice
        if palette==1:
            gene_color=gene_dict[gn]

        if palette==2:
            gene_color=colors[gene_count%10]
        
        gene_count+=1
        for m in range(len(d.gene)):
            #A boxplot is generated if there are 2 or more mutations of that gene
            if len(d.gene)==1:
                count+=1

            if m==0 and len(d.gene)>1:
                #Manual generation of the boxplot
                xm=[count/sep,(count+1)/sep] #boxplot width
                ym=[np.max(list(d.fitness)),np.min(list(d.fitness)),np.median(list(d.fitness))] #look for boxplot marks
                yq=list(np.percentile(list(d.fitness),q=[25,75])) #Find the 1st and 3rd quantiles

                #Draw all the lines of the boxplot
                fig.add_trace(
                    go.Scatter(x=[xm[0],xm[1],xm[1],xm[0],xm[0]],y=[yq[1],yq[1],yq[0],yq[0],yq[1]],fill='toself',mode='lines',line=dict(color=gene_color,width=bxp_w))
                )
                fig.add_trace(
                    go.Scatter(x=[xm[0],xm[1]],y=[ym[2],ym[2]],mode='lines',line=dict(color=gene_color,width=bxp_w))
                )
                fig.add_trace(
                    go.Scatter(x=[xm[0],xm[1]],y=[ym[0],ym[0]],mode='lines',line=dict(color=gene_color,width=bxp_w))
                )
                fig.add_trace(
                    go.Scatter(x=[xm[0],xm[1]],y=[ym[1],ym[1]],mode='lines',line=dict(color=gene_color,width=bxp_w))
                )
                fig.add_trace(
                    go.Scatter(x=[(xm[0]+xm[1])/2,(xm[0]+xm[1])/2],y=[ym[0],yq[1]],mode='lines',line=dict(color=gene_color,width=bxp_w))
                )
                fig.add_trace(
                    go.Scatter(x=[(xm[0]+xm[1])/2,(xm[0]+xm[1])/2],y=[ym[1],yq[0]],mode='lines',line=dict(color=gene_color,width=bxp_w))
                )
                count+=3
            
            #Look for the values of y that represent the 90% confidence interval limits
            ymin=np.argmin(abs(list(d.x_mpd)[m]-list(d.percentiles_minus)[m]))
            ymax=np.argmin(abs(list(d.x_mpd)[m]-list(d.percentiles_plus)[m]))

            #Normalise the height of the distributions so that they all have the same height at the maximum value
            xdist=list(d.y_mpd)[m]/(sep*np.max(list(d.y_mpd)[m]))

            #Create lists with the values for the distributions
            xvals=[count/sep]
            xvals.extend(-xdist[ymin:(ymax+1)]+count/sep)
            xvals.append(count/sep)
            yvals=[list(d.x_mpd)[m][ymin]]
            yvals.extend(list(d.x_mpd)[m][ymin:(ymax+1)])
            yvals.append(list(d.x_mpd)[m][ymax])

            #Plot the distribution (with opacity 50% if there is more than 1 mutation in that clone)
            fig.add_trace(
                go.Scatter(x=xvals,y=yvals,fill="toself",mode='lines',line=dict(color=gene_color,width=1),name=gn,opacity=1/list(d.cs)[m])
            )
            #Add the line of the fitness estimate for that mutation (dashed if there is more than 1 mutation in that clone)
            if list(d.cs)[m]==1:
                fig.add_trace(
                go.Scatter(x=[count/sep,-np.max(xdist)+count/sep],y=[list(d.fitness)[m],list(d.fitness)[m]],mode='lines',line=dict(dash='solid',color=gene_color),name=gn)
            )
            if list(d.cs)[m]==2:
                fig.add_trace(
                    go.Scatter(x=[count/sep,-np.max(xdist)+count/sep],y=[list(d.fitness)[m],list(d.fitness)[m]],mode='lines',line=dict(dash='dash',color=gene_color,width=0.5),name=gn,opacity=1/list(d.cs)[m])
                )
            count+=1
        count+=1
        x_markers.append(count/sep)

    #Add x-labels and set the axis range
    x_markers=list(np.array(x_markers)+1/(2*sep)) #Center the labels to the middle of the firsts distributions/boxplots
    fig.update_xaxes(tickangle=-45, tickvals=x_markers[:-1],ticktext=median_fitness.gene_and_count,range=[-1/sep,x_markers[-1]])
    fig.update_yaxes(title='Fitness')
    fig.update_yaxes(showgrid=True)
    fig.update_layout(showlegend=False)

    #Show the plot
    return fig



def gradient_cohort_plot(cohort,
                    xrange_2=[0.01, 0.03],yrange_2=[-0.01, 0.01]):
    """ Plot the selection of filtered variants.
    threshold: Default = 0.02.
    x_axis: VAF.
    y_axis: gradient.
    color: filter attribute.
    """


    # Find fit and netural trajectories to create group
    fit_traj = find_fit(cohort)
    artifact_traj = find_neutral(cohort)

    def legend_group(filter):
        """Legend group assigns a legend group based on filter attribute"""
        if filter is True:
            return 'fit mutation'
        else:
            return 'artifact'

    def color_group(filter):
        """ returns a different color based on filter attribute"""
        return colors[color_dict[legend_group(filter)]]

    
    def opacity_group(filter):
        """ returns a different color based on filter attribute"""
        return opacity_dict[legend_group(filter)]


    # Create figure with inset
    fig = go.Figure()

    # Create trajectories for legend groups
    # Fit trajectories
    fig.add_trace(
        go.Scatter(x=[fit_traj.data.AF.max()],
                    y=[fit_traj.gradient],
                    name='Fit mutation',
                    mode='markers',
                    marker_color=color_group(fit_traj.filter),
                    marker_size=5,
                    legendgroup=legend_group(fit_traj.filter),
                    showlegend=True))

    # Neutral trajectories
    fig.add_trace(
        go.Scatter(x=[artifact_traj.data.AF.max()],
                    y=[artifact_traj.gradient],
                    name='<2%',
                    mode='markers',
                    marker_color=color_group(artifact_traj.filter),
                    marker_size=5,
                    showlegend=True))

    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is False:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=5,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)


    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=5,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)


    # Add vertical delimitation of threshold
    # fig.add_vline(x=VAF_threshold,
    #               line_width=3, line_dash="dash", line_color=colors[2])
    fig.add_vrect(x0=0.01, x1=0.02, line_width=0, fillcolor="red", opacity=0.2)

    fig.update_layout(
        title=None,
        template='simple_white',
        xaxis=dict(title='VAF', type='log', tickvals=[0.01, 0.02, 0.05, 0.1, 0.2]),
        yaxis=dict(title='Gradient',
                   range=[-0.02, 0.02]
                  )
    )

    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)


    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0)
    )

    return fig


def gradient_LiFT_plot(cohort, syn):
    """ Plot the selection of filtered variants.
    threshold: Default = 0.02.
    x_axis: VAF.
    y_axis: gradient.
    color: filter attribute.
    """


    # Find fit and netural trajectories to create group
    fit_traj = find_fit(cohort)
    artifact_traj = find_neutral(cohort)

    def legend_group(filter):
        """Legend group assigns a legend group based on filter attribute"""
        if filter is True:
            return 'fit mutation'
        else:
            return 'artifact'

    def color_group(filter):
        """ returns a different color based on filter attribute"""
        return colors[color_dict[legend_group(filter)]]

    
    def opacity_group(filter):
        """ returns a different color based on filter attribute"""
        return opacity_dict[legend_group(filter)]

    # Create figure with inset
    fig = go.Figure()

    # Add vertical delimitation of threshold
    fig.add_vrect(x0=0.01, x1=0.02, line_width=0, fillcolor="red", name='2% filter' ,opacity=0.05)

    # Create trajectories for legend groups
    # Fit trajectories
    fig.add_trace(
        go.Scatter(x=[fit_traj.data.AF.max()],
                    y=[fit_traj.gradient],
                    name='Fit mutation',
                    mode='markers',
                    marker_color=color_group(fit_traj.filter),
                    marker_size=5,
                    legendgroup=legend_group(fit_traj.filter),
                    showlegend=True))


    # Neutral trajectories
    fig.add_trace(
        go.Scatter(x=[artifact_traj.data.AF.max()],
                    y=[artifact_traj.gradient],
                    name='Filtered',
                    mode='markers',
                    marker_color=color_group(artifact_traj.filter),
                    marker_size=5,
                    showlegend=True))

    # Plot synonymous variants
    syn_traj = syn[0].trajectories[0]
    fig.add_scatter(
        x=[syn_traj.data.AF.max()],
        y=[syn_traj.gradient],
        mode='markers',
        name='Synonymous',
        marker_color=colors[2],
        marker_size=5,
        marker_opacity=0.3,
        legendgroup=legend_group('synonymous'),
        showlegend=True)

    for part in syn:
        for traj in part.trajectories:
            fig.add_scatter(
                x=[traj.data.AF.max()],
                y=[traj.gradient],
                mode='markers',
                marker_color=colors[2],
                marker_size=5,
                marker_opacity=0.3,
                legendgroup=legend_group('synonymous'),
                showlegend=False)

    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is False:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=5,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)


    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                fig.add_scatter(x=[traj.data.AF.max()],
                                y=[traj.gradient],
                                mode='markers',
                                marker_color=color_group(traj.filter),
                                marker_size=5,
                                marker_opacity=opacity_group(traj.filter),
                                legendgroup=legend_group(traj.filter),
                                showlegend=False)

    fig.update_layout(
        title=None,
        template='simple_white',
        xaxis=dict(title='VAF', type='log', tickvals=[0.01, 0.02, 0.05, 0.1, 0.2]),
        yaxis=dict(title='Gradient',
                   range=[-0.02, 0.02]
                  )
    )

    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)


    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0)
    )
    fig.update_layout(margin=dict(pad=4))

    return fig

def gene_trajectories_legend(cohort):
    # Load variant color dictionary
    with open("../Resources/var_dict.json") as json_file:
        var_dict = json.load(json_file)

    default_color="Grey"

    genes = ['DNMT3A', 'JAK2', 'TET2', 'ASXL1']

    cohort_variants = set([traj.variant_class
                           for part in cohort
                             for traj in part.trajectories
                               if traj.gene in genes])

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(x=[0], y=[0], mode='markers', marker_color='grey',
                name='Data points',
                showlegend=True))

    fig.add_trace(
        go.Scatter(x=[0,1], y=[0,1], mode='lines', marker_color='grey',
                    name='Fitted trajectories',
                showlegend=True))

    for variant in cohort_variants:
        fig.add_trace(
            go.Scatter(x=[0], y=[0], mode='markers', name=variant, marker_color=var_dict[variant]))
    fig.update_layout(margin=dict(pad=4))
    return fig