# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================
"""Module to associate LiFT variants with structural damage predictions.
"""
# =============================================================================
#region: Import packages
# =============================================================================
# Import local packages
# =============================================================================
# Append root directory to system's path

import sys
sys.path.append('../ARCH_package')

import basic

# =============================================================================
# Import python packages
# =============================================================================
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio

default_color="Grey"
colors = px.colors.qualitative.D3
pio.templates.default = "simple_white"

#endregion
# =============================================================================
#region: Load damage prediction data
# =============================================================================

damage = pd.read_csv("../Datasets/LiFT_variants_damage.csv")

# Create likely damaging prediction
damage = damage.where(pd.notnull(damage), None)

# We base our categorical prediction based 
# on the mean of ClinVar mean and gnomAD mean
damage['mean'] = damage[['clinVar mean', 'gnomAD mean']].mean(axis=1)

# create categorical classes based on thresholds 
# [0.33, 0.66]
damage_mean = np.array(damage['mean'])
damage_mean = np.where(damage_mean < 0.33, -1, damage_mean)
damage_mean = np.where((damage_mean > 0) & (damage_mean < 0.66), 0, damage_mean)
damage_mean = np.where(damage_mean > 0.66, 1, damage_mean)

# Create likely damaging column
damage['likely_damaging'] = damage_mean

damage['likely_damaging'] = damage['likely_damaging'].replace(-1, 'likely benign')
damage['likely_damaging'] = damage['likely_damaging'].replace(0, 'possibly damaging')
damage['likely_damaging'] = damage['likely_damaging'].replace(1, 'likely damaging')

damage['p_key'] = damage['Gene'] + ' p.' + damage['Variant']

def make_damage_plot(cohort):
    """Create a box plot comparing fitness and predicted damage
    in a fitted cohort"""

    # Update damage class attribute for each trajectory in the cohort
    for part in cohort:
        for traj in part.trajectories:
            traj.damage_class = None
            if traj.p_key in list(damage.p_key):
                traj.damage_class = damage[damage.p_key == traj.p_key].likely_damaging.values[0]

            elif ('fs' in str(traj.p_key)) or ('*' in str(traj.p_key)):
                traj.damage_class = 'fs - ter'

    # filter trajectories in model without damage_class
    damage_model = [traj for part in cohort for traj in part.trajectories
                    if traj.damage_class is not None]

    # Create dataframe for plotting
    df = pd.DataFrame(columns=['damage_class', 'fitness'])
    for traj in damage_model:
        df = df.append({'damage_class': traj.damage_class,
                        '2%': max(traj.data.AF)>0.02,
                        'VAF':max(traj.data.AF),
                        'fitness': traj.fitness*1}, ignore_index=True)

    df = df[df.fitness > 0.02]
    df['2%'] = np.array(df['2%'], dtype=bool)

    # order df by damage class in a specific order
    df.damage_class = pd.Categorical(df.damage_class,
                                        categories=['likely benign',
                                                    'possibly damaging',
                                                    'likely damaging',
                                                    'fs - ter'],
                                        ordered=True)

    df.sort_values('damage_class', inplace=True)

    # Start subplot figure with shared xaxis
    fig = make_subplots(rows=2, cols=2,
                        row_heights=[0.3, 0.5],
                        column_widths=[0.5,0.1],
                        shared_xaxes=True,
                        shared_yaxes=True,
                        horizontal_spacing = 0.02,
                        vertical_spacing=0.01)

    name_dict ={False:'<2%', True:'>2%'}

    for i in [False, True]:
        df_2 = df[df['2%']==i]
        fig.add_trace(
            go.Histogram(x=df_2['damage_class'],
                            y=df_2['fitness'],
                            marker_color=colors[i],
                            name=name_dict[i],
                            showlegend=True),
                            row=1, col=1)

        fig.add_trace(
            go.Box(x=df_2['damage_class'],
                y=df_2['fitness'],
                boxpoints='all',
                marker_color=colors[i],
                showlegend=False),
            row=2, col=1)

        fig.add_trace(
            go.Histogram(
                y=df.loc[df['2%'] == i]['fitness'],
                marker_color=colors[i],
                showlegend=False),
                col=2,row=2)

    fig.update_layout(
        boxmode='group',
        legend_title='Maximum VAF'

    )

    fig.update_yaxes(title='counts', row=1, col=1)
    fig.update_yaxes(title='fitness',
                        linewidth=2,
                        row=2, col=1
                        )
    fig.update_xaxes(title='counts',
                     linewidth=2,
                     row=2, col=2
                        )

    fig.update_layout(boxgroupgap=0.5)
    fig.update_yaxes(linewidth=2)
    fig.update_xaxes(linewidth=2, row=2, col=1)
  
    fig.update_layout(legend=dict(
    yanchor="top",
    y=0.99,
    xanchor="right",
    x=0.99
    ))

    return fig




def extract_kde(data, plot_range=[0.02, 0.6], resolution=1000):
    """Create kde profile from data.
    Parameters:
    -----------
    data: List. List of data for which we want to extract its kde profile.
    quantiles: list. List of low and top range for kde evaluation.

    Returns:
    kde_profil: Array. 2-d array containing kde profile.
    -----------

    """
    kernel_dist = gaussian_kde(data)

    # Sample uniformly from the range of binomial probabilities
    kernel_sample = np.linspace(plot_range[0],
                                plot_range[1],
                                resolution)

    # Compute associated prior probabilities using kde
    kernel_values = kernel_dist.evaluate(kernel_sample)

    # normalise kernel to unit integral
    kernel_values = kernel_values/max(kernel_values)

    return np.array([kernel_sample, kernel_values])


def make_damage_plot_dist(cohort):
    """Create a box plot comparing fitness and predicted damage
    in a fitted cohort"""

    # Update damage class attribute for each trajectory in the cohort
    for part in cohort:
        for traj in part.trajectories:
            traj.damage_class = None
            if traj.p_key in list(damage.p_key):
                traj.damage_class = damage[damage.p_key == traj.p_key].likely_damaging.values[0]

            elif ('fs' in str(traj.p_key)) or ('*' in str(traj.p_key)):
                traj.damage_class = 'fs - ter'

    # filter trajectories in model without damage_class
    damage_model = [traj for part in cohort for traj in part.trajectories
                    if traj.damage_class is not None]

    # Create dataframe for plotting
    df = pd.DataFrame(columns=['damage_class', 'fitness'])
    for traj in damage_model:
        df = df.append({'damage_class': traj.damage_class,
                        '2%': max(traj.data.AF)>0.02,
                        'VAF':max(traj.data.AF),
                        'fitness': traj.fitness*1}, ignore_index=True)

    df = df[df.fitness > 0.02]
    df['2%'] = np.array(df['2%'], dtype=bool)

    # order df by damage class in a specific order
    df.damage_class = pd.Categorical(df.damage_class,
                                        categories=['likely benign',
                                                    'possibly damaging',
                                                    'likely damaging',
                                                    'fs - ter'],
                                        ordered=True)

    df.sort_values('damage_class', inplace=True)

    # Start subplot figure with shared xaxis
    fig = make_subplots(rows=2, cols=2,
                        row_heights=[0.3, 0.5],
                        column_widths=[0.5,0.1],
                        shared_xaxes=True,
                        shared_yaxes=True,
                        horizontal_spacing = 0.02,
                        vertical_spacing=0.01)

    name_dict ={False:'<2%', True:'>2%'}

    for i in [False, True]:
        df_2 = df[df['2%']==i]
        fig.add_trace(
            go.Histogram(x=df_2['damage_class'],
                            y=df_2['fitness'],
                            marker_color=colors[i],
                            name=name_dict[i],
                            showlegend=True),
                            row=1, col=1)

        fig.add_trace(
            go.Box(x=df_2['damage_class'],
                   y=df_2['fitness'],
                   boxpoints='all',
                   marker_color=colors[i],
                   showlegend=False),
            row=2, col=1)
        
        fitness_dist = extract_kde(df.loc[df['2%'] == i]['fitness'])
        fig.add_trace(
            go.Scatter(x=fitness_dist[1],
                        y=fitness_dist[0],
                        marker_color=colors[i],
                        showlegend=False),
                        row=2, col=2)

    fig.update_layout(
        boxmode='group',
        legend_title='Maximum VAF'

    )

    fig.update_yaxes(title='counts', row=1, col=1)
    fig.update_yaxes(title='fitness',
                        linewidth=2,
                        type='log',
                        tickvals=[0.02,0.05 ,0.1, 0.2, 0.4],
                        range=np.log10([0.02, 0.7]),
                        row=2, col=1
                        )
    fig.update_yaxes(
                    linewidth=2,
                    type='log',
                    tickvals=[0.02,0.05 ,0.1, 0.2, 0.4],
                    range=np.log10([0.02, 0.7]),
                    row=2, col=2
                    )


    fig.update_xaxes(title='distribution',
                     linewidth=2,
                     tickvals=[0,0.5,1],
                     row=2, col=2)

    fig.update_layout(boxgroupgap=0.5)
    fig.update_yaxes(linewidth=2)
    fig.update_xaxes(linewidth=2, row=2, col=1)
  
    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="right",
        x=0.99
        ))
    fig.update_layout(margin=dict(pad=4))


    return fig


def dnmt3a_damage(cohort):
    # Update damage class attribute for each trajectory in the cohort
    for part in cohort:
        for traj in part.trajectories:
            traj.damage_class = None
            if traj.p_key in list(damage.p_key):
                traj.damage_class = damage[damage.p_key == traj.p_key]['mean'].values[0]

    # filter trajectories in model without damage_class
    damage_model = [traj for part in cohort for traj in part.trajectories
                    if (traj.damage_class is not None and traj.gene =='DNMT3A')]

    fig = go.Figure()
    damage_below_2 = [traj.damage_class for traj in damage_model if traj.data.AF.max()<0.02]
    fig.add_trace(
        go.Box(y=damage_below_2, name='<2%', boxpoints='all')
    )

    damage_above_2 = [traj.damage_class for traj in damage_model if traj.data.AF.max()>0.02]
    fig.add_trace(
        go.Box(y=damage_above_2, name='>2%', boxpoints='all')
    )
    fig.add_hline(y=0.33, line_width=3, line_dash="dash", line_color="grey")
    fig.add_hline(y=0.66, line_width=3, line_dash="dash", line_color="grey")

    fig.update_yaxes(title='Damage prediction', range=[0.1,1])
    
    return fig


def make_damage_gene_plot_dist(cohort, gene='DNMT3A'):
    """Create a box plot comparing fitness and predicted damage
    in a fitted cohort"""

    # Update damage class attribute for each trajectory in the cohort
    for part in cohort:
        for traj in part.trajectories:
            traj.damage_class = None
            if traj.p_key in list(damage.p_key):
                traj.damage_class = damage[damage.p_key == traj.p_key]['likely_damaging'].values[0]


    # filter trajectories in model without damage_class
    damage_model = [traj for part in cohort for traj in part.trajectories
                    if (traj.damage_class is not None and traj.gene == gene)]

    # Create dataframe for plotting
    df = pd.DataFrame(columns=['damage_class', 'fitness'])
    for traj in damage_model:
        df = df.append({'damage_class': traj.damage_class,
                        '2%': max(traj.data.AF)>0.02,
                        'VAF':max(traj.data.AF),
                        'fitness': traj.fitness*1}, ignore_index=True)

    df = df[df.fitness > 0.02]
    df['2%'] = np.array(df['2%'], dtype=bool)

    # order df by damage class in a specific order
    df.damage_class = pd.Categorical(df.damage_class,
                                        categories=['likely benign',
                                                    'possibly damaging',
                                                    'likely damaging',
                                                    'fs - ter'],
                                        ordered=True)

    df.sort_values('damage_class', inplace=True)

    # Start subplot figure with shared xaxis
    fig = make_subplots(rows=2, cols=2,
                        row_heights=[0.3, 0.5],
                        column_widths=[0.5,0.1],
                        shared_xaxes=True,
                        shared_yaxes=True,
                        horizontal_spacing = 0.02,
                        vertical_spacing=0.01)

    name_dict ={False:'<2%', True:'>2%'}

    for i in [False, True]:
        df_2 = df[df['2%']==i]
        fig.add_trace(
            go.Histogram(x=df_2['damage_class'],
                            y=df_2['fitness'],
                            marker_color=colors[i],
                            name=name_dict[i],
                            showlegend=True),
                            row=1, col=1)

        fig.add_trace(
            go.Box(x=df_2['damage_class'],
                   y=df_2['fitness'],
                   boxpoints='all',
                   marker_color=colors[i],
                   showlegend=False),
            row=2, col=1)
        
        fitness_dist = extract_kde(df.loc[df['2%'] == i]['fitness'])
        fig.add_trace(
            go.Scatter(x=fitness_dist[1],
                        y=fitness_dist[0],
                        marker_color=colors[i],
                        showlegend=False),
                        row=2, col=2)

    fig.update_layout(
        boxmode='group',
        legend_title='Maximum VAF'

    )

    fig.update_yaxes(title='counts', row=1, col=1)
    fig.update_yaxes(title='fitness',
                        linewidth=2,
                        type='log',
                        tickvals=[0.02,0.05 ,0.1, 0.2, 0.4],
                        range=np.log10([0.02, 0.7]),
                        row=2, col=1
                        )
    fig.update_yaxes(
                    linewidth=2,
                    type='log',
                    tickvals=[0.02,0.05 ,0.1, 0.2, 0.4],
                    range=np.log10([0.02, 0.7]),
                    row=2, col=2
                    )


    fig.update_xaxes(title='distribution',
                     linewidth=2,
                     tickvals=[0,0.5,1],
                     row=2, col=2)

    fig.update_layout(boxgroupgap=0.5)
    fig.update_yaxes(linewidth=2)
    fig.update_xaxes(linewidth=2, row=2, col=1)
  
    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="right",
        x=0.99
        ))
    fig.update_layout(margin=dict(pad=4))


    return fig