# Last modification: E. Latorre -- 25-10-2020

# Collection of plotting functions

import pandas as pd
from scipy import stats

import plotly.express as px
import plotly.graph_objects as go
colors = px.colors.qualitative.Plotly


""" Plotting functions for longitudinal trajectories"""


def plot_id(cohort, participant_id, germline=False):
    # Plot trajectories by participant_id
    for part in cohort:
        if part.id == participant_id:
            fig = part.profile(germline=germline)
    return fig


def stack_plot(participant, norm=True):
    # stacked plot of participant
    # norm=True plots % of VAF contribution rather than VAF
    fig = go.Figure()
    if norm is True:
        for traj in participant.trajectories:
            fig.add_trace(go.Scatter(x=traj.data.age,
                                     y=traj.data.AF,
                                     name=f'{traj.mutation}',
                                     hovertemplate='VAF: %{y}',
                                     mode='lines',
                                     stackgroup='one',  # define stack group
                                     groupnorm='percent'))

        fig.update_layout(showlegend=True,
                          yaxis=dict(type='linear',
                                     range=[1, 100],
                                     ticksuffix='%'))
    else:
        for traj in participant.trajectories:
            fig.add_trace(go.Scatter(x=traj.data.age,
                                     y=traj.data.AF,
                                     name=f'{traj.mutation}',
                                     hovertemplate='VAF: %{y}',
                                     mode='lines',
                                     stackgroup='one'))

    return fig


""" Plotting functions for mutation trajectories"""


def mutation(cohort, mutation):
    # plot all trajectories with a mutations
    fig = go.Figure()

    for part in cohort:
        for i, word in enumerate(part.mutation_list):
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


""" Plotting functions for mutation statistics"""


def top_bar(cohort, n_genes=10, all=False):
    gene_list = []
    for part in cohort:
        for traj in part.trajectories:
            gene_list.append(traj.mutation.split()[0])
    gene_dict = {element: 0 for element in set(gene_list)}
    for part in cohort:
        for traj in part.trajectories:
            gene_dict[traj.mutation.split()[0]] = gene_dict[traj.mutation.split()[0]] + 1
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


""" Plotting functions for gradients statistics"""


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


def local_gardients(cohort, mutations):
    # violin plots of local gradients by mutation
    data = pd.DataFrame(columns=['gradient', 'participant', 'mutation'])

    for part in cohort:
        for traj in part.trajectories:
            if traj.mutation.split()[0] in mutations:
                for i, row in traj.data.iterrows():
                    data = data.append({'gradient': row.regularized_gradient,
                                        'participant': part.id,
                                        'mutation': traj.mutation.split()[0]},
                                       ignore_index=True)
    # violin plot of data
    fig = px.box(data, y="gradient", x='mutation',
                 color='mutation', points='all',
                 hover_name="participant", hover_data=["mutation"])

    fig.update_layout(title='Local gradients by mutation')
    return fig


""" Comparing fitness with Protein damage prediction"""


def damage_class(model, fs_damage=1, termination_damage=1):
    """ Box plot of variant predicted protein damage class ~ Fitness.
    Returns:
    * Box plot
    * Modified model with damage_class attribute."""

    # Load dataset with prediction damage
    xls = pd.ExcelFile('Datasets/new_kristina_variants.xlsx')
    df = pd.DataFrame(columns=['p_key', 'damaging'])

    # Access each sheet of xls file
    for name in xls.sheet_names:
        df_temp = pd.read_excel(xls, name)

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
            traj.damage = df.loc[df.p_key == traj.p_key, 'damaging'].values[0]
        elif 'fs' in traj.p_key:
            traj.damage = fs_damage
        elif '*' in traj.p_key:
            traj.damage = termination_damage
        else:
            traj.damage = -2
    damaging_var = []
    possibly_var = []
    benign_var = []
    for traj in damage_model:
        if traj.fitness > 0:
            if traj.damage == 1:
                damaging_var.append(traj.nelder.params['fitness'])
            elif traj.damage == 0:
                possibly_var.append(traj.nelder.params['fitness'])
            elif traj.damage == -1:
                benign_var.append(traj.nelder.params['fitness'])
    fig = go.Figure()
    fig.add_trace(
        go.Box(y=damaging_var, name='Likely damaging', boxpoints='all'))
    fig.add_trace(
        go.Box(y=possibly_var, name='Possibly damaging', boxpoints='all'))
    fig.add_trace(
        go.Box(y=benign_var, name='Likely benign', boxpoints='all'))
    fig.update_yaxes(title='Fitness', type='log', dtick=1)

    return fig, damage_model


def damage_score(model, fs_score=-1, termination_score=-1,
                 unknown_score=-1):
    """ Scatter plot of variant predicted protein damage ~ Fitness
        Returns:
        * Scatter plot
        * Modified model with damage_class attribute."""

    # Load xls file
    xls = pd.ExcelFile('Datasets/new_kristina_variants.xlsx')
    df = pd.DataFrame(columns=['p_key', 'damaging'])
    # For each sheet, compute the mean of gnomad and clinvar
    for name in xls.sheet_names:
        df_temp = pd.read_excel(xls, name)
        df_temp.columns = [col.lower() for col in df_temp.columns]

        # Detect the presence of gnomad mean and clinvar mean columns
        predictors = {'gnomad mean', 'clinvar mean'}
        col = set(df_temp.columns)
        inter = predictors.intersection(col)

        # calculate the mean were possible
        df_temp['score'] = df_temp[inter].mean(axis=1)

        # extract p_keys and score for each variant
        p_keys = name + ' p.' + df_temp.iloc[:, 0]
        score = df_temp['score']
        for i, j in zip(p_keys, score):
            df = df.append({'p_key': i, 'score': j}, ignore_index=True)

    # Exclude trajectories without a p_key
    damage_model = [traj for traj in model if isinstance(traj.p_key, str)]

    # Assign damage to trajectories
    for traj in damage_model:
        if traj.p_key in set(df['p_key']):
            traj.damage = df.loc[df.p_key == traj.p_key, 'score'].values[0]
        elif 'fs' in traj.p_key:
            traj.damage = fs_score
        elif '*' in traj.p_key:
            traj.damage = termination_score
        else:
            traj.damage = unknown_score
        damaging_var = []
        possibly_var = []
        benign_var = []
        for traj in damage_model:
            if traj.fitness > 0:
                if traj.damage == 1:
                    damaging_var.append(traj.nelder.params['fitness'])
                elif traj.damage == 0:
                    possibly_var.append(traj.nelder.params['fitness'])
                elif traj.damage == -1:
                    benign_var.append(traj.nelder.params['fitness'])

    fitness = []
    damage = []
    fig = go.Figure()
    for traj in damage_model:
        if 0 < traj.damage < 1:
            fig.add_trace(
                go.Scatter(x=[traj.fitness], y=[traj.damage],
                           marker_color=colors[0],
                           showlegend=False))
            fitness.append(traj.fitness)
            damage.append(traj.damage)

    tau, p_value = stats.kendalltau(fitness, damage)

    fig.update_layout(
        xaxis_title='Fitness',
        yaxis_title='damage prediction',
        title=(f'Fitness ~ damage prediction: '
               f'tau: {round(tau,2)} p-value: {round(p_value,3)}'))

    return fig, damage_model
