# Last modification: E. Latorre -- 25-10-2020

# Collection of plotting functions

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

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
