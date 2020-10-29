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
            fig.add_trace(go.Scatter(x=traj.data.wave,
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
            fig.add_trace(go.Scatter(x=traj.data.wave,
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
                fig.add_trace(go.Scatter(x=traj.data.wave,
                                         y=traj.data.AF,
                                         mode='lines+markers',
                                         name=traj.mutation,
                                         hovertemplate=f"{part.id}"
                                         ))
    # Edit the layout
    fig.update_layout(title=f'Trajectories containing mutation {mutation}',
                      xaxis_title='Time (years since first wave)',
                      yaxis_title='VAF')

    return fig


def threshold(cohort, threshold, mutation='All'):
    # plot trajectories with gradient > threshold
    fig = go.Figure()

    if mutation == 'All':
        for part in cohort:
            for traj in part.trajectories:
                if traj.gradient > threshold:
                    fig.add_trace(go.Scatter(x=traj.data.wave,
                                             y=traj.data.AF,
                                             mode='lines+markers',
                                             name=part.id,
                                             hovertemplate=f"{part.id}"
                                             ))
    else:
        for part in cohort:
            for traj in part.trajectories:
                if (mutation in traj.mutation) and (traj.gradient > threshold):
                    fig.add_trace(go.Scatter(x=traj.data.wave, y=traj.data.AF,
                                             mode='lines+markers',
                                             name=part.id,
                                             hovertemplate=f"{part.id}"
                                             ))

    # Edit the layout
    fig.update_layout(title=f'Trajectories containing mutation {mutation}',
                      xaxis_title='Time (years since first wave)',
                      yaxis_title='VAF')

    return fig


""" Plotting functions for mutation statistics"""


def top_bar(cohort, n_mutations=5, l_vaf=0.1, u_vaf=0.4):
    m_count = pd.DataFrame(columns=['mutation', 'VAF_average'])
    for part in cohort:
        for traj in part.trajectories:
            # For each trajectory append mutated gene and average VAF
            m_count = m_count.append({'mutation': traj.mutation.split()[0],
                                      'AF_average': traj.data.AF.mean()},
                                     ignore_index=True)

    # Filter trajectories whose average is below VAF=0.1 or above 0.4
    m_count = m_count[m_count.AF_average.between(l_vaf, u_vaf)]

    top_count = pd.DataFrame(columns=['mutation', 'count', 'VAF_average'])

    for i in range(0, n_mutations):
        mutation = m_count['mutation'].value_counts().index[i]
        count = m_count['mutation'].value_counts()[i]
        AF = m_count.groupby('mutation').mean().loc[mutation][0]
        top_count = top_count.append({'mutation': mutation,
                                      'count': count,
                                      'VAF': AF},
                                     ignore_index=True)
    # bar plot of mutations
    fig = px.bar(top_count, x='mutation', y='count',
                 hover_data=['count', 'VAF'], color='VAF')
    fig.update_layout(title=f'Most commonly mutated genes',
                      xaxis_title='Gene',
                      yaxis_title='Count')

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
                    data = data.append({'gradient': row.gradient,
                                        'participant': part.id,
                                        'mutation': traj.mutation.split()[0]},
                                       ignore_index=True)
    # violin plot of data
    fig = px.box(data, y="gradient", x='mutation',
                 color='mutation', points='all',
                 hover_name="participant", hover_data=["mutation"])

    fig.update_layout(title='Local gradients by mutation')
    return fig
