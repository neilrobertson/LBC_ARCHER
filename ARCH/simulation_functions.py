import numpy as np
import random
import pandas as pd
from itertools import combinations
from scipy.integrate import quad
import plotly.express as px
import plotly.graph_objects as go
colors = px.colors.qualitative.Plotly


def bd(t0, lamb=0.125, mu=0.125, tmax=80, x0=1,
       timepoints=[70, 73, 76, 79], extintion=True):
    ''' Stochastic B-D model of a clone
        t0 = initiating time
        lamb = birth rate / year
        mu = death rate / year
        tmax = age of individuals for stopping time
        x0 = initial population of the clone
        HSC = amount of hematopoietic stem cells
    '''

    # set a seed for reproducibility
    # the seed shouldn't be applied for the test scenario of
    # all mutations starting at 0
    if t0 > 0:
        random.seed(t0)

    # initialize time and clone population
    t = t0
    x = x0

    # initialize dataframe to record the trajectory
    history = pd.DataFrame(columns=['count', 'time'])
    checkpoint = pd.DataFrame(columns=['count', 'time'])

    # append first time point to history dataframe
    history = history.append({'count': x, 'time': t}, ignore_index=True)

    # model loop
    while t <= tmax:
        # save current time in previous_t
        previous_t = t

        # update birth and death rates
        lamb_n = lamb * x
        mu_n = mu * x

        #  Timestep to next event
        rate = lamb_n + mu_n
        t += -np.log(random.random())/rate

        # append x to all remaining checkpoints in (previous_t, tmax)
        if t > tmax:
            for time in timepoints:
                if time > previous_t:
                    checkpoint = checkpoint.append({'count': x,
                                                    'time': float(time)},
                                                   ignore_index=True)
            history = history.append({'count': x,
                                      'time': float(tmax)},
                                     ignore_index=True)
            break

        # If the next event happens after a checkpoint then
        # save the population at eatch checkpoint
        for time in timepoints:
            if previous_t <= time < t:
                checkpoint = checkpoint.append({'count': x,
                                                'time': float(time)},
                                               ignore_index=True)

        # Decision of which event happened during time-step
        urv = random.random()  # uniform random variable
        if urv >= lamb_n / rate:
            x -= 1
        else:
            x += 1

        # Append the current population and time of clone in dataframe
        history = history.append({'count': x, 'time': t}, ignore_index=True)

        if x == 0:  # append 0 to all remaining checkpoints in (t, tmax)
            if previous_t > timepoints[0]:
                for time in timepoints:
                    if time > t:
                        checkpoint = checkpoint.append({'count': x,
                                                        'time': float(tmax)},
                                                       ignore_index=True)
            history = history.append({'count': x, 'time': float(tmax)},
                                     ignore_index=True)
            break     # clone death

    # A row at tmax could have been duplicated
    if checkpoint.empty is False:
        checkpoint = checkpoint.drop_duplicates()

    return history, checkpoint


def integrand(x, lamb, n):
    '''Expected number of families of n elements at time t'''

    a = 1 / (1 + lamb * x)
    b = (lamb * x) / (1 + lamb * x)

    return np.power(a, 2)*np.power(b, n-1)


def find_mutation_rate(lamb, N, threshold):
    """ Given a birth rate (lamb), number of HSC (N) and a threshold,
    return the immigration rate so that under the dynamics of a CBM,
    no families are expected to reach threshold (VAF)"""
    for mu in range(1, 100):
        n = 0
        E_n = 2
        while E_n > 1:
            n += 1
            E_n = mu*quad(integrand, 0, 80, args=(lamb, n))[0]

        if n / N > threshold:
            return mu
    return 0


def plot_result(result):
    fig = go.Figure()

    # Add a trajectory for labeling
    fig.add_trace(
        go.Scatter(x=[0], y=[0], mode='lines',
                   legendgroup="Mutation trajectory",
                   name="Mutation trajectory"))
    for traj in result:
        fig.add_trace(
            go.Scatter(x=traj[0]['time'], y=traj[0]['count'],
                       mode='lines', legendgroup="Mutation trajectory",
                       showlegend=False))
    fig.update_layout(
        title="Simulation of neutral mutations during a human lifespan",
        xaxis_title="Time (age)",
        yaxis_title="HSC counts"
    )
    return fig


def plot_VAF(result):
    # Remove unobservable events of extinction
    result[:] = [traj[traj['count'] >= 1] for traj in result]

    fig = go.Figure()
    # Add a trajectory for labeling
    fig.add_trace(
        go.Scatter(x=[70], y=[0], mode='lines',
                   marker_color='rgb(204,204,204)',
                   legendgroup="Neutral trajectory",
                   name="Neutral mutation"))

    for traj in result:
        fig.add_trace(
            go.Scatter(x=traj['time'], y=traj['VAF'],
                       marker_color='rgb(204,204,204)',
                       mode='lines',
                       legendgroup="Neutral trajectory", showlegend=False))
    fig.update_layout(
        title="Simulation of neutral fitness mutations",
        xaxis_title="Time (age)",
        yaxis_title="VAF"
    )
    return fig


def plot_vignette(checkpoint_range):
    # Find the distribution of gradients of all trajectories
    melt_syn = pd.DataFrame(columns=['VAF', 'regularized_gradient'])

    # for all mutations but the fit (last one) compute local reg_gradients
    for traj in checkpoint_range[:-1]:
        for combination in list(combinations(traj.index, 2)):
            data = traj.loc[list(combination)]
            gradient = np.diff(data.VAF) / np.sqrt(np.diff(data.time))
            melt_syn = (
                melt_syn.append({'VAF': traj.loc[combination[0]]['VAF'],
                                 'regularized_gradient': gradient[0]},
                                ignore_index=True))

    no_filter = go.Figure()
    no_filter.add_trace(
        go.Scatter(x=melt_syn['VAF'], y=melt_syn['regularized_gradient'],
                   mode='markers', name='Neutral growth trajectories'))
    no_filter.update_layout(title='Gradient distribution of neutral mutations',
                            xaxis_title="VAF", yaxis_title='Gradient')

    # Find the distribution of gradients
    melt_syn = pd.DataFrame(columns=['VAF', 'regularized_gradient'])

    # for all mutations but the fit (last one) compute local reg_gradients
    for traj in checkpoint_range[:-1]:
        if traj['VAF'].max() > 0.01:
            for combination in list(combinations(traj.index, 2)):
                data = traj.loc[list(combination)]
                gradient = np.diff(data.VAF) / np.sqrt(np.diff(data.time))
                melt_syn = (
                    melt_syn.append({'VAF': traj.loc[combination[0]]['VAF'],
                                     'regularized_gradient': gradient[0]},
                                    ignore_index=True))
    filter = go.Figure()
    filter.add_trace(
        go.Scatter(x=melt_syn['VAF'], y=melt_syn['regularized_gradient'],
                   mode='markers', name='Neutral growth trajectories'))
    filter.update_layout(title='Gradient distribution of neutral mutations',
                         xaxis_title="VAF", yaxis_title='Gradient')

    # Create VAF_bins
    melt_syn['bin'] = pd.cut(melt_syn['VAF'], 25)

    filter_bin = go.Figure()
    filter_bin.add_trace(
        go.Scatter(x=[melt_syn.iloc[0]['VAF']],
                   y=[melt_syn.iloc[0]['regularized_gradient']],
                   mode='markers', showlegend=False,
                   legendgroup='group'))
    for bin in melt_syn.bin.unique():
        data = melt_syn[melt_syn['bin'] == bin]
        filter_bin.add_trace(
            go.Scatter(x=data['VAF'], y=data['regularized_gradient'],
                       mode='markers', showlegend=False,
                       legendgroup='group'))
        filter_bin.update_layout(
            title='Gradient distribution of neutral mutations',
            xaxis_title="VAF", yaxis_title='Gradient')

    filtered_melt_syn = melt_syn[melt_syn['VAF'] > 0.01]
    return filtered_melt_syn, no_filter, filter, filter_bin
