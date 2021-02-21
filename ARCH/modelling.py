# Last modification: E. Latorre -- 19-02-2021

# Module containing a model class definition and functions to fit exponential
# and logistic models of growth.

# import math packages
import numpy as np
import random

# import lmfit packages
from lmfit import Parameters, Minimizer, minimize

# import plotting packages
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
colors = px.colors.qualitative.Plotly  # create a list of plotly colors


class model:
    """Create a model class from a VAF trajectory.

    Each returned model has the following attributes:
    - x: float. Time-coordinate of measurement.
    - y: float. VAF-coordinate of measurement.
    - gene: string. Gene where mutation is present.
    - mutation: string. Mutation variant.
    - id: string. LBC id of participant harbouring this mutation.
    - damage: float. Clinvar damage prediction of variant.
    - damaging_class: bool, default=None. Predicted damaging_class by Joe.
    - known: bool, default=None. Previously reported variant.
    - emcee: bool, default=False. Traing an MCMC.
    """

    # mutation trajectory class for training mp models
    def __init__(self, x, y, gene, mutation, id, damage=None,
                 damaging_class=None, known=None, emcee=False):
        self.x = x
        self.y = y
        self.mutation = mutation
        self.gene = gene
        self.id = id
        self.damage = damage
        self.damaging_class = damaging_class
        self.known = known
        self.emcee = emcee


""" Exponential model"""


def exponential(t, neutral_cells, fitness, origin):
    """exponential model lineshape"""
    clone_size = np.exp(fitness*(t-origin))
    total_cells = neutral_cells + clone_size
    return clone_size / (2*total_cells)


def exponential_fit(self):
    """Fit an exponentially growing population of cells to a VAF
    trajectory.

    This model uses 3 parameters:
    - fitness   (prolifertive advantage conferred by a mutation)
    - origin    (time of mutation acquisition)
    - neutral_cells   (constant population of neutral cells)

    Note that origin is fully correlated with neutral_cells
    """

    self.model_type = '3-parameter exponential model'

    # Initialize model parameters
    p = Parameters()

    # Create a fit attribute to detect growing trajectories
    # if fit = True append fitness parameter for fitting else fix fitness=-0.1
    if self.y[-1]-self.y[0] > 0:
        self.fit = True
        p.add('fitness', value=0.5, min=0)
    else:
        self.fit = False
        p.add('fitness', value=-0.1, vary=False)

    p.add('neutral_cells', value=50, min=1)
    p.add('origin', value=50, min=0, max=80)

    # step 1: Crete Minimizer object model
    self.model = Minimizer(exponential_residual, params=p,
                           fcn_args=(self.x, self.y))

    # step 2: minimize mini with Nelder-Mead algorithm (robust)
    self.nelder = self.model.minimize(method='nelder')

    # step 3: fine tune minimization with Levenberg-Marquardt algorithm ()
    # This allows the computation of confidence intervals
    self.ls = self.model.minimize(method='leastsq',
                                  params=self.nelder.params.copy())

    if self.emcee is True:
        # add sigma parameter (~exp(variance) of the error distribution)
        self.nelder.params.add('__lnsigma',
                               value=np.log(0.1),
                               min=np.log(0.001), max=np.log(2))

        self.emcee = minimize(exponential_residual, method='emcee', seed=1,
                              nan_policy='omit', burn=300, steps=1000,
                              nwalkers=100, thin=1,
                              params=self.nelder.params,
                              args=(self.x, self.y),
                              is_weighted=False, progress=False)

    return self


def exponential_residual(p, x, y):
    """Calculate total residual for fits of logistics to several data sets."""
    # translate parameter values to dict and compute predicted xalue for x
    v = p.valuesdict()
    prediction = exponential(x, v['neutral_cells'], v['fitness'], v['origin'])
    return prediction - y


""" Exponential fit with only 2 parameters"""


def exponential_2(t, aux, fitness):
    """exponential model lineshape"""
    clone_size = np.exp(t*fitness)
    total_cells = aux + clone_size
    return clone_size / (2*total_cells)


def exponential_fit_2(self):
    """Fit an exponentially growing population of cells to a VAF
    trajectory.

    This model uses 3 parameters:
    - fitness     (prolifertive advantage conferred by a mutation)
    - auxiliary   (neutral_cells*np.exp(fitness*origin))
    """

    self.model_type = '2-parameter exponential model'

    # Create model parameters
    p = Parameters()
    p.add('auxiliary', value=100, min=10)

    # Create a fit attribute to detect growing trajectories
    # if fit = True append fitness parameter for fitting else fix fitness=-0.1
    if self.y[-1]-self.y[0] > 0:
        self.fit = True
        p.add('fitness', value=0.5, min=0)
    else:
        self.fit = False
        p.add('fitness', value=-0.1, vary=False)

    # step 1: Crete Minimizer object model
    self.model = Minimizer(exponential_residual_2, params=p,
                           fcn_args=(self.x, self.y))
    # step 2: minimize mini with Nelder-Mead algorithm (robust)
    self.nelder = self.model.minimize(method='nelder')

    if self.emcee is True:
        # step 3: fine tune minimization with Levenberg-Marquardt algorithm ()
        # This allows the computation of confidence intervals
        self.ls = self.model.minimize(method='leastsq',
                                      params=self.nelder.params.copy())

        # add sigma parameter (~exp(variance) of the error distribution)
        self.nelder.params.add('__lnsigma',
                               value=np.log(0.1),
                               min=np.log(0.001), max=np.log(2))
        self.emcee = minimize(exponential_residual_2, method='emcee', seed=1,
                              nan_policy='omit', burn=300, steps=1000,
                              nwalkers=100, thin=1,
                              params=self.nelder.params,
                              args=(self.x, self.y),
                              is_weighted=False, progress=False)

    return self


def exponential_residual_2(p, x, y):
    """Calculate total residual for fits of logistics to several data sets."""

    # translate parameter values to dict and compute predicted xalue for x
    v = p.valuesdict()
    VAF_value = exponential_2(np.array(x), v['auxiliary'], v['fitness'])
    return VAF_value - np.array(y)


""" Auxiliary plotting functions"""


def participant_filter(part):
    """ Plot a participant's profile using different colors for
    fit and non-fit mutations"""

    fig = go.Figure()

    counter = 1  # Counter used for changing colors in non-fit trajectories

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
                   marker=dict(color='rgba(255,127,0, 0.2)'),
                   legendgroup='non-fit',
                   showlegend=True))

    for traj in part.trajectories:
        # Plot trajectories  using different colors according to traj.filter
        if traj.filter is False:
            fig.add_trace(
                go.Scatter(x=traj.data.age, y=traj.data.AF,
                           name='non-fit',
                           mode='lines',
                           marker=dict(color='rgba(255, 127, 0, 0.2)'),
                           legendgroup='non-fit',
                           showlegend=False))

        if traj.filter is True:
            counter += 1
            fig.add_trace(
                go.Scatter(x=traj.data.age, y=traj.data.AF,
                           name=traj.mutation,
                           mode='lines',
                           marker=dict(color=colors[counter % 10]),
                           legendgroup='fit',
                           showlegend=True))

    # Update layout
    fig.update_layout(title=f'Participant {part.id}',
                      xaxis_title='Age',
                      yaxis_title='VAF')

    return fig


def mutation_plot(trajectories, mutation):
    # plot fitted trajectories containing a mutation

    fig = go.Figure()
    x_line = np.linspace(0, 90, 1000)
    i = 0  # color counter
    for traj in trajectories:
        if traj.mutation == mutation or \
           traj.mutation.split()[0] == mutation:
            x_line = np.linspace(traj.nelder.params['origin'], 95, 1000)
            y_fit = exponential(x_line, traj.nelder.params['neutral_cells'],
                                traj.nelder.params['fitness'],
                                traj.nelder.params['origin'])

            fig.add_trace(go.Scatter(x=traj.x, y=traj.y,
                                     name=traj.mutation,
                                     mode='markers',
                                     line=dict(color=colors[i % 10-1])))

            fig.add_trace(go.Scatter(x=x_line, y=y_fit, name='fit',
                                     line=dict(color=colors[i % 10-1]),
                                     showlegend=False))
            i += 1  # update color counter
    fig.update_layout(title=f' {mutation} mutation fitted trajectories')
    fig.update_xaxes(title='Age in years')
    fig.update_yaxes(title='VAF')
    return fig


def extended_plot(model_list, mutation):
    """Produce a side-by-side plot of inferred trajectories
    and fitness distribution for all for all trajectories with a mutated gene.
    """

    model_filtered = [model for model in model_list
                      if model.gene == mutation
                      and model.nelder.params['fitness'] > 0
                      and model.nelder.aic < -10]

    if model_filtered[0].model_type == '3-parameter exponential model':
        def prediction(params, x):
            return exponential_residual(params, x, 0)

    elif model_filtered[0].model_type == '2-parameter exponential model':
        def prediction(params, x):
            return exponential_residual_2(params, x, 0)

    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=(model_filtered[0].model_type,
                                        'Fitness distribution'))

    # Create legend groups
    traj = model_filtered[0]
    x_label = np.linspace(70, 71, 2)
    y_label = prediction(traj.nelder.params, x_label)

    fig.add_trace(
        go.Scatter(x=[traj.x[0]], y=[traj.y[0]],
                   mode='markers',
                   marker=dict(color='rgba(0,0,0,0)'),
                   legendgroup="group",
                   name="Data points"),
        row=1, col=1)

    fig.add_trace(
        go.Scatter(x=x_label, y=y_label,
                   mode='lines', line=dict(color='rgba(0,0,0,0)'),
                   legendgroup="group",
                   name="Fitted trajectories"),
        row=1, col=1)

    fitness = []  # create the list of fitness to use in boxplot
    # Create a trace for each traject
    for i, traj in enumerate(model_filtered):
        x_line = np.linspace(65, 95)
        y_fit = prediction(traj.nelder.params, x_line)

        fig.add_trace(
            go.Scatter(x=traj.x, y=traj.y,
                       name=traj.mutation,
                       legendgroup="group2",
                       mode='markers',
                       line=dict(color=colors[i % 10-1])),
            row=1, col=1)

        fig.add_trace(
            go.Scatter(x=x_line, y=y_fit,
                       name='fit', legendgroup='group2',
                       showlegend=False,
                       line=dict(color=colors[i % 10-1])),
            row=1, col=1)

        fitness.append(traj.nelder.params['fitness'])  # append fitness to list

    # box plot with fitness estimates
    # Frist create the box
    fig.add_trace(go.Box(y=fitness, boxpoints=False,
                         name=mutation, marker_color=colors[0],
                         showlegend=False),
                  row=1, col=2)

    # Seccond create the dots adding jitter
    for i, item in enumerate(fitness):
        fig.add_trace(go.Box(y=[item], boxpoints='all', name=mutation,
                             jitter=random.uniform(0, 0.5),
                             marker_color=colors[i % 10-1],
                             line=dict(color='rgba(0,0,0,0)'),
                             fillcolor='rgba(0,0,0,0)',
                             showlegend=False),
                      row=1, col=2)

    fig.update_layout(title=f'{mutation} trajectories')
    fig.update_xaxes(title_text='Age (in years)', row=1, col=1)
    fig.update_yaxes(title_text='VAF', row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=1, col=2)
    fig.update_yaxes(title_text='Fitness', row=1, col=2)

    return fig
