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
    - p_key: string. Amino acid change.
    - emcee: bool, default=False. Traing an MCMC.
    """

    # mutation trajectory class for training mp models
    def __init__(self, x, y, gene, mutation, id, p_key, emcee=False):
        self.x = x
        self.y = y
        self.mutation = mutation
        self.gene = gene
        self.id = id
        self.p_key = p_key
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
        value = np.random.normal(-0.1, 0.01)
        p.add('fitness', value=value, vary=False)

    p.add('neutral_cells', value=500, min=1)
    p.add('origin', value=50, min=0, max=80)

    # step 1: Crete Minimizer object model
    self.model = Minimizer(exponential_residual, params=p,
                           fcn_args=(self.x, self.y))

    # step 2: minimize mini with Nelder-Mead algorithm (robust)
    self.nelder = self.model.minimize(method='nelder')
    self.fitness = self.nelder.params['fitness']
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


def extended_plot_param(model_list, mutation, col_param=0):
    """Produce a side-by-side plot of inferred trajectories
    and fitness distribution for all for all trajectories with a mutated gene.
    """

    if col_param == 'prediction':
        def color(traj):
            return traj.damage
    else:
        def color(traj):
            return 0

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
                   marker=dict(color='Grey'),
                   legendgroup="group",
                   name="Data points"),
        row=1, col=1)

    fig.add_trace(
        go.Scatter(x=x_label, y=y_label,
                   mode='lines', line=dict(color='Grey'),
                   legendgroup="group",
                   name="Fitted trajectories"),
        row=1, col=1)

    fitness = []  # create the list of fitness to use in boxplot
    # Create a trace for each traject
    for i, traj in enumerate(model_filtered):
        x_line = np.linspace(65, 95, 100)
        y_fit = prediction(traj.nelder.params, x_line)

        fig.add_trace(
            go.Scatter(x=traj.x, y=traj.y,
                       name=traj.mutation,
                       legendgroup="group2",
                       mode='markers',
                       line=dict(color=colors[color(traj)])),
            row=1, col=1)

        fig.add_trace(
            go.Scatter(x=x_line, y=y_fit,
                       name='fit', legendgroup='group2',
                       showlegend=False,
                       line=dict(color=colors[color(traj)])),
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
                             marker_color=colors[color(model_filtered[i])],
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
                        column_widths=[0.7, 0.3],
                        subplot_titles=(model_filtered[0].model_type,
                                        'Fitness distribution'))

    # Create legend groups
    traj = model_filtered[0]
    x_label = np.linspace(70, 71, 2)
    y_label = prediction(traj.nelder.params, x_label)

    fig.add_trace(
        go.Scatter(x=[traj.x[0]], y=[traj.y[0]],
                   mode='markers',
                   marker=dict(color='Grey'),
                   legendgroup="group",
                   name="Data points"),
        row=1, col=1)

    fig.add_trace(
        go.Scatter(x=x_label, y=y_label,
                   mode='lines', line=dict(color='Grey'),
                   legendgroup="group",
                   name="Fitted trajectories"),
        row=1, col=1)

    fitness = []  # create the list of fitness to use in boxplot
    # Create a trace for each traject
    for i, traj in enumerate(model_filtered):
        x_line = np.linspace(65, 95, 100)
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
                             width=0.1,
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


def gene_box(cohort, order='mean'):
    """Box plot with counts of filtered mutations by gene. Returns a figure."""

    cohort = [item for item in cohort if item.fitness > 0]
    # Create a dictionary with all filtered genes
    gene_list = []
    for traj in cohort:
        gene_list.append(traj.gene)
    gene_dict = {element: [] for element in set(gene_list)}

    # update the counts for each gene
    for traj in cohort:
        fitness = traj.nelder.params['fitness']*1
        gene_dict[traj.gene].append(fitness)
    # sort dictionary in descending order
    if order == 'mean':
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: np.mean(item[1]),
                                reverse=True))

    if order == 'max':
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: np.max(item[1]),
                                reverse=True))
    # Bar plot
    fig = go.Figure()
    for key in gene_dict:
        fig.add_trace(
            go.Box(y=gene_dict[key],
                   name=key, boxpoints='all', showlegend=False))
    fig.update_layout(title='Gene distribution of filtered mutations',
                      yaxis_title='Fitness',
                      template="simple_white")
    return fig


#
# """Logistic model fit functions"""
#
#
# def logistic_dataset(params, i, x):
#     """Calculate logistic lineshape from parameters for data set."""
#     amp = params['amp_%i' % (i+1)]
#     fit = params['fit_%i' % (i+1)]
#     dis = params['dis_%i' % (i+1)]
#     return logistic(x, amp, fit, dis)
#
#
# def logistic_model(self, l2, common_fit):
#     # l2 = weight of the regularization of amplitude parameter around 0.5
#     # initialize parameters using default values and boundary values
#     fit_params = Parameters()
#     for iy, y in enumerate(self.y):
#         fit_params.add('amp_%i' % (iy+1),
#                        value=0.4, min=0.2, max=0.5)
#         fit_params.add('fit_%i' % (iy+1), value=0)
#         fit_params.add('dis_%i' % (iy+1), value=0)
#
#     if common_fit is True:
#         # force all fit parameters equal
#         for iy in range(2, len(self.y)+1):
#             fit_params['fit_%i' % iy].expr = 'fit_1'
#
#     # fit the data using lmfit 'minimize' function
#     return minimize(logistic_objective, fit_params, args=(self.x, self.y, l2))
#
#
# def logistic_objective(params, x, y, a):
#     """Calculate total residual for fits of logistics to several data sets."""
#
#     total_res = []   # keeps track of the cumulative residual
#
#     # compute the residual for each trajectory in the dataset
#     for i, points in enumerate(zip(x, y)):
#         res = points[1] - logistic_dataset(params, i, points[0])
#         total_res.append(res)
#
#     # each entry of total_res is the list of residuals for each trajectory
#     # as minimize() needs a 1D array -> flatten
#     total_res = np.concatenate(total_res)
#     for i, j in enumerate(y):
#         total_res = np.append(total_res,
#                               np.asarray([a*(0.5-params['amp_%i' % (i+1)])]))
#
#     return total_res
#
#
# """ Logistic model emcee"""
#
#
# def logistic(x, amplitude, fitness, origin, N0=0.004):
#     """logistic lineshape."""
#     sig_center = origin + (1/fitness) * np.log((amplitude-N0)/N0)
#     return amplitude / (1. + np.exp(- fitness * (x - sig_center)))
#
#
# def residual(p, x, y, initial_pop=0.004):
#     """Calculate total residual for fits of logistics to several data sets."""
#     v = p.valuesdict()
#     res = y - logistic(x, v['amplitude'],
#                        v['fitness'], v['origin'], initial_pop)
#     res = np.append(res, v['regularization']*(0.5-v['amplitude']))
#
#     return res
#
#
# def competitive_fit(self,  regularization=0.05):
#     self.model_type = 'Commpetitive'
#
#     # Create model parameters
#     p = Parameters()
#     p.add('amplitude', value=0.4, min=0.2, max=0.5)
#     p.add('fitness', value=0.1, min=0)
#     p.add('origin', value=30, min=0, max=70)
#     p.add('regularization', value=regularization, vary=False)
#
#     # step 1: Crete Minimizer object model
#     self.model = Minimizer(residual, params=p,
#                            fcn_args=(self.x, self.y))
#     # step 2: minimize mini with Nelder-Mead algorithm (robust)
#     self.nelder = self.model.minimize(method='nelder')
#     # step 3: fine tune minimization with Levenberg-Marquardt algorithm ()
#     # This allows the computation of confidence intervals
#     self.ls = self.model.minimize(method='leastsq',
#                                   params=self.nelder.params.copy())
#
#     if self.emcee is True:
#         # add sigma parameter (~exp(variance) of the error distribution)
#         self.nelder.params.add('__lnsigma',
#                                value=np.log(0.1),
#                                min=np.log(0.001), max=np.log(2))
#
#         self.emcee = minimize(residual, method='emcee', seed=1,
#                               nan_policy='omit', burn=300, steps=1000,
#                               nwalkers=100, thin=1,
#                               params=self.nelder.params,
#                               args=(self.x, self.y),
#                               is_weighted=False, progress=False)
#
#     return self
#
#
