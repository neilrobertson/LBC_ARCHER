# Last modification: E. Latorre -- 19-02-2021

# Module containing a model class definition and functions to fit exponential
# and logistic models of growth.

# import math packages
import numpy as np
import pandas as pd
import random
from scipy import stats
from itertools import combinations

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
    def __init__(self, x, y, gene, mutation, variant_class, id, p_key,
                 emcee=False):
        self.x = x
        self.y = y
        self.mutation = mutation
        self.variant_class = variant_class
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


def exponential_fit(self, neutral_cells=10000):
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

    p.add('neutral_cells', value=neutral_cells, min=1)
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


""" Exponential model fitted in bulk by participant"""


def participant_fit(part, init, rate=1.3, fit_method='least_squares',
                    vary_cells=True):
    """ Fits an exponential model of clonal growth in bulk for all
    fit trajectories in a participant. Trajectories share a
    neutral_cell parameter, and the total number of cells is updated
    with the growth of each trajectory.
    Inputs:
    - part: List. List of all trajectories in a participant.
    - init: initialization list. Contains 2 sublists and one element.
            init[0]: List of initialisations for all fitness parameters.
            init[1]: List of initialisation for all origin parameters.
            init[2]: Float initialising neutral_cells parameter.

    The following attributes are added for each trajectory in part:
    - nelder: optimallmfit Nelder model.
    - fitness, origin and neutral cells parameters.
    - data_vaf and model_vaf: dictionaries containing data and model
      predictions as a function of time.
    - vaf_plot: dictionary contianing model predictions for plotting.

    Output:
        """

    fitness_init = init[0]
    origin_init = init[1]
    neutral_cells_init = init[2]
    # Initialize model parameters
    p = Parameters()
    for i, traj in enumerate(part):
        p.add('fitness_%i' % i, value=fitness_init[i],
              min=0, max=1.3)
        p.add('origin_%i' % i, value=origin_init[i],
              min=0, max=traj.x[0])
        p.add('neutral_cells_%i' % i, value=neutral_cells_init,
              min=500, max=200_000, vary=vary_cells)

    # force all neutral_cells parameters equal
    if len(part) > 1:
        for i in range(1, len(part)):
            p['neutral_cells_%i' % i].expr = 'neutral_cells_0'

    # fit the data using lmfit 'minimize' function
    model = Minimizer(joint_exponential_residual, p, fcn_args=(part,))
    fit = model.minimize(method=fit_method)

    # update trajetories
    for i, traj in enumerate(part):
        traj.fitness = fit.params['fitness_%i' % i]*1
        traj.fitness_percentage = 100*traj.fitness / rate
        traj.neutral_cells = int(fit.params['neutral_cells_%i' % i]*1)
        traj.origin = fit.params['origin_%i' % i]*1
        traj.fit = fit
        data_vaf = np.array(list(traj.data_vaf.values()))
        model_vaf = np.array(list(traj.model_vaf.values()))
        traj.r2 = np.linalg.norm(data_vaf-model_vaf)

    # update residual
    residual = sum([traj.r2 for traj in part])

    time = np.linspace(65, 95, 1000)
    for traj in part:
        traj.vaf_plot = dict.fromkeys(time)

    for x in time:
        traj_cells = []
        for traj in part:
            traj_cells.append(exponential_cells(x, traj.fitness, traj.origin))
        total_cells = part[0].neutral_cells + sum(traj_cells)
        for i, traj in enumerate(part):
            traj.vaf_plot[x] = traj_cells[i]/(2*total_cells)

    return fit, residual


def joint_exponential_residual(p, participant):
    """ Computes a residual for minimization.
    For each timepoint:
    Step 1 - compute the number of cells at timepoint in each trajectory:
             traj_cells
    Step 2 - update total number of cells in the participant at time point as
             total_cells = neutral_cells + cells in each trajectory.
    Step 3 - compute the VAF for each trajectory at timepoint as
             traj_vaf = traj_cells / (2*total_cells)
    returns: sum of squared errors for each data point in each trajectory.
    """

    # Initialize dictionaries for data and model_predictions
    # in each time point
    for traj in participant:
        traj.data_vaf = dict(zip(traj.x, traj.y))
        traj.model_vaf = dict.fromkeys(traj.x)

    # Extract all fitted time_points in trajectories
    # Note, there can be time_points present only i one trajectory
    time_points = []
    for traj in participant:
        for x in traj.x:
            time_points.append(x)
    # remove duplicate values
    time_points = list(set(time_points))

    # compute VAF of trajectories at each time_point
    for time in time_points:
        # Create a list of the number of cells in each variant at time x
        traj_cells = []
        for i, traj in enumerate(participant):
            traj_cells.append(
                exponential_cells(time,
                                  p['fitness_%i' % i],
                                  p['origin_%i' % i]))

        # compute the total number of cells present in the system
        # note that we take the integer of parameter neutral cells
        total_cells = int(p['neutral_cells_0']) + sum(traj_cells)

        # compute the vaf of each trajectory at time x
        traj_vaf = np.array(traj_cells) / (2*total_cells)

        for i, traj in enumerate(participant):
            if time in traj.model_vaf.keys():
                traj.model_vaf[time] = traj_vaf[i]

    # compute the error for each trajectory on each time_point
    residual = []
    for traj in participant:
        for key in traj.data_vaf.keys():
            residual.append(abs(traj.data_vaf[key] - traj.model_vaf[key]))

    # return the array that needs minimizing
    return residual


def exponential_cells(time, fitness, origin):
    """ Returns the number of cell of an exponentially growing trajectory.
    The function inputs are
    time: time point for evaluation.
    fitness: growth speed of the family of cells.
    origin: time of origin of the family of cells."""
    return np.exp(fitness*(time-origin))


def init_fit(part, method='least_squares', n_iterations=500):
    '''Runs modelling.participant_fit with a grid of
    initialisation values for all model parameters.
    Selects best fit model based on total R2 error.
    Inputs:
    - part: List. List of trajectories in a participant.
    - method: String. Any minmization algorithm as implemented in lmfit.
    Returns:
    - part_id: String. Participant's id for fit.
    - part: List. List of trajectories in fitted participant.
    - plot_data: Pandas DataFrame. History of all fits using different
                                   initialisations.
    - '''

    # Set number of initialisation iterations
    iterations = len(part)*500

    # Create parameter initialisation list
    init_list = []
    for i in range(iterations):
        # For each iteration append parameters for fit initialisation
        fitness_init = np.random.uniform(0, 1.3, len(part))
        origin_init = [np.random.uniform(0, traj.x[0]) for traj in part]
        neutral_cells_init = np.random.uniform(500, 200_000)
        init_list.append([fitness_init, origin_init, neutral_cells_init])

    # Create fitted fit and residual list to track fitted models
    # and their respective R2 error.
    fit_list = []
    residual_list = []

    # Fit model using all parameter initialisations in init_list.
    for init in init_list:
        fit, residual = participant_fit(part=part, init=init,
                                        fit_method=method)

        # append  fit and residual to lists.
        fit_list.append(fit)
        residual_list.append(residual)

    # plot optimal trajectory
    min_residual = 10
    for residual, init in zip(residual_list, init_list):
        if residual < min_residual:
            min_residual = residual
            min_init = init

    optimal_fit, optimal_residual = participant_fit(
                                init=min_init, part=part, fit_method=method)

    # Create dataframe for plotting distributions of parameter fits
    # by initialisation.
    plot_data = pd.DataFrame(columns=['Gene', 'Gene name', 'Fitness',
                                      'Origin', 'Cells', 'error', 'aic'])
    for fit, residual in zip(fit_list, residual_list):
        for i in range(0, len(part)):
            plot_data = plot_data.append(
                            {'Gene': i + np.random.normal(0, 0.1),
                             'Gene name': part[i].mutation,
                             'Fitness': fit.params['fitness_%i' % i]*1,
                             'Origin': fit.params['origin_%i' % i]*1,
                             'Cells': fit.params[
                             'neutral_cells_%i' % i]*1,
                             'error': residual,
                             'aic': fit.aic}, ignore_index=True)

    # print participant id to keep track of progress
    print(part[0].id)

    return part[0].id, part, plot_data, optimal_fit, optimal_residual


def fit_postprocessing(init_fit_output):
    '''Post-processing of 'init_fit' output when function is mapped to
    a list of participants.
    Input:
    list_participants: List. List of outputs after
                             mapping a list of participants to init_fit.
    Returns:
    model: List. Flat list of all fitted trajectories.
    fit_dict: Dictionary.
        keys: String. Participant id.
        items: List. Fit history DataFrame, Lmfit optimal model,
                     optimal fit residual.
    '''
    model = []
    fit_dict = dict()
    for fit in init_fit_output:
        fit_dict[fit[0]] = fit[2:]
        model.append(fit[1])

    model = [item for sublist in model for item in sublist]

    return model, fit_dict


def fitting_error_plot(key, part_dict, fit_dict, filter_quantile=0.1):
    '''Creates a plot of all obtained fitted parameters by gene,
    coloured according to the R2 error.
    Inputs:
    - key: part_id
    - part_dict: dictionary of all participants trajectories.
    - fit_dict: dictionary of all participants fit histories.
    - filter_quantile: Float. Quantile filtering threshold.
    Returns:
    fig: plotly graphical object.'''

    # Extract participants trajectories, and fit initialisation history.
    part = part_dict[key]
    fit = fit_dict[key]

    # Extract DataFrame from particpant initialisation history
    data = fit[0]

    # Filter data using filter_quantile
    if filter_quantile is not None:
        data = data[data.error < np.quantile(data.error, filter_quantile)]

    # plot optimal fitness
    fig = px.scatter(data, x='Gene', y='Fitness',
                     color='error', hover_data=['aic'])
    tickvals = []
    ticktext = []
    for i, traj in enumerate(part):
        ticktext.append(traj.gene)
        tickvals.append(i)
    # Set custom x-axis labels
    fig.update_xaxes(title=None, ticktext=ticktext, tickvals=tickvals)

    return fig


def cells_origin(fit_dict, part_id, filter_quantile=0.1):
    '''Returns scatter plot of origin vs neutral cells for all trajectories in
    participant.'''

    data = fit_dict[part_id][0]

    # Filter data using filter_quantile
    if filter_quantile is not None:
        data = data[data.error < np.quantile(data.error, filter_quantile)]

    # Plot the relation between origin and cells for the first variant
    fig = px.scatter(data, x='Origin', y='Cells', color='Gene name')

    return fig


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
                      and model.fit.params['fitness'] > 0
                      and model.fit.aic < -10]

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
    y_label = prediction(traj.fit.params, x_label)

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
        y_fit = prediction(traj.fit.params, x_line)

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

        fitness.append(traj.fit.params['fitness'])  # append fitness to list

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
                      and model.fit.params['fitness'] > 0]

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
    y_label = prediction(traj.fit.params, x_label)

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
        y_fit = prediction(traj.fit.params, x_line)

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

        fitness.append(traj.fit.params['fitness'])  # append fitness to list

    # box plot with fitness estimates
    # Frist create the box
    fig.add_trace(go.Box(y=fitness, boxpoints=False,
                         name=mutation, marker_color=colors[0],
                         showlegend=False),
                  row=1, col=2)

    # Seccond create the dots adding jitter
    for i, item in enumerate(fitness):
        fig.add_trace(go.Box(y=[item], pointpos=-2, boxpoints='all',
                             name=mutation,
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
                       marker_color=colors[i % 10],
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
                       marker_color=colors[i % 10],
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


def extended_plot_variant(model_list, mutation, var_dict=None):
    if var_dict is None:
        var_types = ["3'Flank", "5'Flank", "5'UTR", "Frame_Shift_Del",
                     "Frame_Shift_Ins", "In_Frame_Ins", "Missense_Mutation",
                     "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region",
                     "Splice_Site"]
        var_colors = ["#00BBDA", "#E18A00", "#BE9C00", "#8CAB00", "#24B700",
                      "#00BE70", "#00C1AB", "#F8766D", "#8B93FF", "#D575FE",
                      "#F962DD"]

        var_zip = zip(var_types, var_colors)
        var_dict = dict(var_zip)

    model_filtered = [traj for traj in model_list
                      if traj.gene == mutation
                      and traj.fit.params['fitness'] > 0]

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
    y_label = prediction(traj.fit.params, x_label)

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
    fitness_variant = []
    # Create a trace for each traject
    for i, traj in enumerate(model_filtered):
        x_line = np.linspace(65, 95, 100)
        y_fit = prediction(traj.fit.params, x_line)
        fitness_variant.append(traj.variant_class)
        fig.add_trace(
            go.Scatter(x=traj.x, y=traj.y,
                       name=traj.mutation,
                       legendgroup="group2",
                       showlegend=False,
                       mode='markers',
                       line=dict(color=var_dict[traj.variant_class])),
            row=1, col=1)

        fig.add_trace(
            go.Scatter(x=x_line, y=y_fit,
                       name='fit', legendgroup='group2',
                       showlegend=False,
                       line=dict(color=var_dict[traj.variant_class])),
            row=1, col=1)
        fitness.append(traj.fit.params['fitness'])  # append fitness to list

    # box plot with fitness estimates
    # Frist create the box
    fig.add_trace(
        go.Box(y=fitness, boxpoints=False,
               name=mutation, marker_color=colors[0],
               showlegend=False),
        row=1, col=2)

    # Seccond create the dots adding jitter
    for i, item in enumerate(fitness):
        fig.add_trace(
            go.Box(y=[item], pointpos=-2, boxpoints='all',
                   name=mutation,
                   jitter=random.uniform(0, 0.5),
                   marker_color=var_dict[fitness_variant[i]],
                   line=dict(color='rgba(0,0,0,0)'),
                   fillcolor='rgba(0,0,0,0)',
                   showlegend=False),
            row=1, col=2)

    fitness_variant = list(set(fitness_variant))
    for key in fitness_variant:
        for traj in model_filtered:
            if traj.variant_class == key:
                fig.add_trace(
                    go.Scatter(x=traj.x,
                               y=traj.y,
                               name=key,
                               mode='markers',
                               line=dict(color=var_dict[key]),
                               legendgroup="group2"))
                break
    fig.update_layout(title=f'{mutation} trajectories')
    fig.update_xaxes(title_text='Age (in years)', row=1, col=1)
    fig.update_yaxes(title_text='VAF', row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=1, col=2)
    fig.update_yaxes(title_text='Fitness', row=1, col=2)

    return fig


def gene_trajectories(model_list, gene):

    var_types = ["3'Flank", "5'Flank", "5'UTR", "Frame_Shift_Del",
                 "Frame_Shift_Ins", "In_Frame_Ins", "Missense_Mutation",
                 "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region",
                 "Splice_Site"]
    var_colors = ["#00BBDA", "#E18A00", "#BE9C00", "#8CAB00", "#24B700",
                  "#00BE70", "#00C1AB", "#F8766D", "#8B93FF", "#D575FE",
                  "#F962DD"]

    var_zip = zip(var_types, var_colors)
    var_dict = dict(var_zip)

    fig = make_subplots(rows=1, cols=2,
                        column_widths=[0.7, 0.3],
                        subplot_titles=(f'{gene} trajectories',
                                        'Fitness distribution'))
    fitness_color = []
    fitness = []
    for traj in model_list:
        if traj.gene == gene:
            x = list(traj.data_vaf.keys())
            y = list(traj.data_vaf.values())
            fig.add_trace(
                go.Scatter(x=x, y=y,
                           mode='markers',
                           marker_color=var_dict[traj.variant_class],
                           name=traj.mutation,
                           showlegend=False),
                row=1, col=1)
            x_prediction = list(traj.vaf_plot.keys())
            y_prediction = list(traj.vaf_plot.values())
            fig.add_trace(
                go.Scatter(x=x_prediction,
                           y=y_prediction, mode='lines',
                           marker_color=var_dict[traj.variant_class],
                           showlegend=False,
                           text=(f'fitness: {round(traj.fitness,3)}<br>'
                                 f'origin: {round(traj.origin,3)}<br>'
                                 f'r2: {round(traj.r2,3)}'),
                           name=traj.mutation),
                row=1, col=1)
            fitness.append(traj.fitness)
            fitness_color.append(var_dict[traj.variant_class])

    # box plot with fitness estimates
    # Frist create the box
    fig.add_trace(
        go.Box(y=fitness, boxpoints=False,
               name=gene,
               marker_color='Grey',
               showlegend=False),
        row=1, col=2)

    # Seccond create the dots adding jitter
    for i, item in enumerate(fitness):
        fig.add_trace(
            go.Box(y=[item], pointpos=-2, boxpoints='all',
                   name=gene,
                   jitter=random.uniform(0, 0.2),
                   marker_color=fitness_color[i],
                   line=dict(color='rgba(0,0,0,0)'),
                   fillcolor='rgba(0,0,0,0)',
                   showlegend=False),
            row=1, col=2)

    fig.update_xaxes(title_text='Age (in years)', row=1, col=1)
    fig.update_yaxes(title_text='VAF', row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=1, col=2)
    fig.update_yaxes(title_text='Fitness', row=1, col=2)

    return fig


def gene_box(cohort, order='median', percentage=False):
    """Box plot with counts of filtered mutations by gene. Returns a figure."""

    # Create a dictionary with all filtered genes
    gene_list = []
    for traj in cohort:
        gene_list.append(traj.gene)
    gene_dict = {element: [] for element in set(gene_list)}

    # update the counts for each gene
    if percentage is False:
        y_label = 'Fitness'
        for traj in cohort:
            fitness = traj.fitness
            gene_dict[traj.gene].append(fitness)
    if percentage is True:
        y_label = 'fitness_percentage'
        for traj in cohort:
            fitness = traj.fitness_percentage
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
    color_dict = dict()
    for i, key in enumerate(gene_dict):
        color_dict[key] = colors[i % 10]
        fig.add_trace(
            go.Box(y=gene_dict[key],
                   marker_color=color_dict[key],
                   name=key, boxpoints='all', showlegend=False))

    fig.update_layout(title='Gene distribution of filtered mutations',
                      yaxis_title=y_label,
                      template="simple_white")
    if percentage is False:
        fig.update_yaxes(type='log', tickvals=[0.05, 0.1, 0.2, 0.4])
    fig.update_layout(xaxis_tickangle=-45)
    return fig, gene_dict, color_dict


def gene_statistic(gene_dict, statistic='kruskal-walllis', filter=True):
    """ compute a statistical test to find significant differences in the
    distribution of fitness by gene.
    statistic parameter accepts: 'kruskal' or 'anova'.
    Returns:
    * heatmap with significant statistical differences.
    * dataframe."""
    # extract all possible gene combinations
    gene_list = []
    for gene in gene_dict.keys():
        if len(gene_dict[gene]) > 2:
            gene_list.append(gene)

    # Create dataframe to store statistics
    test_df = pd.DataFrame(index=gene_list, columns=gene_list)
    for gene1, gene2 in combinations(gene_list, 2):
        # compute statistic for each possible comination of genes
        if statistic == 'kruskal-wallis':
            stat, pvalue = stats.kruskal(gene_dict[gene1], gene_dict[gene2])
        if statistic == 'anova':
            stat, pvalue = stats.f_oneway(gene_dict[gene1], gene_dict[gene2])
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

    return fig, test_df


def damage_class(model, fs_ter='fs - ter',
                 color_list=None, plot_color=colors[0]):
    """ Box plot of variant predicted protein damage class ~ Fitness.
    Returns:
    * Box plot
    * Modified model with damage_class attribute."""

    # Load dataset with prediction damage
    xls = pd.ExcelFile('Datasets/new_kristina_variants.xlsx')
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
            traj.damage_class = fs_ter
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
                     marker_color=plot_color,
                     showlegend=False),
        row=1, col=1)

    fig.add_trace(
        go.Box(x=df['damage_class'],
               y=df['fitness'],
               marker_color=plot_color,
               boxpoints=False,
               showlegend=False),
        row=2, col=1)

    fig.update_layout(
        template='simple_white',
        boxmode='group'
    )
    fig.update_yaxes(title='fitness', type='log', dtick=[0.05,0.1,0.2,0.4], row=2, col=1)
    fig.update_yaxes(title='trajectory counts', row=1, col=1)
    return fig, df
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
