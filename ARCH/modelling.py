# Last modification: E. Latorre -- 29-10-2020


import numpy as np
from lmfit import Parameters, minimize
import plotly.figure_factory as ff
import plotly.express as px
import plotly.graph_objects as go
colors = px.colors.qualitative.Plotly


class model:
    def __init__(self, cohort, gene, l_filter, h_filter, min_points=3,
                 trajectories=None, participants=None):
        self.cohort = cohort
        self.gene = gene
        self.l_filter = l_filter
        self.h_filter = h_filter
        self.min_points = min_points
        self.x, self.y, self.participants, \
            self.mutation_list, self.trajectories \
            = preprocess(cohort, gene=self.gene,
                         l_filter=self.l_filter, h_filter=self.h_filter,
                         min_points=self.min_points)

    def fit(self, model='logistic', l2=0.1, common_fit=False):
        self.model = model           # set model attribute
        if model == 'logistic':
            self.out = logistic_model(self, l2, common_fit)
        elif model == 'exponential':
            self.out = exponential_model(self, common_fit)
        else:
            print('model name not recognized')

        # create the list of all estimated fitness parameters
        self.fitness = [self.out.params[s].value for s in self.out.params
                        if "fit" in s]

    def fitness_plot(self, bin_size=.01):
        return fitness_distribution(self, bin_size)

    def plot(self, min_x=0, max_x=10):
        return fit_plot(self, min_x, max_x, model=self.model)


def preprocess(cohort, gene, l_filter, h_filter, min_points):
    # set the trajectories that we want to fit

    y = []                # initializing list of y_axis data
    x = []                # initializing list of x_axis data
    participant = []      # initializing list participants
    mut = []
    for part in cohort:
        for traj in part.trajectories:
            if gene in traj.mutation.split() and \
               l_filter <= traj.gradient <= h_filter and \
               len(traj.data) >= min_points and \
               traj.germline is False:

                y.append(traj.data.AF)
                x.append(traj.data.wave)
                participant.append(part.id)
                mut.append(traj.mutation)

    # plot trajectories preprocessed for fitting
    fig = go.Figure()
    for i, points in enumerate(zip(x, y)):
        fig.add_trace(go.Scatter(x=points[0],
                                 y=points[1],
                                 mode='lines+markers',
                                 name=mut[i]))
    # Edit the layout
    fig.update_layout(title=f'Mutated {gene} trajectories for fitting',
                      xaxis_title='years',
                      yaxis_title='VAF')

    return x, y, participant, mut, fig


""" Logistic model fit functions"""


def logistic_model(self, l2, common_fit):
    # l2 = weight of the regularization of amplitude parameter around 0.5
    # initialize parameters using default values and boundary values
    fit_params = Parameters()
    for iy, y in enumerate(self.y):
        fit_params.add('amp_%i' % (iy+1),
                       value=0.4, min=0.2, max=0.5)
        fit_params.add('fit_%i' % (iy+1), value=0)
        fit_params.add('dis_%i' % (iy+1), value=0,   min=-25,  max=25)

    if common_fit is True:
        # force all fit parameters equal
        for iy in range(2, len(self.y)+1):
            fit_params['fit_%i' % iy].expr = 'fit_1'

    # fit the data using lmfit 'minimize' function
    return minimize(logistic_objective, fit_params, args=(self.x, self.y, l2))


def logistic(x, amp, fit, dis):
    """logistic lineshape."""
    return amp / (1.+np.exp(-fit*(x-dis)))


def logistic_dataset(params, i, x):
    """Calculate logistic lineshape from parameters for data set."""
    amp = params['amp_%i' % (i+1)]
    fit = params['fit_%i' % (i+1)]
    dis = params['dis_%i' % (i+1)]
    return logistic(x, amp, fit, dis)


def logistic_objective(params, x, data, a):
    """Calculate total residual for fits of logistics to several data sets."""

    total_res = []   # keeps track of the cumulative residual

    # compute the residual for each trajectory in the dataset
    for i, points in enumerate(zip(x, data)):
        res = points[1] - logistic_dataset(params, i, points[0])
        total_res.append(res)

    # now flatten this to a 1D array, as minimize() needs
    # return np.concatenate(total_res).ravel()
    total_res = np.concatenate(total_res).ravel()
    for i, j in enumerate(data):
        total_res = np.append(total_res,
                              np.asarray([a*(0.5-params['amp_%i' % (i+1)])]))

    return total_res


""" Exponential model"""


def exponential_model(self, common_fit):
    # initialize parameters
    fit_params = Parameters()
    for iy, y in enumerate(self.y):
        fit_params.add('fit_%i' % (iy+1), value=1)
        fit_params.add('dis_%i' % (iy+1), value=0,   min=-25,  max=25)

    if common_fit is True:
        # force all fit parameters equal
        for iy in range(2, len(self.y)+1):
            fit_params['fit_%i' % iy].expr = 'fit_1'

    # fit the data using lmfit 'minimize' function
    return minimize(exp_objective, fit_params, args=(self.x, self.y))


def exponential(x, fit, dis):
    """Exponential lineshape."""
    return np.exp(fit*(x-dis))


def exponential_dataset(params, i, x):
    """Calculate exponential lineshape from parameters for data set."""
    fit = params['fit_%i' % (i+1)]
    dis = params['dis_%i' % (i+1)]
    return exponential(x, fit, dis)


def exp_objective(params, x, data):
    """Calculate total residual for fits of exponents to several data sets."""

    total_res = []    # cumulative residual

    # compute the residual for each trajectory in the dataset
    for i, points in enumerate(zip(x, data)):
        res = points[1] - exponential_dataset(params, i, points[0])
        total_res.append(res)

    total_res = np.concatenate(total_res).ravel()

    return total_res


""" Auxiliary plotting functions"""


def fit_plot(self, min_x, max_x, model):
    # round up fitness
    fitness = [self.out.params[s].value for s in self.out.params if "fit" in s]
    fitness = np.mean(fitness)
    chi = round(self.out.chisqr*1, 3)
    if model == 'logistic':
        evaluate = logistic_dataset
    elif model == 'exponential':
        evaluate = exponential_dataset

    fig = go.Figure()
    x_line = np.linspace(min_x, max_x, 1000)

    y_fit = evaluate(self.out.params, 0, x_line)

    fig.add_trace(go.Scatter(x=x_line, y=y_fit,
                             mode='lines', line=dict(color='#222A2A'),
                             legendgroup="group",
                             name="Fitted trajectories"))

    fig.add_trace(go.Scatter(x=self.x[0], y=self.y[0],
                             mode='markers', line=dict(color='#222A2A'),
                             legendgroup="group",
                             name="Data points"))

    for i, points in enumerate(zip(self.x, self.y)):
        x_line = np.linspace(min_x, max_x, 1000)
        y_fit = evaluate(self.out.params, i, x_line)

        fig.add_trace(go.Scatter(x=points[0], y=points[1], mode='markers',
                                 line=dict(color=colors[i % 10]),
                                 legendgroup="group2",
                                 name=f'{self.participants[i]}'))

        fig.add_trace(go.Scatter(x=x_line, y=y_fit, mode='lines',
                                 line=dict(color=colors[i % 10]),
                                 legendgroup="group2",
                                 name='Fitted trajectories',
                                 showlegend=False))

        # Edit the layout
        fig.update_layout(title=f'{model} fit of {self.gene} mutations <br>'
                                f'total squared-error: {chi}, '
                                f'fitness={fitness}',
                                xaxis_title='years',
                                yaxis_title='VAF')
    return fig


def fitness_distribution(self, bin_size):
    hist_data = [self.fitness]
    group_labels = [self.gene]

    fig = ff.create_distplot(hist_data, group_labels,
                             curve_type='normal',
                             show_hist=False,
                             bin_size=bin_size)
    # Add title
    fig.update_layout(title='Distribution of estimated fitness parameters')
    return fig
