# Last modification: E. Latorre -- 4-11-2020

# Module containing a model class definition and functions to fit exponential
# and logistic models of growth.

import numpy as np
from lmfit import Parameters, Minimizer, minimize
from lmfit.models import GaussianModel

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

    def fit(self, model='logistic', reg=0.1, common_fit=False):
        self.model = model           # set model attribute
        if model == 'logistic':
            self.out = logistic_model(self, reg, common_fit)
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


class model_emcee:
    # mutation trajectory class for training mp models
    def __init__(self, x, y, gene, mutation, id, damage=None,
                 damaging_class=None, known=None, out=None):
        self.x = x
        self.y = y
        self.mutation = mutation
        self.gene = gene
        self.id = id
        self.damage = damage
        self.damaging_class = damaging_class
        self.known = known
        self.out = out


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

                # append an entry for each trajectory
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


def logistic(x, amplitude, fitness, origin, N0=0.004):
    """logistic lineshape."""
    sig_center = origin + (1/fitness) * np.log((amplitude-N0)/N0)
    return amplitude / (1. + np.exp(- fitness * (x - sig_center)))


def logistic_dataset(params, i, x):
    """Calculate logistic lineshape from parameters for data set."""
    amp = params['amp_%i' % (i+1)]
    fit = params['fit_%i' % (i+1)]
    dis = params['dis_%i' % (i+1)]
    return logistic(x, amp, fit, dis)


def logistic_model(self, l2, common_fit):
    # l2 = weight of the regularization of amplitude parameter around 0.5
    # initialize parameters using default values and boundary values
    fit_params = Parameters()
    for iy, y in enumerate(self.y):
        fit_params.add('amp_%i' % (iy+1),
                       value=0.4, min=0.2, max=0.5)
        fit_params.add('fit_%i' % (iy+1), value=0)
        fit_params.add('dis_%i' % (iy+1), value=0)

    if common_fit is True:
        # force all fit parameters equal
        for iy in range(2, len(self.y)+1):
            fit_params['fit_%i' % iy].expr = 'fit_1'

    # fit the data using lmfit 'minimize' function
    return minimize(logistic_objective, fit_params, args=(self.x, self.y, l2))


def logistic_objective(params, x, y, a):
    """Calculate total residual for fits of logistics to several data sets."""

    total_res = []   # keeps track of the cumulative residual

    # compute the residual for each trajectory in the dataset
    for i, points in enumerate(zip(x, y)):
        res = points[1] - logistic_dataset(params, i, points[0])
        total_res.append(res)

    # each entry of total_res is the list of residuals for each trajectory
    # as minimize() needs a 1D array -> flatten
    total_res = np.concatenate(total_res)
    for i, j in enumerate(y):
        total_res = np.append(total_res,
                              np.asarray([a*(0.5-params['amp_%i' % (i+1)])]))

    return total_res


""" Logistic model emcee"""


def residual(p, x, y, initial_pop=0.004):
    """Calculate total residual for fits of logistics to several data sets."""
    v = p.valuesdict()
    res = y - logistic(x, v['amplitude'],
                       v['fitness'], v['origin'], initial_pop)
    res = np.append(res, v['regularization']*(0.5-v['amplitude']))

    return res


def logistic_emcee(self, emcee=True, regularization=0.05):
    # set origin_min as participants age.
    if 'LBC0' in self.id:
        origin_min = - 79
    else:
        origin_min = -70

    # Create model parameters
    p = Parameters()
    p.add('amplitude', value=0.4, min=0.2, max=0.5)
    p.add('fitness', value=0.1, min=0)
    p.add('origin', value=-30, min=origin_min, max=0)
    p.add('regularization', value=regularization, vary=False)

    # step 1: Crete Minimizer object model
    self.model = Minimizer(residual, params=p,
                           fcn_args=(self.x, self.y))
    # step 2: minimize mini with Nelder-Mead algorithm (robust)
    self.nelder = self.model.minimize(method='nelder')
    # step 3: fine tune minimization with Levenberg-Marquardt algorithm ()
    # This allows the computation of confidence intervals
    self.ls = self.model.minimize(method='leastsq',
                                  params=self.nelder.params.copy())

    if emcee is True:
        # add sigma parameter (~exp(variance) of the error distribution)
        #self.nelder.params.add('__lnsigma',
        #                       value=np.log(0.1),
        #                       min=np.log(0.001), max=np.log(2))
        self.emcee = minimize(residual, method='emcee', seed=1,
                              nan_policy='omit', burn=300, steps=1000,
                              nwalkers=100, thin=1,
                              params=self.nelder.params,
                              args=(self.x, self.y),
                              is_weighted=False, progress=False)

    return self


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


""" Gaussian model"""


def gauss_fit(data, center=None, bin_number=200):
    # Distribution - bin the data to create a histogram:
    # x = bins, y = counts
    y, x = np.histogram(data, bins=bin_number,
                        range=(min(data), max(data)),
                        density=True)
    y = y/sum(y)
    x = 0.5 * (x[:-1] + x[1:])

    # Fitting a normal distribution using lmfit
    # Design the model
    gauss1 = GaussianModel(prefix='g1_')
    pars = gauss1.guess(data=y, x=x)
    if center is not None:
        pars['g1_center'].value = center
        pars['g1_center'].vary = False
    mod = gauss1

    # Fit the model
    out = mod.fit(y, pars, x=x)

    # Model components
    comps = out.eval_components(x=x)

    # Scatter plot
    fig = go.Figure()

    fig.add_trace(go.Bar(x=x, y=y, marker_color=colors[0],
                         name='total gradients'))
    fig.add_trace(go.Scatter(x=x, y=comps['g1_'],
                             mode='lines', name='Fit'))

    fig.update_layout(height=600, width=800,
                      title_text="Fitted Gaussian distribution")

    return out, fig


""" Auxiliary plotting functions"""


def fit_plot(self, min_x, max_x, model):
    # round up fitness
    fitness = [self.out.params[s].value for s in self.out.params if "fit" in s]
    fitness = np.mean(fitness)
    chi = round(self.out.chisqr*1, 3)

    # define a general evaluate function for both exponential and logistic
    if model == 'logistic':
        evaluate = logistic_dataset
    elif model == 'exponential':
        evaluate = exponential_dataset

    fig = go.Figure()
    x_line = np.linspace(min_x, max_x, 1000)

    # Plot the first trajectory to create labels
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


def mutation_plot(trajectories, mutation):
    # plot fitted trajectories containing a mutation

    fig = go.Figure()
    x_line = np.linspace(-5, 15, 1000)
    i = 0  # color counter
    for traj in trajectories:
        if traj.mutation == mutation or \
           traj.mutation.split()[0] == mutation:

            y_fit = logistic(x_line, traj.nelder.params['amplitude'],
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
    fig.update_xaxes(title='time in years')
    fig.update_yaxes(title='AF')
    return fig
