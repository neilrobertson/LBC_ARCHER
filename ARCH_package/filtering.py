# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================
"""Module to create a filtering of a cohort of participants based on their
genetic trajectories.
"""
# =============================================================================
#region: Import packages
# =============================================================================
# Import local packages
# =============================================================================
# Append root directory to system's path

import sys
sys.path.append('../ARCH_package')

import auxiliary
import distributions
import modelling
import plot

# =============================================================================
# Import python packages
# =============================================================================
import pandas as pd
import numpy as np


from scipy.optimize import minimize
from scipy.stats import beta, betabinom, binom, nbinom
from sklearn.linear_model import LinearRegression
from numpy import concatenate, sort

import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

from functools import partial
from itertools import combinations, product
from multiprocessing import Pool
from tqdm import tqdm

colors = px.colors.qualitative.Plotly
pio.templates.default = "simple_white"

#endregion
# =============================================================================
#region: Define local classes
# =============================================================================


class filter():
    """Class corresponding to a filtering of the full cohort
    Attributes:
    - model: List. Contains a list of model class instances following
                   trajectories that pass the filter. This is used
                   for later model fitting.
    - gradient_plot: Plotly figure. Figure showing data VAF (x-axis) and
                     regularized gradient between first and last time-point
                     (y-axis). Data points are colored according to wether they
                     pass the filter and are selected for model fitting.
    - gene_bar: Plotly figure. Bar plot of trajectory counts by gene passing
                the filter.
    - gene_dict: Dictionary. Counts of trajectories passing the filter by gene.
                     """

    def __init__(self, model=None, gradient_plot=None,
                 neutral_dist=None, gene_bar=None, gene_dict=None):
        self.model = model
        self.gradient_plot = gradient_plot
        self.gene_bar = gene_bar
        self.gene_dict = gene_dict


class neutral_fit():
    """ Class containing linear models fitted to VAF means and variance.
    Attributes:
    - mean_regr: scikit linear model. Linear regression of VAF gradient means.
    - var_regr: scikit linear model. Linear regression of VAF gradient
                variance.
    - figures: Plotly figures."""

    def __init__(self, mean_regr=None, var_regr=None, figures=None):
        self.mean_regr = mean_regr
        self.var_regr = var_regr
        self.figures = figures


class trajectory_model():
    """Model the likelihood of observing a time-series clonal trajectory.

    Attributes
    -----------
    model_type: str. Model used to describe the data.
    fit. scipy.OptimizeResult. Optimized fit of the negative log likelihood
         of the model.
    data. Pandas DataFrame. DataFrame containing the trajectorie's time series
          data.

    Methods:
    compute_IC: Add information criterion attribute to the class object.
    """

    def __init__(self, model_type=None, fit=None, data=None, nll=None, likelihood=None, n_params=None):
        self.model_type = model_type
        self.fit = fit
        self.data = data
        self.nll = nll
        self.likelihood = likelihood
        self.n_params = n_params

    def compute_IC(self, n_params=1):
        """Add information criterion attribute to the class object.
        Parameters
        -----------
        n_params: int. Override the number of parameters in the model.
        """
        self.IC = information_criterion(self, n_params=n_params)

    def plots(self):
        if self.model_type == 'bd_process':
            # self.heatmap = likelihood_heatmap_N(self.brute_force_df)
            self.posterior_dist = conditional_distributions(
                                                self.brute_force_df)
        elif self.model_type == 'binomial_artifact':
            self.posterior_plot = px.line(x=self.prior[0], y=self.posterior)
            self.prior_plot = px.line(x=self.prior[0], y=self.prior[1])


class model_comparison ():
    """class object for model comparison."""
    def __init__ (self, id, data, betabinom_artifact_prob=None,
                  bd_prob=None, binom_artifact_prob=None,
                  bayes_factor=None, optimal_model=None):
        self.id = id
        self.data = data
        self.betabinom_artifact_prob = betabinom_artifact_prob
        self.bd_prob = bd_prob
        self.binom_artifact_prob = binom_artifact_prob
        self.bayes_factor = bayes_factor
        self.optimal_model = optimal_model

    def plot(self):
        """ Plot mutational trajectories associated with mutation."""
        fig = go.Figure()
        for traj in self.data:
            fig.add_trace(
                go.Scatter(
                    x=traj.age,
                    y=traj.AF
                )
            )
        fig.update_layout(title=self.id)
        return fig


    def compute_optimal_model(self, binomial=True, bayes_factor_threshold=100):
        """Compute Bayes factor and set optimal model."""
        if binomial is True:
            # compute bayes factor using binomial sequencing artifact model
            # for mutations with only 1 trajectory
            if len(self.data) > 1:
                self.bayes_factor = self.bd_prob/self.betabinom_artifact_prob
            else:
                self.bayes_factor = self.bd_prob/self.binom_artifact_prob
        else:
            # compute bayes factor using betabinomial sequencing artifact model
            self.bayes_factor = self.bd_prob/self.betabinom_artifact_prob

        if self.bayes_factor > bayes_factor_threshold:
            self.optimal_model = 'fit mutation'
        else:
            self.optimal_model = 'artifact'


#endregion
# =============================================================================
#region: Common filter functions
# =============================================================================


def model_cohort(cohort):
    """Create a list of model class object trajectories for fitting."""
    model_traj = []
    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                model_traj.append(
                    modelling.model(x=traj.data.age.tolist(),
                                    y=traj.data.AF.tolist(),
                                    mutation=traj.mutation,
                                    variant_class=traj.variant_class,
                                    gene=traj.mutation.split()[0],
                                    id=part.id,
                                    p_key=traj.p_key))
    return model_traj

#endregion
# =============================================================================
#region: Threshold filter functions
# =============================================================================


def threshold_filter(cohort, threshold=0.02):
    """ Filter cohort according to a VAF size threshold.
    Parameters:
    - cohort: List of participant class objects.
              List of participants to filter.
    - threshold: float. Size of VAF threshold.
                 Trajectories reaching a VAF > threshold are selected.
    Returns:
    - filter class object for model fitting.
    """

    # For each trajectory in each participant of cohort
    # set filter attribute to True if VAF reaches threshold, otherwise False.
    for part in cohort:
        for traj in part.trajectories:
            if max(traj.data.AF) >= threshold:
                traj.filter = True
            else:
                traj.filter = False

    # Create a model_cohort class using filter attribute of trajectories.
    model = model_cohort(cohort)
    # Figure with threshold filter results.
    fig = plot.CHIP_plot_inset(cohort)
    # Bar plot of genes selected through filteringfigure.
    genes = plot.gene_bar(model)

    # Create a filter class object
    filter_class = filter(model=model,
                          gradient_plot=fig,
                          gene_bar=genes[0],
                          gene_dict=genes[1])

    return filter_class

#endregion
# =============================================================================
#region: Neutral filter functions
# =============================================================================


# Set color variables for plotting
mean_color = colors[2]
var_color = colors[1]
bin_color = 'Grey'
bar_color = 'darkorange'


def neutral_filter(cohort, neutral, mean_intercept=True, n_std=2):
    """ Filter cohort using a drift induced fluctuations (DIF) filter.
    Parameters:
    - cohort: List of participant class objects.
              List of participants to filter.
    - neutral: List of participant class objects.
               List of participants' neutral trajectories.
    - mean_intercept: Boolean. Choice of using an intercept when fitting
                      a linear regression model to VAF means.
    - n_std: int. Standard deviations from the distribution of
             neutral fluctuations used for filtering.
    Returns:
    - filter class object for model fitting.
    """

    # Infer the distribution of DIF
    mean_regr, var_regr, figures = neutral_distribution(neutral,
                                                        mean_intercept,
                                                        n_std)

    # Create a neutral_fit class object storing the information of the DIF fit
    neutral_gradient_filter = neutral_fit(mean_regr=mean_regr,
                                          var_regr=var_regr,
                                          figures=figures)

    # For each trajectory in each participant of cohort set filter attribute
    # based on fluctuations between first and last timepoint.
    for part in cohort:
        for traj in part.trajectories:
            # Extract VAF at first time point
            data = traj.data['AF'].iloc[0]
            # Predict maximum expecteed fluctuation based on
            # n_std of the infered fluctuations of DIF.
            threshold = mean_regr.predict(
                        np.c_[data]) + n_std*np.sqrt(
                                            var_regr.predict(np.c_[data]))
            # set filter attribute to True if total gradient
            # exceeds the maximum expected fluctuation.
            if traj.gradient > threshold:
                traj.filter = True
            else:
                traj.filter = False

    # Create a model_cohort class using filter attribute of trajectories.
    model = model_cohort(cohort)
    # Figure with neutral filter results.
    fig = plot.neutral_plot_inset(cohort, neutral, mean_regr, var_regr, n_std)
    # Bar plot of genes selected through filteringfigure.
    genes = plot.gene_bar(model)

    # Create a filter class object
    filter_class = filter(model=model,
                          gradient_plot=fig,
                          gene_bar=genes[0],
                          gene_dict=genes[1])

    # Create neutral_dist attribute containing the output of DIF fit.
    filter_class.neutral_dist = neutral_gradient_filter
    filter_class.neutral_fit = figures

    return filter_class


def neutral_distribution(syn, mean_intercept=False, n_std=2,):
    # Create a pandas dataframe with all possible combinations of
    # trajectory timepoints
    # Step 1: Initialize dataframe.
    melt_syn = pd.DataFrame(columns=['AF', 'regularized_gradient'])

    for part in syn:
        for traj in part.trajectories:
            # For each trajectory extract possible combinations of timepoints
            for combination in list(combinations(traj.data.index, 2)):
                # create subset of 2-timepoint combination
                data = traj.data.loc[list(combination)]
                # compute gradient for the pair of data points
                gradient = np.diff(data.AF) / np.sqrt(np.diff(data.age))
                # append VAF at first timepoint and gradient to dataframe
                melt_syn = (
                    melt_syn.append({'AF': traj.data.loc[combination[0]]['AF'],
                                     'regularized_gradient': gradient[0]},
                                    ignore_index=True))

    # Step 2: Exclude all mutations with VAF < 0.01
    # lack of data below this threshold alters the distribution of gradients
    melt_syn = melt_syn[melt_syn['AF'] > 0.01]

    # Bin the dataframe according to VAF.
    melt_syn['bin'] = pd.cut(melt_syn['AF'], 25)

    # Step 3: Compute mean and variance in each bin using bootstrap.

    bin_centers = []     # track bin centers
    mean_dist_bins = []      # mean distribution in each bin using bootstrap
    var_dist_bins = []  # variance distribution in each bin using bootstrap
    bin_size = []       # Number of trajectories in each bin

    # Set seed for later bootstrap
    np.random.seed(1)

    for bin in melt_syn['bin'].unique():
        # filter dataset to obtain rows in bin
        VAF_bin = melt_syn[melt_syn['bin'] == bin]
        bin_size.append(len(VAF_bin))               # points in bin
        bin_centers.append(VAF_bin['AF'].mean())     # compute bin center

        # Bootstrap to compute gradient mean and variance distributions
        # This is used to plot confidence intervals
        n = len(VAF_bin)
        reps = 200   # number of drawn samples
        subsets = np.random.choice(VAF_bin['regularized_gradient'], (n, reps))
        mean_dist_bins.append(subsets.mean(axis=0))
        var_dist_bins.append(subsets.var(axis=0))

    # Compute mean and quantile of each bin's distribution
    # for mean and variance.
    mean, mean_quant = bin_statistics(mean_dist_bins)

    variance, var_quant = bin_statistics(var_dist_bins)

    # Step 4: Linear regressions.
    # Set weights for each bin used in weighted linear regression
    size_weight = np.array(bin_size)

    # Linear regression fit for Gradient ~ VAF
    mean_regr, mean_fig = linear_regression(
            data_type='mean', data=mean, bin_centers=bin_centers,
            weights=size_weight, quantiles=mean_quant)

    # Linear regression fit for Gradient variance ~ VAF
    var_regr, var_fig = linear_regression(
        data_type='variance', data=variance, bin_centers=bin_centers,
        weights=size_weight, quantiles=var_quant)

    # Stack plot showing variance and mean linear regressions and bin sizes
    stack_fig = plot.stack_plot(melt_syn, bin_centers,
                                mean,  mean_regr, mean_quant,
                                variance, var_regr, var_quant)

    # Check plot filter on synonymous trajectories
    check_fig = plot.check_plot(syn, mean_regr, var_regr,
                                mean_color, var_color, 2)

    return mean_regr, var_regr, [check_fig, mean_fig, var_fig, stack_fig]


def bin_statistics(dist_bins):
    """Given a array of a distribution of points by bins
    returns the mean and quantiles for each bin."""
    # compute the mean for each bin
    mean = np.mean(dist_bins, axis=1)

    # compute the distance from quantiles to the mean for each bin
    quantiles = mean - np.quantile(dist_bins, [0.1, 0.9], axis=1)
    quantiles = np.abs(quantiles)

    return mean, quantiles


def linear_regression(data_type, data, bin_centers, weights, quantiles):
    if data_type == 'mean':
        color = mean_color
        intercept = True
        yaxis = 'Gradient'
    elif data_type == 'variance':
        color = var_color
        intercept = False
        yaxis = 'Gradient\'s variance'
    else:
        return 'Wrong data type'

    # Mean weighted model
    regr = LinearRegression(fit_intercept=intercept)
    regr.fit(np.c_[bin_centers], data, weights)
    score = regr.score(np.c_[bin_centers], data, weights)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=bin_centers, y=data,
                   name='Data', mode='markers',
                   error_y=dict(
                       type='data',
                       symmetric=False,
                       array=quantiles[1],
                       arrayminus=quantiles[0])))

    # Trace linear fit
    x = [0, max(bin_centers)]
    fig.add_trace(
        go.Scatter(x=x, y=regr.predict(np.c_[x]),
                   marker_color=color,
                   name='linear regression fit'))

    fig.update_layout(title=(f'Linear regression {yaxis} ~ VAF <br>'
                             f'R2 score: {round(score,3)}'),
                      template='plotly_white',
                      legend=dict(
                          y=0.95,
                          x=1,
                      ))

    fig.update_xaxes(title='VAF')
    fig.update_yaxes(title=yaxis)

    return regr, fig


#endregion
# =============================================================================
#region : Load prior distributions
# =============================================================================



flat_fitness_prior = auxiliary.prior_load(
    '../Resources/flat_fitness_prior.npy')
exponential_fitness_prior = auxiliary.prior_load(
    '../Resources/exponential_fitness_prior.npy')
flat_N_w_prior = auxiliary.prior_load(
    '../Resources/flat_N_w_prior.npy')
binom_p_prior = auxiliary.prior_load(
    '../Resources/binom_p_prior.npy')
betabinom_p_prior =  auxiliary.prior_load(
    '../Resources/betabinom_p_prior.npy')
betabinom_beta_prior =  auxiliary.prior_load(
    '../Resources/betabinom_beta_prior.npy')

#endregion
# =============================================================================
#region : Artifact binomial model
# =============================================================================


def binom_artifact_fit(trajectories, params=None):
    """Compute the minimium negative log-likelihood of observing a time-series
    data, assuming the data was produced by gaussian noise.

    Parameters
    -----------
    trajectories: list of Pandas DataFrames. Each dataset must includes columns
                  'AO' and 'DP'.
    params: tuple (float, float). Tuple of floats determining the mean and
            standard deviation of the gaussian distribution used to initiate
            fitting proccess.
            By default this is set to mean: 0.01 and std_dev: 0.01.

    Returns
    -----------
    model: tuple (scipy.OptimizeResult, float). Optimized fit of a gaussian
           distribution to describe the time-series data. Bayesian information
           criterion.
     """

    if params is None:
        p = 0.01

    else:
        p = params

    # initial parameters
    x0 = np.array([p])
    bounds = [(0,0.5)]
    model = minimize(binom_artifact_nll, x0=x0, bounds=bounds,
                     args=trajectories, method='Nelder-Mead')

    return model


def binom_artifact_nll(params, trajectories, conditional=False):
    """Compute negative log-likelihood of observing a time-series
    data conditional on the previous time point, assuming the data was produced
    by gaussian noise with a given mean and variance.

    Parameters
    -----------
    params: tuple (float, float). Tuple of floats determining the mean and
            standard deviation of the gaussian distribution used to initiate
            fitting proccess.
    trajectories: list of Pandas DataFrames. Each dataset must includes columns
                  'AO' and 'DP'.

    returns
    -----------
    nll: float. negative log-likelihood of observing a given data.
     """

    if conditional is False:
        index = 0
    else:
        index = 1

    # Extract parameters
    p = params[0]
    ind_nll = []

    for ind_traj in trajectories:
        likelihood = binom.pmf(k=ind_traj[index:].AO,
                                  n=ind_traj[index:].DP,
                                  p=p)
        ind_nll.append(np.sum(np.log(likelihood)))

    return -np.sum(ind_nll)


def artifact_binom_init(trajectory, rep=100, mp=True):
    """Initiate artifact_fit function several times and extract optimal fit,
    to avoid local minima.

    Parameters
    -----------
    trajectory: Pandas DataFrame. Dataset including column 'AF' and 'age'
                with the allele frequency of a trajectory for each age.
    rep: float. Number of initiations.
    mp: Bool. Use multiprocessing. Default is False.


    returns
    -----------
    optimal_model: tuple (scipy.OptimizeResult, float). Optimized fit of a
                   gaussian distribution to describe the time-series data.
                   Bayesian information criterion.
     """
    # Set reproducibility seed
    np.random.seed(1)

    # Make sure that there are no DP values = 0 in dataset
    trajectory.loc[trajectory['DP'] == 0, 'DP'] = (
        trajectory[trajectory.DP != 0].DP.mean())
    trajectory['DP'] = trajectory['DP'].astype(int)

    # Select random parameter initiations
    p_init = np.linspace(artifact_size_bounds[0],
                         artifact_size_bounds[1],
                         rep+1)

    # Create partial fit function, fixing trajectory, for fit mapping
    artifact_binom_fit_partial = partial(artifact_binom_fit, trajectory)

    # Fit function using map with multiprocessing or not
    if mp is True:
        with Pool(8) as p:
            model_list = list(p.map(artifact_binom_fit_partial, p_init))
    else:
        model_list = list(map(artifact_binom_fit_partial, p_init))

    # Determine optimal model
    optimal_model = model_list[0]
    for model in model_list:
        if model.fit.fun < optimal_model.fit.fun:
            optimal_model = model

    optimal_model.data = trajectory
    return optimal_model


def binom_artifact_model_probability(data,
                                     prior_p=binom_p_prior,
                                     mutation_object=True):
    """Compute the marginalised likelihood of observing a time-series
    assuming that the data was produced by a binomial sequencing
    artifact.
    
    Parameters:
    -----------
    trajectories: list of Pandas DataFrames. Each dataset must include
                  columns 'AO' and 'DP' with the allele observations and
                  read depth of a genetic clone for each age.
    prior_p: array. 2D numpy array containing the prior distribution of
             binomial means."""

    if mutation_object is True:
        trajectories = data.data
    else:
        trajectories = data

    ind_likelihoods = []
    for ind_traj in trajectories:
        # initialise list of samples for p integration
        int_p = []
        # integral over betabinom p parameter
        for p_sample in prior_p[0, :]:
            # For each p compute the likelihood of observing
            # a given time-series conditional on
            # initial time-point.
            cond_p_likelihood = binom.pmf(k=ind_traj[1:].AO,
                                          n=ind_traj[1:].DP,
                                          p=p_sample)
            int_p.append(np.product(cond_p_likelihood))
        
        # compute the model probability for each individual
        marginal_ind_likelihood = np.trapz(x=prior_p[0,:], y=int_p*prior_p[1,:])
        ind_likelihoods.append(marginal_ind_likelihood)
    
    # Compute the product of likelihoods
    mutation_prob = np.product(ind_likelihoods)

    if mutation_object is True:
        # return updated model_comparison object 
        data.binom_artifact_prob = mutation_prob
        return data
    else:
        # return marginalised likelihood.
        return mutation_prob

#endregion
# =============================================================================
#region : Artifact beta-binomial model
# =============================================================================

def betabinom_artifact_nll(params, trajectories, conditional=False):

    # For each combination of p and beta parameters compute
    # lieklihood of observing a given time-series conditional
    # on initial time-point.
    # The likelihood for each individual omits first time point

    if conditional is False:
        index = 0
    else:
        index = 1
    
    # Extract beta binomial parameters
    p = params[0] 
    b = params[1]
    a = b*p/(1-p)

    ind_nll = []
    for ind_traj in trajectories:
        likelihood = betabinom.pmf(k=ind_traj[index:].AO,
                                   n=ind_traj[index:].DP,
                                   a=a,
                                   b=b)
        ind_nll.append(np.sum(np.log(likelihood)))

    return -np.sum(ind_nll)


def betabinom_artifact_fit(trajectory):
    p = 0.01
    b = 10e5
    bounds = ((0, 0.05), (0,10e9))

    # initial parameters
    x0 = np.array([p, b])

    model = minimize(betabinom_artifact_nll, x0=x0,
                     bounds=bounds,
                     args=trajectory,
                     method='Nelder-Mead')

    return model


def betabinom_artifact_model_probability(data,
                                         prior_p=betabinom_p_prior[0],
                                         prior_beta=betabinom_beta_prior[0],
                                         mutation_object=True):
    """Compute the marginalised likelihood of observing a time-series
    assuming that the data was produced by a beta binomial sequencing
    artifact.
    
    Parameters:
    -----------
    trajectories: list of Pandas DataFrames. Each dataset must include
                  columns 'AO' and 'DP' with the allele observations and read depth
                  of a genetic clone for each age.
    prior_p: array. 2D numpy array containing the prior distribution of
             beta binomial means.
    prior_beta: array. 2D numpy array containing the prior distribution of
                beta binomial beta parameter."""

    if mutation_object is True:
        trajectories = data.data
    else:
        trajectories = data

    # initialise list of samples for p integration
    int_p = []
    # integral over betabinom p parameter
    for p_sample in prior_p[0, :]:
        # initialise list of samples for beta integration
        int_beta = []
        # integrate over beta parameter
        for beta_sample in prior_beta[0, :]:
            # For each combination of p and beta parameters compute
            # lieklihood of observing a given time-series conditional
            # on initial time-point.
            alpha_sample = beta_sample*p_sample/(1-p_sample)
            # compute likelihood for each individual omiting first time point
            ind_likelihoods = []
            for ind_traj in trajectories:
                likelihood = betabinom.pmf(k=ind_traj[1:].AO,
                                           n=ind_traj[1:].DP,
                                           a=alpha_sample,
                                           b=beta_sample)
                ind_likelihoods.append(np.product(likelihood))
            # for each beta append the total likelihood (product individuals).
            int_beta.append(np.product(ind_likelihoods))
        # For each p, compute the likelihood marginalised over beta
        int_p.append(
            np.trapz(x=prior_beta[0, :],
                     y=int_beta*prior_beta[1, :]))
    
    # marginalise likelihood over p
    mutation_prob = np.trapz(x=prior_p[0, :],
                             y=int_p*prior_p[1, :])

    if mutation_object is True:
        # return updated model_comparison object 
        data.betabinom_artifact_prob = mutation_prob
        return data
    else:
        # return marginalised likelihood.
        return mutation_prob



#endregion
# =============================================================================
#region: BD model
# =============================================================================

N_w_bounds = (5_000, 200_000)
fitness_bounds = (0, 0.5)

def bd_nll(params, trajectory, return_params=False):
    """Compute negative log-likelihood of observing a time-series
    data conditional on the previous time point, assuming the data was produced
    by a birth and death process.

    Parameters
    -----------
    params: tuple (float, float). Tuple of floats determining the fitness
            (s_fit), total number (N_w) of stem cells and time of mutation
            acquisition (t_fit) used compute the nll.
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    return_params: Bool. Returns params. This is used during asynchronous
                   multiprocessing fitting, where parameters are not used
                   in order.
    returns
    -----------
    nll: float. negative log-likelihood of observing a given data.
    """
    # Extract parameters
    s_fit = params[0]
    t_fit = params[2]
    N_w_fit = int(params[1])

    # Exatract inferred clone_sizes from AO/DP ratio
    mean_size, size_range = observations_to_clone_size(AO=trajectory.AO,
                                                       DP=trajectory.DP,
                                                       N_w=N_w_fit)

    # Set inferred size as the mean of all possible ranges of observations
    trajectory['inferred_size'] = mean_size

    # Compute time_steps
    trajectory['delta_t'] = np.insert(np.diff(trajectory.age),
                                      0, trajectory.iloc[0].age - t_fit)
    # Initialise negative log-likelihood computation
    nll = 0
    for i, time_point in trajectory.iterrows():
        # Extract initial clone_size and time difference between observations
        if i == 0:
            init_size = 1
        else:
            init_size = max(trajectory.iloc[i-1].inferred_size, 1)

        # Compute AO/DP observation probability
        prob = AO_prob_value(AO=time_point.AO,
                             DP=time_point.DP,
                             init_size=init_size,
                             s=s_fit,
                             delta_t=time_point.delta_t,
                             N_w=N_w_fit)

        # Avoid divide by zero encountered in log warning
        if prob < 1.0e-100:
            prob = 1.0e-100

        # Compute negative log likelihood
        nll -= np.log(prob)

    if return_params is True:
        return nll, params
    else:
        return nll


def bd_csmargin_nll(params, trajectory,
                    return_params=False,
                    AO_random_samples=50):
    """Compute negative log-likelihood of observing a time-series
    data conditional on the previous time point, assuming the data was produced
    by a birth and death process.

    Parameters
    -----------
    params: tuple (float, float). Tuple of floats determining the fitness
            (s_fit), total number (N_w) of stem cells and time of mutation
            acquisition (t_fit) used compute the nll.
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    cs_subset_init: float. Number of initial clone sizes used to sample the
                    probability marginalisation over possible initial sizes.
    return_params: Bool. Returns params. This is used during asynchronous
                   multiprocessing fitting, where parameters are not used
                   in order.
    returns
    -----------
    nll: float. negative log-likelihood of observing a given data.
    """

    # Extract parameters
    s_fit = params[0]
    N_w_fit = int(params[1])
    t_fit = params[2]

    # Compute time_steps
    trajectory['delta_t'] = np.insert(np.diff(trajectory.age),
                                      0, trajectory.iloc[0].age - t_fit)

    nll = 0
    for i, time_point in trajectory.iterrows():
        if i == 0:
            init_clone_size_range = [1]
            init_clone_size_prob = [1]
        else:
            init_time_point = trajectory.iloc[i-1]
            # Sample AO values given p=AO/DP
            AO_sample = binom.rvs(n=int(init_time_point.DP),
                                  p=init_time_point.AO/init_time_point.DP,
                                  size=AO_random_samples)
            AO_sample = np.unique(AO_sample)

            init_clone_size_range = distributions.vaf_to_clone_size(
                                     AO_sample/init_time_point.DP,
                                     N_w_fit).astype(int)

            p = init_clone_size_range / (2*N_w_fit + 2*init_clone_size_range)
            init_clone_size_prob = binom.pmf(init_time_point.AO.astype('int'),
                                             n=init_time_point.DP,
                                             p=p)

        # Compute probability for each initial clone_size and marginalise
        # over all possib
        AO_prob = 0
        for init_size, init_size_prob in zip(init_clone_size_range,
                                             init_clone_size_prob):

            # Compute AO probability conditional on init_size
            prob = AO_prob_value(AO=time_point.AO,
                                 DP=time_point.DP,
                                 init_size=init_size,
                                 s=s_fit,
                                 delta_t=time_point.delta_t,
                                 N_w=N_w_fit)

            # Marginalise over the probability of init_size
            AO_prob += prob*init_size_prob

        # avoid divide by zero encountered in log warning
        if AO_prob < 1.0e-10:
            AO_prob = 1.0e-10

        nll -= np.log(AO_prob)
    if return_params is True:
        return nll, params
    else:
        return nll


def bd_fit(params=None, trajectory=None, method='Nelder-Mead', fit_func=bd_nll):
    """Compute the minimium negative log-likelihood of observing a time-series
    data, assuming the data was produced by a birth and death process.

    Parameters
    -----------
    params: tuple (float, float). Tuple of floats determining the fitness
            and total number (N_w) of stem cells in used to iniate the fitting
            proccess. By default this is set to their bounds' mean.
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    method: str. Solver algorithm used during minimization process.
            A full list of allowed strings can be found in the documentation
            of scipy.optimize.minimize:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

    returns
    -----------
    model: trajectory_model class. Optimized fit of a birth and
           death distribution to describe the time-series data.
     """

    if trajectory is None:
        print('No trajectory error')
        return
    # parameter bounds
    t0_bounds = (0, trajectory.iloc[0].age-1)
    if params is None:
        s_init = np.mean(fitness_bounds)
        N_w_init = np.mean(N_w_bounds)
        t0_init = np.mean(t0_bounds)

    else:
        s_init, N_w_init, t0_init = params

    # Make sure that there are no DP values = 0 in dataset
    trajectory.loc[trajectory['DP'] == 0, 'DP'] = (
        trajectory[trajectory.DP != 0].DP.mean())
    trajectory['DP'] = trajectory['DP'].astype(int)

    bounds = (fitness_bounds, N_w_bounds, t0_bounds)

    # initial parameters
    x0 = np.array([s_init, N_w_init, t0_init])

    model = minimize(fit_func, x0=x0, bounds=bounds,
                     args=trajectory, method=method)

    model_class = trajectory_model(model_type='bd_process',
                                   fit=model,
                                   data=trajectory)

    return model_class


def bd_init(trajectory, rep=10, method='Nelder-Mead', mp=True,
            fix_N=None, brute=True, fit_func=bd_nll):
    """Initiate artifact_fit function several times and extract optimal fit,
    to avoid local minima.

    Parameters
    -----------
    trajectory: Pandas DataFrame. Dataset including column 'AF' and 'age'
                with the allele frequency of a trajectory for each age.
    rep: float. Number of initiations.
    mp: Bool. Use multiprocessing. Default is False.
    fit: Bool. Choose to evaluate nll in given parameters or to
         optimize parameters.


    returns
    -----------
    optimal_model: tuple (scipy.OptimizeResult, float). Optimized fit of a
                   birth and death process to describe the time-series data.
                   Bayesian information criterion.
     """

    # Make sure that there are no DP values = 0 in dataset
    trajectory.loc[trajectory['DP'] == 0, 'DP'] = (
        trajectory[trajectory.DP != 0].DP.mean())
    trajectory['DP'] = trajectory['DP'].astype(int)

    # Select random parameter initiations
    fitness_range = np.linspace(fitness_bounds[0], fitness_bounds[1], rep+1)
    t0_range = np.linspace(0, trajectory.iloc[0].age-1, rep+1).astype('int')

    # the number of wild type stem cells can be fix using fix_N parameter
    if fix_N is None:
        # if fix_N is not set, then compute a range of valid stem cell counts
        N_w_range = np.linspace(N_w_bounds[0], N_w_bounds[1], rep+1, dtype=int)
    else:
        N_w_range = [fix_N]

    # Create all possible combinations of initial parameters
    params_init = list(product(fitness_range, N_w_range, t0_range))

    # Set fitting function
    if brute is False:
        partial_func = partial(bd_fit,
                               trajectory=trajectory,
                               method=method,
                               fit_func=fit_func)
    else:
        # Set fit_func function to return parameters
        fit_func_return_params = partial(fit_func,
                                         return_params=True)
        partial_func = partial(fit_func_return_params,
                               trajectory=trajectory)

    if mp is True:
        with Pool(8) as p:
            model_list = list(p.map(partial_func, params_init))
    else:
        model_list = list(map(partial_func, params_init))

    if brute is False:
        # Optimal model
        optimal_model = model_list[0]
        for model in model_list:
            if model.fit.fun < optimal_model.fit.fun:
                optimal_model = model

    else:
        # Unpack nll and parameter values from model_list
        nll, params = zip(*model_list)
        s, N_w, t0 = zip(*params)

        # Convert negative log-likelihoods to likelihood
        likelihood = np.exp(-np.array(nll))
        # Create DataFrame with the nll for each combination of parameters
        brute_force_df = pd.DataFrame({'fitness': s,
                                       'N_w': N_w,
                                       't0': t0,
                                       'likelihood': likelihood})
        # find dataframe row of optimal nll
        optimal_idx = brute_force_df['likelihood'].idxmax()

        # Fit new model with optimal parameters as initial parameters
        model_fit = bd_fit(params=[brute_force_df.iloc[optimal_idx].fitness,
                                   brute_force_df.iloc[optimal_idx].N_w,
                                   brute_force_df.iloc[optimal_idx].t0],
                           trajectory=trajectory,
                           method='Nelder-Mead',
                           fit_func=fit_func)

        # Create model class object from optimal model
        optimal_model = trajectory_model(model_type='bd_process',
                                         fit=model_fit.fit,
                                         data=trajectory)

        # Create heatmaps for each combination of 2 parameters
        heatmaps = likelihood_heatmaps(brute_force_df)
        conditional_distribution_plots = conditional_distributions(
                                            brute_force_df)

        # Add dataframe and heatmap as attribute of the class object.
        optimal_model.brute_force_df = brute_force_df
        optimal_model.heatmap = heatmaps
        optimal_model.conditional_dist = conditional_distribution_plots

    return optimal_model


def AO_prob_value(AO, DP, init_size, s, delta_t, N_w, size_sample=1_000):
    """Probability of observing a certain allele observation from time-series
    data, assuming a birth and death process.
    """

    # If init_size = 0 only AO=0 can be observed
    if (init_size == 0) & (AO == 0):
        return 1
    if (init_size == 0) & (AO != 0):
        return 0

    # extract mean and variance of BD process
    mean, variance = distributions.BD_stats(init_size=init_size,
                                            s=s,
                                            delta_t=delta_t)

    # negative binomial probability and number of successes
    p = mean / variance
    n = np.power(mean, 2) / (variance-mean)
    # Extract a 1_000 clone_sizes according to their probability
    random_size_states = nbinom.rvs(n, p,
                                    size=size_sample,
                                    random_state=123)

    unique_random_size_states = np.unique(random_size_states)

    # Compute the associated probability for each clone size
    # of being observed in a BD process
    clone_size_bd_prob = distributions.nb_approx_pmf(
                    size=unique_random_size_states,
                    init_size=init_size,
                    s=s,
                    delta_t=delta_t)

    # Compute the allele proportion for each resulting clone_size
    total_alleles = 2*(N_w + unique_random_size_states)
    allele_proportion_range = unique_random_size_states / total_alleles

    # Compute the probability of observing AO given DP and the allele proportion
    # after sequencing using binomial distribution
    sequencing_prob = binom.pmf(AO, DP, p=allele_proportion_range)

    # marginalise the sequencing probability by clone size
    sequencing_prob = np.trapz(y=clone_size_bd_prob*sequencing_prob,
                               x=unique_random_size_states)

    return sequencing_prob


def observations_to_clone_size(AO, DP, N_w):
    """Find all clone sizes leading to the same observed AO/DP ratio.
    We assume that the binomial draw selects exactly
    AO = DPx/(2N+2x), where x is the true clone size.
    We then use the inverse function of
    f(x) = DPx/(2N+2x), f-1(x) = 2Nx/(DP-2x)
    to estimate upper and lower bounds of all possible clone sizes.

    Parameters:
    -----------
    AO: int. Alternate allele depth.
    DP: int. Total deepth at position.
    N_w: int. Wild-type stem cells count.

    Returns:
    -----------
    mean_size, size_range. Tuple. Estimated average clone size and range
                           of inferred clone sizes.
    """

    # Compute extreme approximate sizes
    min_size = 2*N_w*AO/(DP-2*AO)
    max_size = 2*N_w*(AO+1)/(DP-2*(AO+1))

    # Deal with the case AO = 0, min_size = 0
    min_size = np.where(min_size == 0, 1, min_size)

    # Retrieve clone sizes in the extremes.
    size_range = np.array([np.ceil(min_size), np.floor(max_size)]).astype(int)
    # Note: unless N_w = 0 , min_size_range = 1.

    # Compute average size
    mean_size = size_range.mean(axis=0).astype(int)

    return mean_size, size_range


def AO_DP_to_init_clone_size_range(DP, p, N_w, std=4):
    min_AO = np.array(DP*p - std*np.sqrt(DP*p*(1-p)), dtype=int)
    max_AO = np.array(DP*p + std*np.sqrt(DP*p*(1-p)), dtype=int)

    min_AO = np.where(min_AO < 0, 0, min_AO)
    max_AO = np.where(max_AO > DP/2-1, DP/2-1, max_AO)

    min_clone_size = distributions.vaf_to_clone_size(min_AO/DP, N_w)
    max_clone_size = distributions.vaf_to_clone_size(max_AO/DP, N_w)

    clone_size_range = np.array([min_clone_size, max_clone_size])
    clone_size_range = clone_size_range.astype('int')

    return clone_size_range


def inverse_min_max(param, param_bounds):
    return param * (param_bounds[1]-param_bounds[0]) + param_bounds[0]


def heatmap(df, combination, condition):
    columns = list(combination)
    columns.append('likelihood')

    df_heatmap = df[columns]
    # Create a heatmap of nll as a function of parameters
    pivot_df = df_heatmap.pivot(columns[0], columns[1], 'likelihood')
    pivot_df = pivot_df.reindex(index=pivot_df.index[::-1])
    heatmap_fig = px.imshow(pivot_df,
                            title=(f'Posterior distribution of'
                                   f'{columns[0]} and {columns[1]}<br>'
                                   f'conditional on '
                                   f'{condition[0]} = '
                                   f'{round(condition[1],2)}'),
                            labels=dict(x=columns[1], y=columns[0],
                                        color="Likelihood"),
                            x=np.round(pivot_df.columns, 2).astype('str'),
                            y=np.round(pivot_df.index, 2).astype('str'))
    return heatmap_fig


def likelihood_heatmaps(df):
    heatmaps = []
    # For each combination, the remaining parameter is set to the optimal
    # fit.

    # Find optimal_fit
    optimal_idx = df['likelihood'].idxmax()

    parameter_set = set(df.columns[:-1])
    for comb in combinations(parameter_set, 2):
        # find conditional parameter from a combination of parameters
        comb = set(comb)
        conditional_param = list(parameter_set - comb)[0]
        optimal_conditional_param = df.iloc[optimal_idx][conditional_param]

        # Store condition
        condition = (conditional_param, optimal_conditional_param)

        # Extract rows of dataframe with optimal dropped parameter
        df_heat = df[df[conditional_param]
                     == optimal_conditional_param]
        # Create heatmap plot of nll as a function of comb parameters
        heatmap_fig = heatmap(df_heat, comb, condition)
        heatmaps.append(heatmap_fig)

    return heatmaps


def conditional_distributions(df):
    """Plot distributions for each parameter, conditional on all
    other parameters set to optimal.

    Parameters
    -----------
    df: pandas DataFrame. DataFrame as returned by bd_init when
               brute_force=True

    returns:
    -----------
    conditional_distributions: list of plotly Figures. List containing a plot
                               for each possible conditional distribution.
    -----------
    """

    conditional_distributions = []
    # For each combination, the remaining parameter is set to the optimal
    # fit.

    # Find optimal_fit
    optimal_idx = df['likelihood'].idxmax()

    parameter_set = set(df.columns) - {'likelihood'}
    for param in parameter_set:
        # find dropped parameter from a combination
        conditional_params = list(parameter_set - {param})
        optimal_cond_params = (
            df.iloc[optimal_idx][conditional_params].values)

        conditional_df = df.copy()
        for cond_param, optimal in zip(conditional_params, optimal_cond_params):
            conditional_df = conditional_df[conditional_df[cond_param] == optimal]

        conditional_distributions.append(
            px.line(conditional_df,
                    x=param, y='likelihood',
                    title=(f'{param} distribution <br>'
                           f'conditional on {conditional_params}'
                           f' = {optimal_cond_params}')))

    return conditional_distributions


#endregion
# =============================================================================
#region: Full BD model
# =============================================================================

# Set the minimum value of a log used to avoid divide by zero
min_nll = np.log(np.nextafter(0, 1))

fitness_bounds = (0, 0.5)

def simple_conditional_N_nll(params, trajectory, N_w=50_000,
                             return_params=False):
    """Compute negative log-likelihood of observing a time-series
    data conditional on the previous time point, assuming the data was produced
    by a birth and death process and conditional on a fixed number of stem
    cells N_w.

    This function assumes a unique initial clone size at each time-step.

    Parameters
    -----------
    params: tuple (float, float). Tuple of floats determining the fitness
            (s_fit) and time of mutation acquisition (t_fit) used compute
            the nll.
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    N_w: int. Fixed number of wild type stem cells.
    return_params: Bool. Returns params. This is used during asynchronous
                   multiprocessing fitting, where parameters are not used
                   in order.
    returns
    -----------
    nll: float. negative log-likelihood of observing a given data.
    """
    # Extract parameters
    s_fit = params[0]
    # Exatract inferred clone_sizes from AO/DP ratio
    mean_size, size_range = observations_to_clone_size(AO=trajectory.AO,
                                                       DP=trajectory.DP,
                                                       N_w=N_w)

    # Set inferred size as the mean of all possible ranges of observations
    trajectory['inferred_size'] = mean_size

    # Compute time_steps
    trajectory['delta_t'] = np.insert(np.diff(trajectory.age),
                                      0, 0)
    # Initialise negative log-likelihood computation
    nll = 0
    for i, time_point in trajectory.iterrows():
        # Extract initial clone_size and time difference between observations
        if i == 0:
            # init_size = 1
            continue
        else:
            init_size = max(trajectory.iloc[i-1].inferred_size, 1)

        # Compute AO/DP observation probability
        prob = AO_prob_value(AO=time_point.AO,
                             DP=time_point.DP,
                             init_size=init_size,
                             s=s_fit,
                             delta_t=time_point.delta_t,
                             N_w=N_w)

        # Compute negative log likelihood
        # Avoiding divide by zero
        if prob != 0:
            nll -= np.log(prob)

        else:
            # use min_nll = np.log(nextafter(0,1))
            nll -= min_nll

    if return_params is True:
        return nll, params
    else:
        return nll


def clone_size_prior(trajectory, N_w, prior_resolution=10, alpha=0.01):
    """Computes the probability mass functiondistributions for inital
    clone_sizes"""
    AO, DP = trajectory.AO, trajectory.DP
    # Find 95% confidence intervals for binomial p using a beta function
    p_conf_int = beta.ppf(q=[alpha, 1-alpha],
                          a=AO+1, b=DP-AO+1)

    # List of binomial p ranges and clone_size_ranges
    p_range = np.linspace(p_conf_int[0], p_conf_int[1], prior_resolution)

    # Extract binomial proportion probabilities
    p_prob = beta.pdf(x=p_range, a=AO+1, b=DP-AO+1)

    init_clone_size_range = distributions.vaf_to_clone_size(p_range, N_w)
    init_clone_size_range = init_clone_size_range.astype('int')

    # Create prior
    prior = [init_clone_size_range, p_prob]

    return prior


def binomial_sequencing_probability(AO, DP, init_clone_size, s, delta_t, N_w,
                                    resulting_clone_size_resolution=10):
    if init_clone_size == 0:
        # if init_clone_size = 0 the only possible outcom is AO=0 
        if AO == 0:
            return 1
        else:
            return 0


    # extract mean and variance of BD process
    mean, variance = distributions.BD_stats(init_size=init_clone_size,
                                            s=s,
                                            delta_t=delta_t)

    # negative binomial probability and number of successes
    p = mean / variance
    n = np.power(mean, 2) / (variance-mean)

    # extract 99% confidence interval of resulting clone sizes
    resulting_clone_size_interval = nbinom.interval(0.99, n=n, p=p)
    # Extract range of clone_size
    resulting_clone_size_range = np.linspace(resulting_clone_size_interval[0],
                                             resulting_clone_size_interval[1],
                                             resulting_clone_size_resolution, dtype='int')

    # Compute resulting clone size probabilities
    resulting_clone_size_prob = nbinom.pmf(
        k=resulting_clone_size_range, n=n, p=p)

    # Compute the allele proportion for each resulting clone_size
    total_alleles = 2*(N_w + resulting_clone_size_range)
    allele_proportion_range = resulting_clone_size_range / total_alleles

    # Compute the probability of observing AO given DP and the allele proportion
    # in a sequencing experiment using the binomial distribution
    sequencing_prob = binom.pmf(AO, DP,
                                p=allele_proportion_range)

    # marginalise the sequencing probability by clone size
    AO_prob_cond_init = np.trapz(y=resulting_clone_size_prob*sequencing_prob,
                                 x=resulting_clone_size_range)

    return AO_prob_cond_init


def bd_conditional_N_nll(params, N_w, trajectory, return_params=False,
                         init_clone_size_resolution=10):
    """Compute negative log-likelihood of observing a time-series
    data conditional on the previous time point, assuming the data was produced
    by a birth and death process.

    Same as simple_conditional_N_nll function, but marginalising over the
    distribution of initial clone sizes.

    Parameters
    -----------
    params: tuple (float, float). Tuple of floats determining the fitness
            (s_fit), total number (N_w) of stem cells and time of mutation
            acquisition (t_fit) used compute the nll.
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    return_params: Bool. Returns params. This is used during asynchronous
                   multiprocessing fitting, where parameters are not used
                   in order.
    returns
    -----------
    nll: float. negative log-likelihood of observing a given data.
    """
    s = params[0]
    
    # Compute time_steps
    trajectory['delta_t'] = np.insert(np.diff(trajectory.age),
                                      0, 0)

    init_clone_size_prior = clone_size_prior(trajectory, N_w=N_w,
                                             prior_resolution=init_clone_size_resolution)

    # Initialise nll
    nll = 0
    #select time_point
    for i, row in trajectory[1:].iterrows():
        # Create list of AO probabilities conditional on initial clone size
        cond_init_prob = []
        # fix initial_clone_size
        for j in range(init_clone_size_resolution):
            cond_init_prob.append(
                binomial_sequencing_probability(AO=row.AO,
                                                DP=row.DP,
                                                init_clone_size=(
                                                    init_clone_size_prior[0][j, i-1]),
                                                s=s,
                                                delta_t=row.delta_t,
                                                N_w=N_w))

        # marginalise AO_prob_cond_init with respect to initial clone_size
        AO_prob = np.trapz(y=cond_init_prob*init_clone_size_prior[1][:, i-1],
                           x=init_clone_size_prior[0][:, i-1])

        # Compute negative log likelihood
        # Avoiding divide by zero
        if AO_prob != 0:
            nll -= np.log(AO_prob)

        else:
            # use min_nll = np.log(nextafter(0,1))
            nll -= min_nll

    if return_params is False:
        return nll
    else:
        return nll, params


def bd_marginal_N_nll(params, trajectory,
                      nll_conditional_N_function=bd_conditional_N_nll,
                      N_w_prior=flat_N_w_prior, return_params=False):

    cond_likelihood = []
    for N_w in N_w_prior[0, :]:
        cond_nll = nll_conditional_N_function(params=params,
                                              trajectory=trajectory,
                                              N_w=N_w)
        cond_likelihood.append(np.exp(-cond_nll))

    marginal_likelihood = np.trapz(
        y=N_w_prior[1, :]*cond_likelihood, x=N_w_prior[0, :])

    # Compute negative log likelihood
    # Avoiding divide by zero
    if marginal_likelihood != 0:
        marginal_nll = -np.log(marginal_likelihood)
    else:
        # use min_nll = np.log(nextafter(0,1))
        marginal_nll = -min_nll

    if return_params is True:
        return marginal_nll, params
    else:
        return marginal_nll


def bd_marginal_likelihood(trajectory,
                           nll_conditional_N_function=bd_conditional_N_nll,
                           N_w_prior=flat_N_w_prior,
                           fitness_prior=flat_fitness_prior):

    cond_s_likelihood = []
    for fitness in fitness_prior[0, :]:
        cond_s_N_likelihood = []
        for N_w in N_w_prior[0, :]:
            # for each N_w append the conditional likelihood
            cond_nll = nll_conditional_N_function(params=[fitness],
                                                trajectory=trajectory,
                                                N_w=N_w)
            cond_s_N_likelihood.append(np.exp(-cond_nll))

        # Integrate conditional likelihoods to compute the marginal
        cond_s_likelihood.append(
            np.trapz(y=cond_s_N_likelihood*N_w_prior[1, :],
                     x=N_w_prior[0, :]))
        
    marginal_likelihood = np.trapz(y=cond_s_likelihood*fitness_prior[1, :],
                                   x=fitness_prior[0, :])
    
    # Create trajectory class object
    model = trajectory_model(model_type='bd_process',
                                   data=trajectory,
                                   likelihood=marginal_likelihood,
                                   nll=-np.log(marginal_likelihood),
                                   n_params=0)

    # Add prior and posteerior attributes
    model.prior = [fitness_prior, N_w_prior]
    return model
   

def bd_fitness_fit(params=None, trajectory=None,
                   method='Nelder-Mead',
                   minimize_func=bd_marginal_N_nll,
                   nll_conditional_N_function=bd_conditional_N_nll,
                   ):
    """Compute the minimium negative log-likelihood of observing a time-series
    data, assuming the data was produced by a birth and death process and a
    fixed number of stem cells.

    Parameters
    -----------
    params: tuple (float, float). Tuple of floats determining the fitness
            and time of mutation initiaion used to iniate the fitting
            proccess. By default this is set to their bounds' mean.
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    N_w: int. Number of wild type stem cells.
    method: str. Solver algorithm used during minimization process.
            A full list of allowed strings can be found in the documentation
            of scipy.optimize.minimize:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

    returns
    -----------
    model: trajectory_model class. Optimized fit of a birth and
           death distribution to describe the time-series data.
     """

    if trajectory is None:
        print('No trajectory error')
        return

    # parameter bounds
    if params is None:
        s_init = np.mean(fitness_bounds)

    else:
        s_init = params

    # Make sure that there are no DP values = 0 in dataset
    trajectory.loc[trajectory['DP'] == 0, 'DP'] = (
        trajectory[trajectory.DP != 0].DP.mean())
    trajectory['DP'] = trajectory['DP'].astype(int)

    bounds = (fitness_bounds,)

    # initial parameters
    x0 = [s_init]

    model = minimize(minimize_func, x0=x0, bounds=bounds,
                     args=(trajectory, nll_conditional_N_function), method=method)

    model_class = trajectory_model(model_type='bd_process',
                                   fit=model,
                                   data=trajectory,
                                   nll=model.fun,
                                   n_params=1)

    return model_class


def bd_fitness_init(trajectory, rep=10, method='Nelder-Mead', mp=True):
    """Initiate artifact_fit function several times and extract optimal fit,
    to avoid local minima.

    Parameters
    -----------
    trajectory: Pandas DataFrame. Dataset including column 'AF' and 'age'
                with the allele frequency of a trajectory for each age.
    rep: float. Number of initiations.
    mp: Bool. Use multiprocessing. Default is False.
    fit: Bool. Choose to evaluate nll in given parameters or to
         optimize parameters.


    returns
    -----------
    optimal_model: tuple (scipy.OptimizeResult, float). Optimized fit of a
                   birth and death process to describe the time-series data.
                   Bayesian information criterion.
     """

    # Make sure that there are no DP values = 0 in dataset
    trajectory.loc[trajectory['DP'] == 0, 'DP'] = (
        trajectory[trajectory.DP != 0].DP.mean())
    trajectory['DP'] = trajectory['DP'].astype(int)

    # Select random parameter initiations
    fitness_range = np.linspace(fitness_bounds[0], fitness_bounds[1], rep)
    fitness_range = [[s] for s in fitness_range]

    # Set fit_func function to return parameters
    partial_func = partial(bd_marginal_N_nll,
                           trajectory=trajectory,
                           return_params=True)

    if mp is True:
        with Pool(8) as p:
            model_list = list(p.map(partial_func, fitness_range))
    else:
        model_list = list(map(partial_func, fitness_range))

    # Unpack nll and parameter values from model_list
    nll, fitness = zip(*model_list)
    fitness = [s[0] for s in fitness]

    # Convert negative log-likelihoods to likelihood
    likelihood = np.exp(-np.array(nll))
    # Create DataFrame with the nll for each combination of parameters
    brute_force_df = pd.DataFrame({'fitness': fitness,
                                   'likelihood': likelihood})
    # find dataframe row of optimal nll
    optimal_idx = brute_force_df['likelihood'].idxmax()

    # Fit new model with optimal parameters as initial parameters
    optimal_model = bd_fitness_fit(params=brute_force_df.iloc[optimal_idx].fitness,
                                   trajectory=trajectory,
                                   method=method)

    # Add dataframe and heatmap as attribute of the class object.
    optimal_model.brute_force_df = brute_force_df
    return optimal_model


def likelihood_heatmap_N(df):

    columns = df.columns
    pivot_df = df.pivot(columns[0], columns[1], 'likelihood')

    # Create a heatmap of nll as a function of parameters
    pivot_df = pivot_df.reindex(index=pivot_df.index[::-1])
    heatmap_fig = px.imshow(pivot_df,
                            title=(f'Posterior distribution of'
                                   f'{columns[0]} and {columns[1]}'),
                            labels=dict(x=columns[1], y=columns[0],
                                        color="Likelihood"),
                            x=np.round(pivot_df.columns, 2).astype('str'),
                            y=np.round(pivot_df.index, 2).astype('str'))

    return heatmap_fig


#endregion
# =============================================================================
#region: Final BD model
# =============================================================================

def bd_process_conditional_likelihood_s_N(trajectory, s=0, N_w=50_000, alpha=0.01, resolution=25):
    """Compute the likelihood of observing a time-series given by the birth
    and death model, conditional on a given fitness, s, and a given number
    wild type stem cells, N_W.
    
    Parameters:
    -----------
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    s: float. Fitness.
    N_w: float. Number of wild-type stem cells.
    alpha: float. quantile used for distribution approximation.
    resolution: float. Number of points used to approximate integral over
                prior distribution.
    
    Returns:
    -----------
    float. Time series likelihood conditional on s and N_w."""
    time_point_likelihood = []
    for i in range(len(trajectory[:-1])):

        init_time_point = trajectory.iloc[i]
        next_time_point = trajectory.iloc[i+1]

        # Beta function
        beta_p_conf_int = beta.ppf(q=[alpha, 1-alpha],
                            a=init_time_point.AO+1,
                            b=(init_time_point.DP
                                - init_time_point.AO
                                + 1))

        # List of binomial p ranges and clone_size_ranges
        beta_p_range = np.linspace(beta_p_conf_int[0],
                            beta_p_conf_int[1],
                            resolution)

        # We only want to integrate over beta_p_range < 0.5
        beta_p_range = beta_p_range[beta_p_range < 0.5]

        # Compute the pdf associated to each p value
        beta_p_prob = beta.pdf(x=beta_p_range,
                            a=init_time_point.AO+1,
                            b=(init_time_point.DP
                                - init_time_point.AO
                                + 1))

        # Transform p values to clone sizes
        init_cs_range = distributions.vaf_to_clone_size(v=beta_p_range,
                                                        N_w=N_w)

        total_AO_prob = []
        for init_cs in init_cs_range:
            # Special case of initial clone size 0
            if init_cs == 0:
                if next_time_point.AO == 0:
                    total_AO_prob.append(1)
                else:
                    total_AO_prob.append(0)
                
                # continue to next initital clone size 
                continue
            
            # Compute next time point clone size range

            # Extract proportion and trials NB approximating BD process
            nb_p, nb_n = distributions.BD_parametrization(
            init_size=init_cs, s=s, delta_t=init_time_point.delta_t)

            # Compute next clone size range interval
            next_cs_int = nbinom.ppf(q=[alpha, 1-alpha],
                                    n=nb_n, p=nb_p)
            # Create next time point clone size range
            next_cs_range = np.linspace(next_cs_int[0], next_cs_int[1], resolution, dtype='int')
            
            # Compute associated BD probabilities
            next_cs_prob = nbinom.pmf(k=next_cs_range, n=nb_n, p=nb_p)

            # Transform clone sizes to VAFs
            next_binom_p_range = distributions.clone_size_to_vaf(next_cs_range, N_w=N_w)

            next_cs_AO_prob = binom.pmf(k=next_time_point.AO,
                                        n=next_time_point.DP,
                                        p=next_binom_p_range)

            total_AO_prob.append(
                np.trapz(x=next_cs_range, y=next_cs_prob*next_cs_AO_prob))

        time_point_likelihood.append(
            np.trapz(x=beta_p_range, y=total_AO_prob*beta_p_prob))

    time_series_likelihood = np.product(time_point_likelihood)

    return time_series_likelihood


def bd_process_model_probability(data,
                                 fitness_prior=flat_fitness_prior,
                                 N_w_prior=flat_N_w_prior,
                                 mutation_object=True):
    """Marginalise the conditional likelihood of observing a
    given time-series data.

    Parameters:
    -----------
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    fitness_prior: array. 2-D numpy array containing x and y values of
                   a sampled fitness prior.
    N_w_prior: array. 2-D numpy array containing x and y values of
               a sampled wild-type stem cell counts prior.
    """

    if mutation_object is True:
        trajectories = data.data
    else:
        trajectories = data

    ind_likelihood = []
    for traj in trajectories:
        int_s = []
        for s in fitness_prior[0, :]:
            int_N_w = []
            for N_w in N_w_prior[0, :]:
                int_N_w.append(
                    bd_process_conditional_likelihood_s_N(traj,
                                                        s=s, N_w=N_w)
                    )
            int_s.append(np.trapz(x=N_w_prior[0, :],
                                y=int_N_w*N_w_prior[1, :]))

        marginalised_likelihood = np.trapz(x=fitness_prior[0, :],
                                        y=int_s*fitness_prior[1, :])
        ind_likelihood.append(marginalised_likelihood)
    
    mutation_prob = np.product(ind_likelihood)

    if mutation_object is True:
        # return updated model_comparison object 
        data.bd_prob = mutation_prob
        return data
    else:
        # return marginalised likelihood.
        return mutation_prob

#endregion
# =============================================================================
#region: Hidden Markov
# =============================================================================
def y_1_cond(x_range, data, N_w):

    prob_init = binom.pmf(k=data.AO,
                    n=data.DP,
                    p=distributions.clone_size_to_vaf(x_range, N_w=N_w))

    normalisation = np.trapz(x=x_range, y=prob_init)

    normalised_prob_init = prob_init/normalisation

    return normalised_prob_init

def y_k_cond(x, data, recursive_prob, previous_x_range, fitness, N_w):
    
    data_observation_prob = binom.pmf(k=data.AO,
                                      n=data.DP,
                                      p=distributions.clone_size_to_vaf(x, N_w=N_w))

    # Compute the transition probability of the hidden process for all initial clone sizes
    # Extract proportion and trials NB approximating BD process
    nb_p, nb_n = distributions.BD_parametrization(
        init_size=previous_x_range, s=fitness, delta_t=3)
        
    # Compute associated BD probabilities
    next_cs_prob = nbinom.pmf(k=x, n=nb_n, p=nb_p)
    marginal = np.trapz(x=previous_x_range, y=recursive_prob*next_cs_prob)

    return data_observation_prob*marginal

def find_range (traj, k, bd_init_ranges, fitness, N_w, alpha = 0.01, resolution = 25):

    # BD RANGE
    # Find clone sizes with non-zero probability according to the hidden process
    # Extract proportion and trials NB approximating BD process
    hidden_range_total = [] 
    for i, cs_range in enumerate(bd_init_ranges):
        # compute time-steps since observation
        time_steps = k-i
        nb_p, nb_n = distributions.BD_parametrization(
            init_size=cs_range, s=fitness, delta_t=time_steps*3)

        # Compute minimum probable clone size
        min_clone_size = nbinom.ppf(q=alpha,
                                n=nb_n[0], p=nb_p[0])

        # Compute maximum probable clone size
        max_clone_size = nbinom.ppf(q=1-alpha,
                                n=nb_n[1], p=nb_p[1])

        # create a list of sampled clone sizes
        hidden_range = list(np.linspace(min_clone_size, max_clone_size, resolution, dtype='int'))
        hidden_range_total += hidden_range
    
    # Merge list of sampled clone sizes
    hidden_range_total = np.array(hidden_range_total)
    
    # OBSERVATION RANGE
    # Use Beta function to find range of clone sizes near observation
    beta_p_conf_int = beta.ppf(q=[alpha, 1-alpha],
                               a=traj.iloc[k].AO+1,
                               b=(traj.iloc[k].DP
                                  - traj.iloc[k].AO
                                  + 1))

    # List of binomial p ranges and clone_size_ranges
    beta_p_range = np.linspace(beta_p_conf_int[0],
                               min(beta_p_conf_int[1], 0.5),
                               resolution)

    # Transform p values to clone sizes
    observation_range = distributions.vaf_to_clone_size(v=beta_p_range,
                                                        N_w=N_w)
    
    bd_init_ranges.append((observation_range[0], observation_range[-1]))

    next_range = concatenate((observation_range, hidden_range_total))
    next_range.sort(kind='mergesort')

    return next_range       


def hidden_markov_conditional_s_N(traj, fitness=0, N_w=50_000, resolution=25, alpha=0.01):
    """Compute the conditional probability of the hidden markov chain model, marginalised
    over all possible hidden paths"""
    
   # Step 1: compute range of reasonable initial clone sizes based on data alone
   # Beta function
    beta_p_conf_int = beta.ppf(q=[alpha, 1-alpha],
                              a=traj.iloc[0].AO+1,
                              b=(traj.iloc[0].DP
                                 - traj.iloc[0].AO
                                 + 1))

    # List of binomial p ranges and clone_size_ranges
    beta_p_range = np.linspace(beta_p_conf_int[0],
                        min(beta_p_conf_int[1], 0.5),
                        resolution)

    # Transform p values to clone sizes
    init_range = distributions.vaf_to_clone_size(v=beta_p_range,
                                                N_w=N_w)

    # initialise clone size ranges of starting birth and death processes
    bd_init_ranges = [(init_range[0], init_range[-1])]

    # Compute first term of recursion
    recursive_prob = y_1_cond(init_range, traj.iloc[0], N_w)

    previous_x_range = init_range

    for i in range(1,len(traj)):
        
        next_x_range = find_range (traj, i, bd_init_ranges, fitness, N_w, alpha, resolution)

        recursive_prob = np.array([y_k_cond(x, traj.iloc[i], recursive_prob, previous_x_range,
                                            fitness, N_w) for x in next_x_range])
        previous_x_range = next_x_range

    probability = np.trapz(x=previous_x_range, y=recursive_prob)
    return probability


def hidden_markov_model_probability(data,
                                 fitness_prior=flat_fitness_prior,
                                 N_w_prior=flat_N_w_prior,
                                 mutation_object=True):
    """Marginalise the conditional likelihood of observing a
    given time-series data.

    Parameters:
    -----------
    trajectory: Pandas DataFrame. Time-series dataset including columns:
                'AO, 'DP' and 'age'.
    fitness_prior: array. 2-D numpy array containing x and y values of
                   a sampled fitness prior.
    N_w_prior: array. 2-D numpy array containing x and y values of
               a sampled wild-type stem cell counts prior.
    """

    if mutation_object is True:
        trajectories = data.data
    else:
        trajectories = data

    ind_likelihood = []
    for traj in trajectories:
        int_s = []
        for s in fitness_prior[0, :]:
            int_N_w = []
            for N_w in N_w_prior[0, :]:
                int_N_w.append(
                    hidden_markov_conditional_s_N(traj,
                                                  fitness=s, N_w=N_w)
                    )
            marginalised_N_w = np.trapz(x=N_w_prior[0, :],
                                y=int_N_w*N_w_prior[1, :])
            int_s.append(marginalised_N_w)

        marginalised_likelihood = np.trapz(x=fitness_prior[0, :],
                                        y=int_s*fitness_prior[1, :])
        ind_likelihood.append(marginalised_likelihood)
    
    mutation_prob = np.product(ind_likelihood)

    if mutation_object is True:
        # return updated model_comparison object 
        data.bd_prob = mutation_prob
        return data
    else:
        # return marginalised likelihood.
        return mutation_prob



#endregion
# =============================================================================
#region: Model comparison
# =============================================================================

def bayesian_model_comparison(data, hidden=True, binomial=True, bayes_factor_threshold=4):

    if len(data.data) > 16:
        
        # Set optimal model to artifact by default 
        data.optimal_model = 'artifact'

        # Probabilities need to be set for future plots
        data.betabinom_artifact_prob = 1
        data.binom_artifact_prob = 1
        data.bd_prob = 0.1

        return data


    if len(data.data) > 1:
        data = betabinom_artifact_model_probability(
                    data,
                    prior_p=betabinom_p_prior[0],
                    prior_beta=betabinom_beta_prior[0])
    else:
        data = betabinom_artifact_model_probability(
                    data,
                    prior_p=betabinom_p_prior[1],
                    prior_beta=betabinom_beta_prior[1])
    
    data = binom_artifact_model_probability(data)
    
    if hidden is True:
        data = hidden_markov_model_probability(data)
    else:
        data = bd_process_model_probability(data)

    data.compute_optimal_model(binomial=binomial,
                               bayes_factor_threshold=bayes_factor_threshold)
    
    return data