# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================
"""Module to create a filtering of a cohort of participants based on their
genetic trajectories.
"""
# =============================================================================
# Import local packages
# =============================================================================
import plot
import modelling
# =============================================================================
# Imports
# =============================================================================
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn.linear_model import LinearRegression

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio

colors = px.colors.qualitative.Plotly
pio.templates.default = "simple_white"

# =============================================================================
# Create local classes
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


# =============================================================================
# Common filter functions
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

# =============================================================================
# Threshold filter functions
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
    fig = plot.CHIP_plot(cohort)
    # Bar plot of genes selected through filteringfigure.
    genes = plot.gene_bar(model)

    # Create a filter class object
    filter_class = filter(model=model,
                          gradient_plot=fig,
                          gene_bar=genes[0],
                          gene_dict=genes[1])

    return filter_class

# =============================================================================
# Threshold filter functions
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
                      template='plotly_white')
    fig.update_xaxes(title='VAF')
    fig.update_yaxes(title=yaxis)

    return regr, fig
