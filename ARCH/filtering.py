from ARCH import modelling
import numpy as np
import pandas as pd
from itertools import combinations

from sklearn.linear_model import LinearRegression

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
colors = px.colors.qualitative.Plotly


class filter():
    """Class object corresponding to a filtering of the full cohort"""

    def __init__(self, model=None, gradient_plot=None,
                 neutral_dist=None, gene_bar=None, gene_dict=None):
        self.model = model
        self.gradient_plot = gradient_plot
        self.neutral_dist = neutral_dist
        self.gene_bar = gene_bar
        self.gene_dict = gene_dict


class gradient_filter():
    """ Class object of filter distribution. Attributes:
    - mean_regr: scikit Linearr"""

    def __init__(self, mean_regr=None, var_regr=None, figures=None):
        self.mean_regr = mean_regr
        self.var_regr = var_regr
        self.figures = figures


def threshold_filter(cohort, threshold=0.02):
    """ Filter cohort according to a VAF threshold.
    Filter attribute is modified for each trajectory object in the cohort
    Returns:
    - filtered cohort for model fitting.
    - Table report of filter.
    """

    for part in cohort:
        for traj in part.trajectories:
            if max(traj.data.AF) >= threshold:
                traj.filter = True
            else:
                traj.filter = False

    model = model_cohort(cohort)
    fig = CHIP_plot(cohort)
    genes = gene_bar(model)

    filter_class = filter(model=model,
                          gradient_plot=fig,
                          gene_bar=genes[0],
                          gene_dict=genes[1])

    return filter_class


def neutral_filter(cohort, neutral, n_std=2):

    mean_regr, var_regr, figures = neutral_distribution(neutral, n_std)
    neutral_gradient_filter = gradient_filter(mean_regr=mean_regr,
                                              var_regr=var_regr,
                                              figures=figures)
    for part in cohort:
        for traj in part.trajectories:
            data = traj.data['AF'].iloc[0]
            threshold = mean_regr.predict(
                        np.c_[data]) + n_std*np.sqrt(
                                            var_regr.predict(np.c_[data]))
            if traj.gradient > threshold:
                traj.filter = True
            else:
                traj.filter = False

    model = model_cohort(cohort)
    fig = neutral_plot_inset(cohort, neutral, mean_regr, var_regr, n_std)
    genes = gene_bar(model)
    filter_class = filter(model=model,
                          gradient_plot=fig,
                          neutral_dist=neutral_gradient_filter,
                          gene_bar=genes[0],
                          gene_dict=genes[1])
    filter_class.neutral_fit = figures

    return filter_class


def neutral_distribution(syn, n_std=2,
                         mean_color=colors[2], var_color=colors[1],
                         bin_color='rgb(102,102,102)', bar_color='darkorange'):
    """ Find the distribution of gradients of neutral mutations
    Returns:
    - mean regressor
    - variance regressor
    - plots"""

    # Create a dataset with all possible combinations of timepoints
    # in each participant
    melt_syn = pd.DataFrame(columns=['AF', 'regularized_gradient'])

    for part in syn:
        for traj in part.trajectories:
            for combination in list(combinations(traj.data.index, 2)):
                data = traj.data.loc[list(combination)]
                gradient = np.diff(data.AF) / np.sqrt(np.diff(data.age))
                melt_syn = (
                    melt_syn.append({'AF': traj.data.loc[combination[0]]['AF'],
                                     'regularized_gradient': gradient[0]},
                                    ignore_index=True))

    # Exclude all mutations with VAF < 0.01
    # lack of data below this threshold alters the distribution of gradients
    filtered_melt_syn = melt_syn[melt_syn['AF'] > 0.01]
    # Create bins
    filtered_melt_syn['bin'] = pd.cut(filtered_melt_syn['AF'], 25).copy()
    # Variance and mean with bootstrap for linear regression and plotting

    bin_center = []     # track bin centers
    mean_dist = []      # mean distribution in each bin using bootstrap
    variance_dist = []  # variance distribution in each bin using bootstrap
    bin_size = []       # Number of trajectories in each bin

    for bin in filtered_melt_syn['bin'].unique():
        # filter dataset to obtain rows in bin
        VAF_bin = filtered_melt_syn[filtered_melt_syn['bin'] == bin]
        bin_size.append(len(VAF_bin))               # points in bin
        bin_center.append(VAF_bin['AF'].mean())     # compute bin center

        # Bootstrap to compute gradient mean and variance distributions
        n = len(VAF_bin)
        reps = 200   # samples drawn
        subsets = np.random.choice(VAF_bin['regularized_gradient'], (n, reps))
        mean_dist.append(subsets.mean(axis=0))
        variance_dist.append(subsets.var(axis=0))

    # compute mean and quantiles for each bin's mean
    mean = [x.mean() for x in mean_dist]
    full = [np.abs(np.quantile(x, [0.1, 0.9]) - np.mean(x)) for x in mean_dist]
    mean_low = [x[0] for x in full]
    mean_high = [x[1] for x in full]
    # points in bin
    # compute mean and quantiles for each bin's variance
    variance = [x.mean() for x in variance_dist]
    full = [np.abs(np.quantile(x, [0.1, 0.9]) - np.mean(x))
            for x in variance_dist]
    var_low = [x[0] for x in full]
    var_high = [x[1] for x in full]

    size_weight = np.array(bin_size)

    # Mean weighted model
    mean_regr = LinearRegression(fit_intercept=True)
    mean_regr.fit(np.c_[bin_center], mean, size_weight)
    mean_score = mean_regr.score(np.c_[bin_center], mean, size_weight)

    # Mean model figure
    mean_fig = px.scatter(x=bin_center, y=mean, hover_data=[bin_size],
                          error_y=mean_high, error_y_minus=mean_low)

    # trace mean linear fit
    x = [0, max(bin_center)]
    mean_fig.add_trace(
        go.Scatter(x=x, y=mean_regr.predict(np.c_[x]),
                   marker_color=mean_color,
                   name='linear regression fit'))
    mean_fig.update_layout(title=('Linear regression gradient ~ VAF <br>'
                                  f'R2 score: {round(mean_score,3)}'),
                           template='plotly_white')
    mean_fig.update_xaxes(title='VAF')
    mean_fig.update_yaxes(title='Gradient')

    # Variance weighted model
    var_regr = LinearRegression(fit_intercept=False)
    var_regr.fit(np.c_[bin_center], variance, size_weight)
    var_score = var_regr.score(np.c_[bin_center], variance, size_weight)

    # Variance model figure
    var_fig = px.scatter(x=bin_center, y=variance, hover_data=[bin_size],
                         error_y=var_high, error_y_minus=var_low)

    x = [0, max(bin_center)]
    var_fig.add_trace(
        go.Scatter(x=x, y=var_regr.predict(np.c_[x]),
                   marker_color=var_color,
                   name='Linear regression fit'))
    var_fig.update_layout(title=('Linear regression of variance ~ VAF <br>'
                                 f'R2 score: {round(var_score,3)}'),
                          template='plotly_white')
    var_fig.update_xaxes(title='VAF')
    var_fig.update_yaxes(title='Gradients Variance')

    # Create a dataframe with the regularized gradient
    # between first and last timepoint
    total_gradient = pd.DataFrame(columns=['AF', 'gradient'])
    for part in syn:
        for traj in part.trajectories:
            if len(traj.data) > 1:
                total_gradient = (
                    total_gradient.append({'AF': traj.data.iloc[0]['AF'],
                                           'gradient': traj.gradient},
                                          ignore_index=True))

    # Compute the percentage of succesfully filtered trajectories:

    x = total_gradient['AF']
    total_gradient['filter'] = total_gradient['gradient'] \
        - (mean_regr.predict(np.c_[x])
           + n_std*np.sqrt(var_regr.predict(np.c_[x])))
    percentage = 100 * len(total_gradient[total_gradient['filter'] < 0]) \
        / (len(total_gradient)-1)

    # Plot filter
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(x=total_gradient['AF'],
                   y=total_gradient['gradient'],
                   name='Synonymous mutations',
                   mode='markers'))
    # Confidence line
    x = np.linspace(0, max(total_gradient['AF']), 1000)
    y_std = mean_regr.predict(np.c_[x]) \
        + n_std*np.sqrt(var_regr.predict(np.c_[x]))
    y_mean = mean_regr.predict(np.c_[x])

    fig.add_trace(
        go.Scatter(x=x, y=y_mean,
                   name='Mean linear regression',
                   marker_color=mean_color))

    fig.add_trace(
        go.Scatter(x=x, y=y_std,
                   name='Neutral growth filter',
                   marker_color=var_color))

    fig.update_layout(title=(f'Neutral growth filter <br>Filters '
                             f'{round(percentage,1)}% of all '
                             f'synonymous mutations'),
                      template='plotly_white')
    fig.update_xaxes(title='VAF')
    fig.update_yaxes(title='Regularized gradient')

    # NGF stack plot
    # Start subplot figure with shared xaxis
    NGF_stack = make_subplots(rows=3, cols=1,
                              row_heights=[0.4, 0.4, 0.2],
                              shared_xaxes=True,
                              vertical_spacing=0.05)

    # Row 1 variance plot
    NGF_stack.add_trace(
        go.Scatter(
            x=bin_center,
            y=variance,
            mode='markers',
            marker_color=bin_color,
            showlegend=False,
            error_y=dict(
                type='data',
                symmetric=False,
                array=var_high,
                arrayminus=var_low)
            ),
        row=1, col=1)

    x = [0, max(bin_center)]
    NGF_stack.add_trace(
        go.Scatter(x=x,
                   y=var_regr.predict(np.c_[x]),
                   mode='lines',
                   marker_color=var_color,
                   showlegend=False),
        row=1, col=1)

    # Row 2 mean plot
    NGF_stack.add_trace(
        go.Scatter(
            x=bin_center,
            y=mean,
            mode='markers',
            marker_color=bin_color,
            showlegend=False,
            error_y=dict(
                type='data',
                symmetric=False,
                array=mean_high,
                arrayminus=mean_low)
            ),
        row=2, col=1)

    # trace mean linear fit
    x = [0, max(bin_center)]
    NGF_stack.add_trace(
        go.Scatter(
            x=x,
            y=mean_regr.predict(np.c_[x]),
            mode='lines',
            marker_color=mean_color,
            showlegend=False),
        row=2, col=1)

    # Compute bins for barplot
    bin_center = []
    bin_width = []
    bin_size = []
    for bin in filtered_melt_syn.bin.unique():
        bin_size.append(
            len(filtered_melt_syn[filtered_melt_syn['bin'] == bin]))
        bin_center.append((bin.left+bin.right)/2)
        bin_width.append(bin.right-bin.left)

    # Row 3 Plot bin histogram
    NGF_stack.add_trace(
        go.Bar(
            x=bin_center,
            y=bin_size,
            marker_color='Orange',
            width=bin_width,
            showlegend=False),
        row=3, col=1)

    NGF_stack.update_yaxes(row=1, col=1,
                           title_text='Gradient Variance',
                           range=[0, max(variance)])
    NGF_stack.update_yaxes(title_text='Gradient', row=2, col=1)
    NGF_stack.update_yaxes(title_text='Bin counts',
                           tickmode='linear', dtick=400, row=3, col=1)
    NGF_stack.update_xaxes(title_text="VAF", row=3, col=1)
    NGF_stack.update_layout(template='plotly_white')

    NGF_stack.update_layout(
        template='plotly_white',
        legend=dict(orientation="h",
                    yanchor="bottom", y=1.02,
                    xanchor="right", x=1))

    return mean_regr, var_regr, [fig, mean_fig, var_fig, NGF_stack]


def gene_bar(cohort, split=True, relative=False):
    """Bar plot with counts of filtered mutations by gene."""

    if split is True:
        """Bar plot with counts of filtered mutations by gene."""

        # Create a dictionary with all filtered genes
        gene_list = []
        for traj in cohort:
            gene_list.append(traj.mutation.split()[0])

        gene_dict = {element: 0 for element in set(gene_list)}
        gene_dict_increasing = {element: 0 for element in set(gene_list)}
        gene_dict_decreasing = {element: 0 for element in set(gene_list)}

        # update the counts for each gene
        for traj in cohort:
            gene_dict[traj.mutation.split()[0]] = \
                gene_dict[traj.mutation.split()[0]] + 1

            # Separate counts for increasing and decreasing trajectories
            gradient = traj.y[-1]-traj.y[0] > 0
            if gradient is True:
                gene_dict_increasing[traj.mutation.split()[0]] = \
                    gene_dict_increasing[traj.mutation.split()[0]] + 1
            else:
                gene_dict_decreasing[traj.mutation.split()[0]] = \
                    gene_dict_decreasing[traj.mutation.split()[0]] + 1

        # sort dictionary in descending order
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: item[1], reverse=True))

        # extract ordered list of genes and values
        genes = [key for key in gene_dict.keys()]
        values_increasing = [gene_dict_increasing[gene] for gene in genes]
        if relative is False:
            values_decreasing = [gene_dict_decreasing[gene] for gene in genes]
        else:
            values_decreasing = [-gene_dict_decreasing[gene] for gene in genes]

        total_inc = sum(values_increasing)
        total_dec = np.abs(sum(values_decreasing))

        # Bar plot
        fig = go.Figure()
        fig.add_trace(
            go.Bar(x=genes, y=values_increasing,
                   name=f'{total_inc} Increasing',
                   marker_color=colors[0]))
        fig.add_trace(
            go.Bar(x=genes, y=values_decreasing,
                   name=f'{total_dec} Decreasing',
                   marker_color=colors[0],
                   opacity=0.3))

        fig.update_layout(
            title='Gene distribution of filtered mutations',
            barmode='relative',
            template="simple_white",
            yaxis_title='Trajectory counts',
            xaxis_tickangle=-45)

    else:                   # unsplit barplot

        # Create a dictionary with all filtered genes
        gene_list = []
        for traj in cohort:
            gene_list.append(traj.mutation.split()[0])
        gene_dict = {element: 0 for element in set(gene_list)}

        # update the counts for each gene
        for traj in cohort:
            gene_dict[traj.mutation.split()[0]] = \
                gene_dict[traj.mutation.split()[0]] + 1
        # sort dictionary in descending order
        gene_dict = dict(sorted(gene_dict.items(),
                                key=lambda item: item[1], reverse=True))
        genes = [key for key in gene_dict.keys()]
        values = [value for value in gene_dict.values()]

        # Bar plot
        fig = go.Figure()
        fig.add_trace(go.Bar(x=genes, y=values))
        fig.update_layout(title='Gene distribution of filtered mutations',
                          template="simple_white",
                          xaxis_tickangle=-45)
    # Labels position
    fig.update_layout(legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95,
        ))
    return fig, gene_dict


def find_fit(cohort):
    """Find the first fit trajectory of the cohort"""
    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                return traj


def find_neutral(cohort):
    """Find the first fit trajectory of the cohort"""
    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is False:
                return traj


def CHIP_plot(cohort, threshold=0.02):
    """ Plot the selection of filtered variants.
    threshold: Default = 0.02.
    x_axis: VAF.
    y_axis: gradient.
    color: filter attribute.
    """

    fit_traj = find_fit(cohort)
    neutral_traj = find_neutral(cohort)

    def legend_group(filter):
        if filter is True:
            return 'Fit'
        else:
            return 'Neutral'

    def color_group(filter):
        if filter is True:
            return colors[0]
        else:
            return colors[4]

    fig = go.Figure()
    # Create trajectories for legend groups
    # Fit trajectories
    if fit_traj is not None:
        fig.add_trace(
            go.Scatter(x=[fit_traj.data.AF.iloc[0]],
                       y=[fit_traj.gradient],
                       name='CHIP mutations',
                       mode='markers',
                       marker=dict(color=color_group(fit_traj.filter),
                                   size=5),
                       legendgroup=legend_group(fit_traj.filter),
                       showlegend=True))

    # Neutral trajectories
    if neutral_traj is not None:
        fig.add_trace(
            go.Scatter(x=[neutral_traj.data.AF.iloc[0]],
                       y=[neutral_traj.gradient],
                       name='non-CHIP mutations',
                       mode='markers',
                       marker=dict(color=color_group(neutral_traj.filter),
                                   size=5),
                       legendgroup=legend_group(neutral_traj.filter),
                       showlegend=True))

    for part in cohort:
        for traj in part.trajectories:
            fig.add_trace(
                go.Scatter(x=[traj.data.AF.iloc[0]],
                           y=[traj.gradient],
                           mode='markers',
                           marker=dict(color=color_group(traj.filter),
                                       size=5),
                           legendgroup=legend_group(traj.filter),
                           showlegend=False))
    # Add vertical delimitation of threshold
    fig.add_trace(
        go.Scatter(x=[0.02, 0.02], y=[-0.05, 0.1],
                   name=f'{threshold} VAF threshold',
                   mode='lines',
                   line=dict(color=colors[2], width=3, dash='dash')))
    fig.update_layout(title=('Filtering trajectories '
                             'achieving VAF > 0.02 over timespan'))
    fig.update_xaxes(title='VAF', range=[0, 0.35])
    fig.update_yaxes(title='Gradient', range=[-0.03, 0.08])

    fig.update_layout(template='plotly_white')

    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.95,
        xanchor="right",
        x=0.95,
    ))

    return fig


def CHIP_plot_inset(cohort, threshold=0.02,
                    xrange_2=[0.01, 0.03], yrange_2=[-0.01, 0.005]):
    """ Plot the selection of filtered variants.
    threshold: Default = 0.02.
    x_axis: VAF.
    y_axis: gradient.
    color: filter attribute.
    """

    # Find fit and netural trajectories to create group
    fit_traj = find_fit(cohort)
    neutral_traj = find_neutral(cohort)

    def legend_group(filter):
        """Legend group assigns a legend group based on filter attribute"""
        if filter is True:
            return 'Fit'
        else:
            return 'Neutral'

    def color_group(filter):
        """ returns a different color based on filter attribute"""
        if filter is True:
            return colors[0]
        else:
            return colors[4]

    # Create figure with inset
    fig = make_subplots(insets=[{'cell': (1, 1), 'l': 0.7, 'b': 0.3}])

    # Create trajectories for legend groups
    # Fit trajectories
    if fit_traj is not None:
        fig.add_trace(
            go.Scatter(x=[fit_traj.data.AF.iloc[0]],
                       y=[fit_traj.gradient],
                       name='CHIP mutations',
                       mode='markers',
                       marker=dict(color=color_group(fit_traj.filter),
                                   size=5),
                       legendgroup=legend_group(fit_traj.filter),
                       showlegend=True))

    # Neutral trajectories
    if neutral_traj is not None:
        fig.add_trace(
            go.Scatter(x=[neutral_traj.data.AF.iloc[0]],
                       y=[neutral_traj.gradient],
                       name='non-CHIP mutations',
                       mode='markers',
                       marker=dict(color=color_group(neutral_traj.filter),
                                   size=5),
                       legendgroup=legend_group(neutral_traj.filter),
                       showlegend=True))

    for part in cohort:
        for traj in part.trajectories:
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            mode='markers',
                            marker=dict(color=color_group(traj.filter),
                                        size=5),
                            legendgroup=legend_group(traj.filter),
                            showlegend=False)

            # add inset scatter plot
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            xaxis='x2', yaxis='y2',
                            mode='markers',
                            marker=dict(color=color_group(traj.filter),
                                        size=3),
                            legendgroup=legend_group(traj.filter),
                            showlegend=False)

    # Add vertical delimitation of threshold
    fig.add_trace(
        go.Scatter(x=[0.02, 0.02], y=[-0.05, 0.1],
                   name=f'{threshold} VAF threshold',
                   mode='lines',
                   line=dict(color=colors[2], width=3, dash='dash')))

    # Add rectangle in zoomed area
    fig.add_shape(type="rect",
                  x0=xrange_2[0], y0=yrange_2[0],
                  x1=xrange_2[1], y1=yrange_2[1],
                  line=dict(color="Black", width=1))

    fig.update_layout(
        title=None,
        template='plotly_white',
        xaxis=dict(
            title='VAF',
            range=[0, 0.35]),
        yaxis=dict(
            title='Gradient',
            range=[-0.017, 0.08]),
        yaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=yrange_2,
            domain=[0.55, 1],
            showline=None,
            linewidth=1,
            linecolor='black',
            mirror=True),
        xaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=xrange_2,
            domain=[0.6, 1],
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True)
    )

    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0
    ))

    return fig


def neutral_plot(cohort, neutral, mean_regr, var_regr, n_std=2):
    # Compute the standard deviation of absolute regularized total gradients
    fig = go.Figure()

    # Add traces for labeling
    fig.add_trace(go.Scatter(x=[0], y=[0],
                             line=dict(color=colors[0]),
                             mode='markers',
                             legendgroup="fit",
                             name="Fit variants"
                             ))

    fig.add_trace(go.Scatter(x=[0], y=[0],
                             line=dict(color=colors[3]),
                             mode='markers',
                             legendgroup="non-fit",
                             name="Neutral variants"
                             ))

    fig.add_trace(go.Scatter(x=[0],
                             y=[0],
                             line=dict(color='Orange'),
                             mode='markers',
                             legendgroup="synonymous",
                             name="Synonymous variants"
                             ))

    # Add traces
    for part in cohort:
        for traj in part.trajectories:
            if traj.filter is True:
                fig.add_trace(
                    go.Scatter(x=[traj.data.AF.iloc[0]],
                               y=[traj.gradient],
                               line=dict(color=colors[0]),
                               legendgroup="fit",
                               showlegend=False,
                               hovertemplate=f'Mutation: {traj.mutation}'
                               ))
            else:
                fig.add_trace(
                    go.Scatter(x=[traj.data.AF.iloc[0]],
                               y=[traj.gradient],
                               line=dict(color=colors[3]),
                               legendgroup="non-fit",
                               showlegend=False,
                               hovertemplate=f'Mutation: {traj.mutation}'
                               ))
    for part in neutral:
        for traj in part.trajectories:
            fig.add_trace(
                go.Scatter(x=[traj.data.AF.iloc[0]],
                           y=[traj.gradient],
                           line=dict(color=colors[4]),
                           marker=dict(size=3),
                           legendgroup="synonymous",
                           showlegend=False,
                           hovertemplate=f'Mutation: {traj.mutation}'
                           ))

    x = np.linspace(0, 0.5, 1000)
    y_std = mean_regr.predict(np.c_[x]) \
        + n_std*np.sqrt(var_regr.predict(np.c_[x]))

    y_mean = mean_regr.predict(np.c_[x])

    fig.add_trace(
        go.Scatter(x=x, y=y_std,
                   name='Neutral growth filter',
                   marker=dict(color=colors[1])))
    fig.add_trace(
        go.Scatter(x=x, y=y_mean,
                   name='Neutral growth mean',
                   marker=dict(color=colors[2])))

    fig.update_layout(title=('Filtering according to '
                             'neutral growth distribution'))
    fig.update_xaxes(title='VAF', range=[0, 0.35])
    fig.update_yaxes(title='Gradient', range=[-0.03, 0.08])

    fig.update_layout(template='plotly_white')
    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0
    ))

    return fig


def neutral_plot_inset(cohort, neutral, mean_regr, var_regr,
                       n_std=2, xrange_2=[0, 0.03], yrange_2=[-0.01, 0.01]):

    # Find fit and netural trajectories to create group
    fit_traj = find_fit(cohort)
    neutral_traj = find_neutral(cohort)

    def legend_group(filter):
        """Legend group assigns a legend group based on filter attribute"""
        if filter is True:
            return 'Fit'
        else:
            return 'Neutral'

    def color_group(filter):
        """ returns a different color based on filter attribute"""
        if filter is True:
            return colors[0]
        else:
            return colors[3]

    # Create figure with inset
    fig = make_subplots(insets=[{'cell': (1, 1), 'l': 0.7, 'b': 0.3}])

    # Create trajectories for legend groups
    # Fit trajectories
    if fit_traj is not None:
        fig.add_trace(
            go.Scatter(x=[fit_traj.data.AF.iloc[0]],
                       y=[fit_traj.gradient],
                       name='Fit',
                       mode='markers',
                       marker=dict(color=color_group(fit_traj.filter),
                                   size=5),
                       legendgroup=legend_group(fit_traj.filter),
                       showlegend=True))

    # Neutral trajectories label
    if neutral_traj is not None:
        fig.add_trace(
            go.Scatter(x=[neutral_traj.data.AF.iloc[0]],
                       y=[neutral_traj.gradient],
                       name='Neutral',
                       mode='markers',
                       marker=dict(color=color_group(neutral_traj.filter),
                                   size=5),
                       legendgroup='synonymous',
                       showlegend=True))

    # Synonymous trajectories label:
    synonymous = neutral[0].trajectories[0]
    fig.add_trace(
        go.Scatter(x=[synonymous.data.AF.iloc[0]],
                   y=[synonymous.gradient],
                   name='Synonymous',
                   mode='markers',
                   marker=dict(color='Orange',
                               size=5),
                   legendgroup='synonymous',
                   showlegend=True))

    for part in cohort:
        for traj in part.trajectories:
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            mode='markers',
                            marker=dict(color=color_group(traj.filter),
                                        size=5),
                            legendgroup=legend_group(traj.filter),
                            showlegend=False)

            # add inset scatter plot
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            xaxis='x2', yaxis='y2',
                            mode='markers',
                            marker=dict(color=color_group(traj.filter),
                                        size=3),
                            legendgroup=legend_group(traj.filter),
                            showlegend=False)
    for part in neutral:
        for traj in part.trajectories:
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            mode='markers',
                            marker=dict(color='Orange',
                                        size=3),
                            legendgroup='synonymous',
                            showlegend=False)

            # add inset scatter plot
            fig.add_scatter(x=[traj.data.AF.iloc[0]],
                            y=[traj.gradient],
                            xaxis='x2', yaxis='y2',
                            mode='markers',
                            marker=dict(color='Orange',
                                        size=3),
                            legendgroup='synonymous',
                            showlegend=False)

    # Add mean and filter lines
    x = np.linspace(0, 0.5, 1000)
    y_std = mean_regr.predict(np.c_[x]) \
        + n_std*np.sqrt(var_regr.predict(np.c_[x]))
    y_mean = mean_regr.predict(np.c_[x])
    fig.add_trace(
        go.Scatter(x=x, y=y_std,
                   name='Neutral growth filter',
                   xaxis='x2', yaxis='y2',
                   legendgroup='filter',
                   marker=dict(color=colors[1])))
    fig.add_trace(
        go.Scatter(x=x, y=y_mean,
                   name='Neutral growth mean',
                   xaxis='x2', yaxis='y2',
                   legendgroup='mean',
                   marker=dict(color=colors[2])))
    # Add mean and filter lines to inset
    fig.add_trace(
        go.Scatter(x=x, y=y_std,
                   name='Neutral growth filter',
                   legendgroup='filter',
                   showlegend=False,
                   marker=dict(color=colors[1])))
    fig.add_trace(
        go.Scatter(x=x, y=y_mean,
                   name='Neutral growth mean',
                   legendgroup='mean',
                   showlegend=False,
                   marker=dict(color=colors[2])))

    # Add rectangle in zoomed area
    fig.add_shape(type="rect",
                  x0=xrange_2[0], y0=yrange_2[0],
                  x1=xrange_2[1], y1=yrange_2[1],
                  line=dict(color="Black", width=1))

    # Update inset layout
    fig.update_layout(
        title=None,
        template='plotly_white',
        xaxis=dict(
            title='VAF',
            range=[0, 0.35]),
        yaxis=dict(
            title='Gradient',
            range=[-0.017, 0.08]),
        yaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=yrange_2,
            domain=[0.55, 1],
            showline=None,
            linewidth=1,
            linecolor='black',
            mirror=True),
        xaxis2=dict(
            title=None,
            tickfont=dict(size=10),
            range=xrange_2,
            domain=[0.6, 1],
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True)
    )

    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0
    ))

    return fig


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


def participant_filter_variant(part):
    """ Plot a participant's profile using different colors for
    fit and non-fit mutations"""

    var_types = ["3'Flank", "5'Flank", "5'UTR", "Frame_Shift_Del",
                 "Frame_Shift_Ins", "In_Frame_Ins", "Missense_Mutation",
                 "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region",
                 "Splice_Site"]
    var_colors = ["#00BBDA", "#E18A00", "#BE9C00", "#8CAB00", "#24B700",
                  "#00BE70", "#00C1AB", "#F8766D", "#8B93FF", "#D575FE",
                  "#F962DD"]

    var_zip = zip(var_types, var_colors)
    var_dict = dict(var_zip)

    fig = go.Figure()

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
                   marker=dict(color='rgb(204,204,204)'),
                   legendgroup='non-fit',
                   showlegend=True))

    for traj in part.trajectories:
        # Plot trajectories  using different colors according to traj.filter
        if traj.filter is False:
            fig.add_trace(
                go.Scatter(x=traj.data.age, y=traj.data.AF,
                           name='non-fit',
                           mode='lines',
                           marker=dict(color='rgb(204,204,204)'),
                           legendgroup='non-fit',
                           showlegend=False))

    for traj in part.trajectories:
        if traj.filter is True:
            fig.add_trace(
                go.Scatter(x=traj.data.age, y=traj.data.AF,
                           name=traj.mutation,
                           mode='lines',
                           marker=dict(color=var_dict[traj.variant_class]),
                           legendgroup='fit',
                           showlegend=True))

    # Update layout
    fig.update_layout(title=f'Participant {part.id}',
                      xaxis_title='Age',
                      yaxis_title='VAF')

    return fig


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
                   marker=dict(color='rgb(204,204,204)'),
                   legendgroup='non-fit',
                   showlegend=True))

    for traj in part.trajectories:
        # Plot trajectories  using different colors according to traj.filter
        if traj.filter is False:
            fig.add_trace(
                go.Scatter(x=traj.data.age, y=traj.data.AF,
                           name='non-fit',
                           mode='lines',
                           marker=dict(color='rgb(204,204,204)'),
                           legendgroup='non-fit',
                           showlegend=False))

    for traj in part.trajectories:
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
