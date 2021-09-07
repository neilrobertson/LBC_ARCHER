import numpy as np
from lifelines import CoxPHFitter
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def survival_analysis(survival_data, keep_columns, cohort):
    # Select cohort and columns
    cox_data = survival_data[
        survival_data.cohort.isin(cohort)][
            keep_columns + ['days_from_wave1', 'dead']]
    # Exclude columns not used as covariates and filter for nan values
    cox_data = cox_data.dropna()

    # normalise columns used for regression
    for column in keep_columns:
        data = cox_data[column] - np.mean(cox_data[column])
        data = data/np.std(data)
        cox_data[column] = data

    # Train Cox proportional hazard model
    cph = CoxPHFitter()
    cph.fit(cox_data, duration_col='days_from_wave1', event_col='dead')

    # access the individual results using cph.summary
    cph.print_summary()

    return cph


def plot_hr_analysis(model, covariate):
    fig = make_subplots(rows=1, cols=2, column_widths=[0.3, 0.7],
                        subplot_titles=('Estimated hazard ratio',
                                        'Survival stratification'))

    fig.add_trace(
        go.Scatter(
            y=[model.hazard_ratios_[0]],
            x=[covariate],
            marker_symbol='diamond',
            marker_size=15,
            showlegend=False,
            error_y=dict(
                type='data',
                symmetric=False,
                array=(np.exp(np.array(model.confidence_intervals_)[:, 1])
                       - model.hazard_ratios_[0]),
                arrayminus=(model.hazard_ratios_[0]
                            - np.exp(np.array(
                                model.confidence_intervals_)[:, 0])))),
        row=1, col=1)

    # Plot covariate effect
    for covariate in model.params_.index:
        values = [-2, 0, 2]
        partial_ax = model.plot_partial_effects_on_outcome(
                        covariates=covariate,
                        values=values, cmap='coolwarm')
        partial_ax.get_figure()

        # add traces to figure
        fig.add_trace(
            go.Scatter(x=partial_ax.lines[1].get_xdata(),
                       y=partial_ax.lines[1].get_ydata(),
                       mode='lines', line=dict(dash='dash', shape='hv'),
                       name='Mean'), row=1, col=2)
        fig.add_trace(
            go.Scatter(x=partial_ax.lines[0].get_xdata(),
                       y=partial_ax.lines[0].get_ydata(),
                       mode='lines', line=dict(shape='hv'),
                       name='-2 SD'), row=1, col=2)
        fig.add_trace(
            go.Scatter(x=partial_ax.lines[2].get_xdata(),
                       y=partial_ax.lines[2].get_ydata(),
                       mode='lines', line=dict(shape='hv'),
                       name='2 SD'), row=1, col=2)

    fig.update_layout(template='simple_white',
                      title=f'Effect of {covariate} on survival',
                      legend_title_text=f'{covariate}')

    y_range_hazards = [np.floor(np.exp(np.array(
                        model.confidence_intervals_)))[0, 0],
                       np.ceil(np.exp(np.array(
                        model.confidence_intervals_)))[0, 1]]
    fig.update_yaxes(title_text="Hazard Ratio (95% CI)", range=y_range_hazards,
                     row=1, col=1, dtick=1)
    fig.update_yaxes(title_text="Survivors (proportion)",
                     row=1, col=2, dtick=0.2)
    fig.update_xaxes(title_text=covariate, showticklabels=False, tickvals=[0],
                     row=1, col=1)
    fig.update_xaxes(title_text="Years", row=1, col=2)
    return fig
