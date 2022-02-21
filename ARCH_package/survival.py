# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2022-02
# =============================================================================
"""Extract data and create functions for survival analysis study of the LBC
"""
# =============================================================================
#region Initialisation
# =============================================================================

# =============================================================================
# Import local packages
# =============================================================================
# Append root directory to system's path

import sys
sys.path.append("../ARCH_package")

# =============================================================================
# Import global packages
# =============================================================================

import dill
import pandas as pd
import numpy as np

from lifelines import CoxPHFitter
import plotly.graph_objects as go
from plotly.subplots import make_subplots

#endregion
# =============================================================================
#region Import datasets
# =============================================================================

# Import lbc metadata
df21 = pd.read_csv('../Datasets/LBC1921_ExtraVariables.csv')
df36 = pd.read_csv('../Datasets/LBC1936_ExtraVariables.csv')
df36['dead'] = df36['dead'].replace('.', 0)
extended = pd.read_csv('../Datasets/lbc_meta.csv')


# Import non-synonymous mutations as exported in LiFT.py
with open('../Exports/LBC_non-synonymous_LiFTed.dill', 'rb') as infile:
    lbc_filtered = dill.load(infile)

lbc_filtered = [part for part in lbc_filtered if len(part.trajectories)>0]

# Import non-synonymous mutations as exported in LiFT.py
with open('../Exports/LBC_non-synonymous_LiFT_fitted.dill', 'rb') as infile:
    lbc_fitted = dill.load(infile)

#endregion
# =============================================================================
#region Define classes and functions
# =============================================================================

class survival_participant():
    def __init__(self, id):
        self.id = id


def from_wave_1(aux):
    if aux.cohort == "21":
        aux.years_from_wave_1 = (aux.age_days - 82*365.25)/365.25
    else:
        aux.years_from_wave_1 = (aux.age_days - 70*365.25)/365.25


def find_sex (aux):
    sex = extended[extended.ID == aux.id].sex.unique()[0]
    
    if sex == 'M':
        return 0
    else:
        return 1


def survival_analysis(keep_columns, cohort):
    # Select cohort and columns
    cox_data = survival_data[survival_data.cohort.isin(cohort)][keep_columns + ['years_from_wave_1', 'dead']]
    # Exclude columns not used as covariates and filter for nan values
    cox_data = cox_data.dropna()

    # normalise columns used for regression
    for column in keep_columns:
        data = cox_data[column] - np.mean(cox_data[column])
        data = data/np.std(data)
        cox_data[column] = data

    # Train Cox proportional hazard model
    cph = CoxPHFitter()
    cph.fit(cox_data, duration_col='years_from_wave_1', event_col='dead') 
    
    return cph


def plot_hr_analysis(model, covariate):
    fig = make_subplots(rows=1, cols=2, column_widths=[0.3, 0.7],
                        subplot_titles=(f'Estimated hazard ratio', f'Survival stratification'))

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
                array=np.exp(np.array(model.confidence_intervals_)[:,1])-model.hazard_ratios_[0],
                arrayminus=model.hazard_ratios_[0]-np.exp(np.array(model.confidence_intervals_)[:,0]))
            ), row=1, col=1)

    # Plot covariate effect
    for covariate in model.params_.index:
        values =[-2, 0 , 2]
        partial_ax = model.plot_partial_effects_on_outcome(covariates=covariate, values=values, cmap='coolwarm')
        partial_ax.get_figure()

        #add traces to figure
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

    y_range_hazards = [np.floor(np.exp(np.array(model.confidence_intervals_)))[0,0], np.ceil(np.exp(np.array(model.confidence_intervals_)))[0,1]] 
    fig.update_yaxes(title_text="Hazard Ratio (95% CI)", range=y_range_hazards, row=1, col=1, dtick=1) 
    fig.update_yaxes(title_text="Survivors (proportion)", row=1, col=2, dtick=0.2) 
    fig.update_xaxes(title_text=covariate, showticklabels=False,tickvals=[0], row=1, col=1)
    fig.update_xaxes(title_text="Years", row=1, col=2)
    return fig

#endregion
# =============================================================================
#region Create Survival analysis dataset
# =============================================================================

part_list = []
for part in lbc_filtered:
    aux = survival_participant(part.id)
    aux.sex = find_sex(aux)
    # Append max vaf
    aux.max_vaf = max([max(traj.data.AF) for traj in part.trajectories])
    
    # append max fitness
    aux.fitness = 0
    aux.fitness_vaf = 0
    for part_fitted in lbc_fitted:
        if part_fitted.id == part.id:
            aux.fitness = max([traj.fitness for traj in part_fitted.trajectories])
            for traj in part_fitted.trajectories:
                aux.fitness_vaf = max(aux.fitness_vaf, max(traj.data.AF*traj.fitness))
    
    if (part.cohort == 'LBC21') and (part.id in list(df21.studyno)):
        
        aux.cohort = '21'
        part_row = df21[df21.studyno == part.id]
        aux.dead = part_row.dead.values[0]
        if aux.dead == 0:
            try:
                aux.age_days = float(part_row.agedaysApx_LastCensor.values[0])
                from_wave_1(aux)
            except:
                aux.age_days=None
        else:
            try:
                aux.age_days =  float(part_row.agedays_death.values[0])
                from_wave_1(aux)

            except:
                aux.age_days=None
        part_list.append(aux)
                
    if (part.cohort == 'LBC36') and (part.id in list(df36.lbc36no)):
        aux.cohort = '36'
        part_row = df36[df36.lbc36no == part.id]
        aux.dead = part_row.dead.values[0]
        if aux.dead == 0:
            try:
                aux.age_days = float(part_row.AgedaysApx_LastCensor.values[0])
                from_wave_1(aux)
            except:
                aux.age_days=None
        else:
            try:
                aux.age_days =  float(part_row.agedays_death.values[0])
                from_wave_1(aux)
            except:
                aux.age_days=None
    
        part_list.append(aux)

survival_data = pd.DataFrame()
for part in part_list:
    survival_data = survival_data.append(part.__dict__, ignore_index=True)
    
survival_data = survival_data.dropna()

#endregion