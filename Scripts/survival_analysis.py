import os
import pandas as pd

import sys
sys.path.append('../ARCH_package')

import survival

# Create path for exporting
path = '../Results/Survival analysis/'
if not os.path.exists(path):
    os.makedirs(path)

survival_data = pd.read_csv('../Datasets/survival_data.csv')
cohort = [21]
survival_columns = ['Max initial vaf']

cph = survival.survival_analysis(survival_data, survival_columns, cohort)
fig = survival.plot_hr_analysis(cph, covariate=survival_columns[0])
fig.write_image(path + 'LBC21_init_vaf.svg', width=1000)

cohort = [21, 36]
survival_columns = ['Gradient']
cph = survival.survival_analysis(survival_data, survival_columns, cohort)
fig = survival.plot_hr_analysis(cph, covariate=survival_columns[0])
fig.write_image(path + 'LBC_gradient.svg', width=1000)
