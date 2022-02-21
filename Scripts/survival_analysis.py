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

import survival

# =============================================================================
# Import global packages
# =============================================================================

import os

# =============================================================================
# Create path for exporting
# =============================================================================

# Create path for exporting
path = '../Results/Survival analysis/'
if not os.path.exists(path):
    os.makedirs(path)

#endregion
# =============================================================================
#region Survival analysis and plots
# =============================================================================

# Create path for exporting
path = '../Results/Survival analysis/'
if not os.path.exists(path):
    os.makedirs(path)

# Cohorts analysed

# Parameters analysed
survival_parameters = [['fitness_vaf'], ['max_vaf']]

for cohort in [["36"],["21"]]:
    for survival_columns in survival_parameters:
        # Fit CoxPh regression model to each collection of parameters 
        cph = survival.survival_analysis(survival_columns, cohort)
        # Export results
        cph.summary.to_csv(path + cohort[0]+ f"_summary_analysis_{survival_columns[0]}.csv")

        fig = survival.plot_hr_analysis(cph, covariate=survival_columns[0])
        fig.update_layout(title=cohort[0])
        fig.write_image(path +cohort[0]+ f'_survival_analysis_{survival_columns[0]}.svg')
        fig.write_image(path +cohort[0]+ f'_survival_analysis_{survival_columns[0]}.png', scale=10)
