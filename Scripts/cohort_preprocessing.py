# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-11
# =============================================================================
"""Use LiFT model comparison to select all fit trajectories.
"""
# =============================================================================
#region Import local packages
# =============================================================================
# Append root directory to system's path
import sys
sys.path.append('../ARCH_package')

import basic

# =============================================================================
# Import general packages
# =============================================================================

import pandas as pd
import os

# =============================================================================
# Check environment being used
# =============================================================================
if os.environ['CONDA_DEFAULT_ENV'] == 'ARCH':
    print('Conda environment ARCH is active')

else:
    print('Conda environment ARCH is not active')
    sys.exit()

#endregion
# =============================================================================
#region Load cohort
# =============================================================================
print('Loading non-synonymous trajectories')
# Load non-synonymous dataset
df = pd.read_csv(r'../Datasets/LBC_ARCHER.1PCT_VAF.13Jan22.non-synonymous.IncludesZeros.tsv', sep='\t')
lbc = basic.load(df, export_name='LBC_non-synonymous', create_dict=True)

print('Loading synonymous trajectories')
# Load synonymous dataset
df = pd.read_csv(r'../Datasets/LBC_ARCHER.1PCT_VAF.13Jan22.synonymous.IncludesZeros.tsv', sep='\t')
syn = basic.load(df, export_name='LBC_synonymous', create_dict=False)

#endregion
