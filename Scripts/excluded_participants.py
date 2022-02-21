# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-11-26
# =============================================================================
"""Create and export in JSON format a list of participants excluded from the
study.
"""
# =============================================================================
# Imports
# =============================================================================
import json
import os
# =============================================================================
# Main
# =============================================================================

# Create a list of participants with a clinical record history
# that might have affected their clonal dynamics.
clinical_records = ['LBC360021', 'LBC360725']

# Create a list of participants with possible sample contamination
sample_contamination = ['LBC361096']

# Create the full list of participants excluded from the study
excluded_samples = clinical_records + sample_contamination

# Create path for exporting
path = '../../Resources/'
if not os.path.exists(path):
    os.makedirs(path)

# convert list into string
# using json.dumps()
with open('../Resources/excluded_samples.json', 'w') as fp:
    json.dump(excluded_samples, fp)
