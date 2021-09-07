# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-08-23
# =============================================================================
"""Create and export in JSON format the 3 to 1 letter amino acid dictionary.
"""
# =============================================================================
# Imports
# =============================================================================
import json
import os
# =============================================================================
# Main
# =============================================================================

# create dictionary assigning a color to each variant type
amino_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
              'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
              'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
              'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M',
              'Ter': '*'}

# Create path for exporting
path = '../../Resources/'
if not os.path.exists(path):
    os.makedirs(path)

# convert dictionary into string
# using json.dumps()
with open('../../Resources/amino_dict.json', 'w') as fp:
    json.dump(amino_dict, fp)
