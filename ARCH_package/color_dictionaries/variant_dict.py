# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-08-23
# =============================================================================
"""Create and export in JSON format a dictionary assigning colors to each
variant type observed in the LBC.
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
var_dict = {"3'Flank": '#00BBDA',
            "5'Flank": '#E18A00',
            "5'UTR": '#BE9C00',
            'Frame_Shift_Del': '#00BDC2',
            'Frame_Shift_Ins': '#24B700',
            'In_Frame_Ins': '#00BE70',
            'Missense_Mutation': '#00C1AB',
            'Nonsense_Mutation': '#F8766D',
            'Nonstop_Mutation': '#8B93FF',
            'Splice_Region': '#D575FE',
            'Splice_Site': '#F962DD',
            'Translation_Start_Site': '#92f962',
            'RNA': '#b5a32f',
            'In_Frame_Del': '#7c4131',
            'Intron': '#88ac53'}

# Create path for exporting
path = '../../Resources/'
if not os.path.exists(path):
    os.makedirs(path)

# convert dictionary into string
# using json.dumps()
with open('../../Resources/var_dict.json', 'w') as fp:
    json.dump(var_dict, fp)
