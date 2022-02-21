# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-08-23
# =============================================================================
"""Create and export in JSON format a dictionary assigning colors to each
gene observed in the LBC.
"""
# =============================================================================
# Imports
# =============================================================================
import json
import os
import numpy as np
import plotly.express as px
# =============================================================================
# Main
# =============================================================================


def create_dict(variant, colors=px.colors.qualitative.D3):
    """ Create and export in JSON format a dictionary assigning a color to
    each variant in the LBC.
    Parameter:
    - seed: int. Seed to randomly order variant
                 and produce different dictionaries.
    - colors: List. List of colors assigned to variant."""

    # # Shuffle gene order.
    # # Changing the seed will result in a different dictionary
    # np.random.seed(seed)
    # np.random.shuffle(variant)

    # Set color palette
    colors = px.colors.qualitative.D3  # create a list of plotly colors

    var_color_dict = dict()
    for i, gene in enumerate(variant):
        var_color_dict[gene] = colors[i % len(colors)]

    # Create path for exporting
    path = '../Resources/'
    if not os.path.exists(path):
        os.makedirs(path)
    # convert dictionary into string
    # using json.dumps()
    with open('../Resources/var_dict.json', 'w') as fp:
        json.dump(var_color_dict, fp)

