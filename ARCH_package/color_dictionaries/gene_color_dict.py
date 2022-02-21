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


def create_dict(genes, seed=2, colors=px.colors.qualitative.D3):
    """ Create and export in JSON format a dictionary assigning a color to
    each gene in the LBC.
    Parameter:
    - seed: int. Seed to randomly order genes
                 and produce different dictionaries.
    - colors: List. List of colors assigned to genes."""

    # Shuffle gene order.
    # Changing the seed will result in a different dictionary
    np.random.seed(seed)
    np.random.shuffle(genes)

    # Set color palette
    colors = px.colors.qualitative.D3  # create a list of plotly colors

    gene_color_dict = dict()
    for i, gene in enumerate(genes):
        gene_color_dict[gene] = colors[i % len(colors)]

    # Create path for exporting
    path = '../Resources/'
    if not os.path.exists(path):
        os.makedirs(path)
    # convert dictionary into string
    # using json.dumps()
    with open('../Resources/gene_color_dict.json', 'w') as fp:
        json.dump(gene_color_dict, fp)
