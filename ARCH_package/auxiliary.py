# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================

import dill
import os
import numpy as np


def dill_export(instance, path, file_name):
    """Export instance using dill.
    Parameters:
    - object: Serializable instance.
    """
    # Create path for exporting
    if not os.path.exists(path):
        os.makedirs(path)

    # export instance using dill
    with open(path + file_name + '.dill', 'wb') as fp:
        dill.dump(instance, fp)


def prior_load(file_path):
    """Load a prior saved in numpy format.
    Parameters:
    -----------
    file_path: str. Path to file for loading.

    Returns:
    -----------
    prior. numpy or None. If path exists returns saved numpy array, else None."""
    try:
        with open(file_path, 'rb') as f:
            prior = np.load(f)

    except IOError:
        print(f"Could not read file: {file_path}")
        prior = None

    return prior