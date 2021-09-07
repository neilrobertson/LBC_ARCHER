# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================

import dill
import os


def dill_export(instance, path, file_name):
    """Export instance using dill.
    Parameters:
    - object: Serializable instance.
    """
    # Create path for exporting

    if not os.path.exists(path):
        os.makedirs(path)
    # convert dictionary into string
    # using json.dumps()
    with open(path + file_name + '.dill', 'wb') as fp:
        dill.dump(instance, fp)
