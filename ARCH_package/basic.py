# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2021-09
# =============================================================================
"""Module to load a dataset into a comprehensive collection of participants.
We define the participant class and define a function to process a cohort.
"""
# =============================================================================
# Imports
# =============================================================================

import json
import numpy as np
from tqdm import tqdm
# =============================================================================
# Import local packages
# =============================================================================

import plot
import auxiliary
import color_dictionaries.gene_color_dict as gene_color_dict
import color_dictionaries.variant_color_dict as variant_color_dict

# =============================================================================
# Create local classes
# =============================================================================


class participant:
    '''Participant class object.
    Attributes:
    - id: String. Participant's id.
    - cohort: str. LBC cohort. 
    - mutation_list: List of strings. Variants detected in participant.
    - trajectories: List of trajectory class objects.
                    Each element of the list corresponds to the trajectory
                    of a variant.
    - data: Pandas dataframe. Slice of the full dataset containing all
            information regarding this participant.

    Methods:
    - profile: Plots all variant trajectories
                         found in this participant.
    '''

    def __init__(self, id=None, cohort=None, mutation_list=None, trajectories=None,
                 data=None):
        self.id = id
        self.cohort=cohort
        self.data = data
        self.mutation_list = mutation_list
        self.trajectories = trajectories

        # Set default mutable attributes
        if trajectories is None:
            trajectories = []
            self.trajectories = trajectories
        if mutation_list is None:
            mutation_list = []
            self.mutation_list = mutation_list

    def profile(self, germline=True):
        """Plot all variant trajectories in a participant class object.
        Parameters:
        - germilne: Boolean. Include trajectories with
                             attribute germline = True.
        Return: Plotly figure.
        """

        profile_fig = plot.profile(self, germline)
        return profile_fig


class trajectory:
    '''Trajectory class object follows the evolution of a genetic trajectory.

    Attributes:
    - mutation: String. Variant followed in this trajectory.
    - p_key: String. Corresponding amino acid change.
    - germline: Boolean. Germline status of the mutation.
    - data: Pandas dataframe. Trajectory data.
    - gradient: float. VAF gradient between first last timepoint
                of this trajectory.
    '''

    def __init__(self, part_id=None, mutation=None, p_key=None, data=None,
                 variant_class=None, germline=False, gradient=None,
                 AF_outlier_pvalue=None, mdaf=None):

        self.part_id = part_id
        self.mutation = mutation
        self.p_key = p_key
        self.variant_class = variant_class
        self.germline = germline
        self.data = data
        self.gradient = gradient
        self.AF_outlier_pvalue = AF_outlier_pvalue
        self.mdaf = mdaf

# =============================================================================
# Cohort loading functions
# =============================================================================


def load(df, export_name=False, create_dict=False):
    """ Transform the full dataset into a list of participant class objects.
    Parameters:
    - df: pandas DataFrame. pandas dataframe in ARCHER_DX format.
    Returns:
    - cohort: list of participants class objects.
    """

    # Basic manipulations

    # Drop missing data points (DP=0)
    df = df.drop(df[df.DP == 0].index)

    # Create a column to store time step information
    df['delta_t'] = 0

    # Set cohort name
    df['cohort'] = 'LBC21'
    df.loc[df.participant_id.str.contains('LBC36') , 'cohort'] = "LBC36"

    # Dropped 1st time point of LBC21 Cohort
    df = df.drop(df[(df.cohort == "LBC21")
                    & (df.wave == 1)].index)

    # Drop participants from study
    # Load excluded participants list
    with open('../Resources/excluded_samples.json') as json_file:
        excluded_samples = json.load(json_file)

    df = df[df.participant_id.isin(excluded_samples) == False]

    # Create 1-letter amino acid change column
    # Load amino acid dictionary from resources
    with open('../Resources/amino_dict.json') as json_file:
        amino_dict = json.load(json_file)

    # Append new columns containing gene + amino acid change
    df['p_key_1'] = (df['PreferredSymbol']
                        + ' ' + df['protein_substitution'].replace(amino_dict,
                                                                regex=True))
    df['p_key'] = df['PreferredSymbol'] + ' ' + df['protein_substitution']

    # Transform the wave into age in years
    # transform wave into years
    df['age'] = df['wave']
    df.age = 3*(df.age-1)

    # add age of participant depending on cohort
    df.loc[df.cohort == "LBC21", 'age'] = (
        df.loc[df.cohort == "LBC21", 'age'] + 79)
    df.loc[df.cohort == "LBC36", 'age'] = (
        df.loc[df.cohort == "LBC36", 'age'] + 70)

    # Create participants and genetic trajectories
    # initialize complete list of participants
    cohort = []
    id_list = list(df.participant_id.unique())
    for part_id in tqdm(id_list):
        part_df = df[df.participant_id == part_id]
        # create a list of participant objects with id ad data attributes.
        new_part = participant(id=part_id,
                               cohort=part_df.cohort.unique()[0],
                               data=part_df)

        # Create trajectories for each key in the participant
        for key in new_part.data.key.unique():
            load_key(new_part, key)

        # Append new participant to cohort
        cohort.append(new_part)

    if create_dict is True:
        # Create dictionary assigning colors to each gene present in the cohort
        gene_color_dict.create_dict(df.PreferredSymbol.unique())
        variant_color_dict.create_dict(df.Variant_Classification.unique())
        

    if export_name is not False:
        # Export cohort.
        auxiliary.dill_export(cohort, path='../Exports/',
                              file_name=export_name)

    return cohort


def load_key(part, key):
    """Creates a trajectory class object associated to a slice of the DataFrame
    associated with a genetic mutation, or key. During this process, the
    participant's trajectory list and mutation_list is updated.
    Parameters:
    - key: string. Name of mutation.
    Returns:
    """
    # Create a copy of the data slice corresponding to a trajectory
    key_data = part.data.loc[part.data['key'] == key].copy()
    key_data = key_data.reset_index()

    # avoid trajectories withh only 1 point
    if len(key_data) < 2:
        return

    # Compute time steps between observations 
    # Notice that we append a 0 at last observation
    key_data['delta_t'] = np.append(
        np.diff(key_data.age), 0)
        
    # avoid trajectories without any observation (due to drop of 1st timepoint)
    if (key_data.AO == 0).all():
         return

    # append key to mutation_list
    part.mutation_list.append(key)
    # create a new trajectory
    key_trajectory = trajectory(
        part_id=part.id,
        mutation=key, data=key_data[['AF', 'age', 'DP', 'AO', 'delta_t']])

    # Set protein change attribute
    key_trajectory.p_key = key_data['p_key_1'].iloc[0]

    # Set Variant classification attribute
    key_trajectory.variant_class = key_data['Variant_Classification'].iloc[0]

    # Set outlier values
    if key_data.AF_Outlier_Pvalue.isnull().all().all() == False:
        key_trajectory.AF_outlier_pvalue = (
            np.nanmean(key_data.AF_Outlier_Pvalue))

    if key_data['95MDAF'].isnull().all().all() == False:
        key_trajectory.mdaf = np.nanmean(key_data['95MDAF'])

    # germline condition
    germline = False
    if np.mean(key_data.AF) > 0.45:
        germline = True
    key_trajectory.germline = germline

    # compute overall gradient and update total_grad
    key_trajectory.gradient = (np.diff(key_data.AF.iloc[[0, -1]])
                               / np.sqrt(np.diff(key_data.age.iloc[[0, -1]])))[0]

    part.trajectories.append(key_trajectory)