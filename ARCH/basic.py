# Last modification: E. Latorre -- 25-10-2020

# Module to load a dataset into a comprehensive collection of participants.
# We define the participant class and define a function to process a cohort.

import pandas as pd
import numpy as np

import plotly.graph_objects as go


class participant:
    '''Definition of the participant class object.
    Attributes:
    - id: a string with the participant's id
    - mutation_list: a list of strings of all mutations detected at any age
    - trajectories: a list of trajectory class objects corresponding to each
                    trajectory present in a participant.
    - data: pandas dataframe slice of the full dataset containing information
            about this participant.
    Methods:
    - plot_trajectories: plots all trajectories found in a participant
    '''

    def __init__(self, id=None, mutation_list=None, trajectories=None,
                 data=None):
        self.id = id
        self.mutation_list = mutation_list
        self.trajectories = trajectories
        self.data = data

    def profile(self, germline=True):
        fig = go.Figure()

        # If germline is True then plot all all trajectories
        if germline is True:
            for traj in self.trajectories:
                fig.add_trace(go.Scatter(x=traj.data.age, y=traj.data.AF,
                                         mode='lines+markers',
                                         name=traj.mutation))
        # If Germline is False only plot trajectories that are not germline
        else:
            for traj in self.trajectories:
                if traj.germline is False:
                    fig.add_trace(go.Scatter(x=traj.data.age, y=traj.data.AF,
                                             mode='lines+markers',
                                             name=traj.mutation))
        fig.update_layout(title=f'Trajectories of participant {self.id}',
                          xaxis_title='Age (in years)',
                          yaxis_title='VAF')
        return fig


class trajectory:
    '''Definition of the trajectory class object.
    Attributes:
    - mutation: string with the name of the mutation this trajectory follows.
    = p_key: amino acid change.
    - germline: boolean with the germline status of the mutation.
    - data: pandas dataframe with the trajectory data.
    - gradient: gradient between the first
                and the last point of the trajectory.
    '''

    def __init__(self, mutation=None, p_key=None, data=None,
                 variant_class=None, germline=False, gradient=None):
        self.mutation = mutation
        self.p_key = p_key
        self.variant_class = variant_class
        self.germline = germline
        self.data = data
        self.gradient = gradient


def load(df):
    """ Transform a dataset into a list of participant class objects.
    Returns list of participants."""

    # Create column with amino acid changes key.
    # Step 1: Create 3 to 1 letter dictionary.
    d = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
         'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
         'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
         'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M',
         'Ter': '*'}

    # Step 2: append new columns containing gene + amino acid change
    df['p_key_1'] = (df['PreferredSymbol']
                     + ' ' + df['protein_substitution'].replace(d, regex=True))
    df['p_key'] = df['PreferredSymbol'] + ' ' + df['protein_substitution']

    # transform the time into years from 1st wave.
    df['age'] = df['wave']
    df.age = 3*(df.age-1)

    # initialize complete list of participants
    cohort = []

    for part in df.participant_id.unique():
        # create a list of participant objects with id ad data attributes.
        cohort.append(participant(id=part,
                                  data=df[df['participant_id'] == part]))

    # create the list of trajectories for each participant.
    for part in cohort:
        # Initialize the lists of mutations and trajectories
        traj = []
        mutation_list = []

        # add a trajectory for each mutation present in a participant
        for key in part.data.key.unique():
            # Create a copy of the data slice corresponding to a trajectory
            data = part.data[part.data['key'] == key].copy()
            if len(data) < 2:   # avoid trajectories withh only 1 point
                continue
            # Extract the 1_letter amino acid change code
            p_key = data['p_key_1'].unique()[0]
            variant = data['Variant_Classification'].unique()[0]

            # Compute participant's age
            if 'LBC0' in part.id:
                data.age = data.age + 79
            else:
                data.age = data.age + 70

            mutation_list.append(key)

            # germline condition
            germline = False
            if np.mean(data.AF) > 0.45:
                germline = True

            # compute overall gradient and update total_grad
            gradient = (np.diff(data.AF.iloc[[0, -1]])
                        / np.sqrt(np.diff(data.age.iloc[[0, -1]])))[0]

            # Compute the regularized gradient: Diff(vaf)/sqrt(Diff(t))
            data['regularized_gradient'] = np.append(
                                          np.diff(data.AF)
                                          / np.sqrt(np.diff(data.age)),
                                          None)

            data.regularized_gradient = data.regularized_gradient.astype(float)

            # append trajectory to the list of trajectories
            traj.append(trajectory(mutation=key,
                                   p_key=p_key,
                                   variant_class=variant,
                                   data=data[['AF', 'age',
                                              'regularized_gradient']],
                                   germline=germline,
                                   gradient=gradient
                                   ))

        # Add all trajectories to the participant
        part.trajectories = traj

        # add the list of all mutations to each individual
        part.mutation_list = mutation_list

    return cohort


def load_id(id, df):
    # process a participant based on its id
    part = participant(id=id, data=df[df['participant_id'] == id])

    # Initialize the lists of mutations and trajectories
    traj = []
    mutation_list = []

    # add a trajectory for each mutation present in a participant
    for key in part.data.key.unique():

        data = part.data[part.data['key'] == key].copy()
        if len(data) < 2:   # avoid trajectories withh only 1 point
            continue
        # Extract the 1_letter amino acid change code
        p_key = data['p_key_1'].unique()[0]
        # transform the time into years since first age
        data['age'] = data['wave']
        data.age = 3*(data.age-1)
        # Correct time for participants age
        if 'LBC0' in part.id:
            data.age = data.age + 79
        else:
            data.age = data.age + 70

        mutation_list.append(key)

        # germline condition
        germline = False
        if np.mean(data.AF) > 0.45:
            germline = True

        # compute overall gradient and update total_grad
        gradient = (np.diff(data.AF.iloc[[0, -1]])
                    / np.sqrt(np.diff(data.age.iloc[[0, -1]])))[0]

        gradient = gradient.astype(float)

        # Compute the regularized gradient: Diff(vaf)/sqrt(Diff(t))
        data['regularized_gradient'] = np.append(
                                      np.diff(data.AF)
                                      / np.sqrt(np.diff(data.age)),
                                      None)
        data.regularized_gradient = data.regularized_gradient.astype(float)

        # append trajectory to the list of trajectories
        traj.append(trajectory(mutation=key,
                               data=data[['AF', 'age',
                                          'regularized_gradient']],
                               p_key=p_key,
                               germline=germline,
                               gradient=gradient))

    # Add all trajectories to the participant
    part.trajectories = traj

    # add the list of all mutations to each individual
    part.mutation_list = mutation_list
    return part


def melt(cohort, filter=1, mutation=None):
    # melts a cohort into one pandas dataframe
    # similar to the original pd.dataframe but with relative gradients computed

    full = pd.DataFrame()

    for part in cohort:
        for traj in part.trajectories:
            if traj.data.AF.mean() < filter:
                new = traj.data
                new['mutation'] = traj.mutation
                new['germline'] = traj.germline
                new['fitness'] = new['regularized_gradient']/new['AF']
                full = full.append(new, ignore_index=True)

    # subset rows containing mutations in a particular gene
    if mutation is not None:
        full = full[full['mutation'].str.contains(mutation)]

    return full


def replace(df):
    '''Append column to dataset with
    1 letter amino acid change description'''
    singleletter = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q',
                    'Lys': 'K', 'Trp': 'W', 'Asn': 'N', 'Pro': 'P',
                    'Thr': 'T', 'Phe': 'F', 'Ala': 'A', 'Gly': 'G',
                    'Ile': 'I', 'Leu': 'L', 'His': 'H', 'Arg': 'R',
                    'Met': 'M', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y',
                    'Ter': '*'}

    df['p_key_1'] = df["protein_substitution"].replace(singleletter, regex=True).str.split('.', expand=True)[1]
    df['p_key_1'] = df['PreferredSymbol'] + ' ' + df['p_key_1']
    return df