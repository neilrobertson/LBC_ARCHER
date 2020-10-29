# Last modification: E. Latorre -- 25-10-2020

# Module to load a dataset into a comprehensive collection of participants.
# We define the participant class and define a function to process a cohort.

import pandas as pd
import numpy as np
import multiprocessing as mp
import plotly.graph_objects as go


class participant:
    '''Definition of the participant class object.
    Attributes:
    - id: a string with the participant's id
    - mutation_list: a list of strings of all mutations detected at any wave
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
                fig.add_trace(go.Scatter(x=traj.data.wave, y=traj.data.AF,
                                         mode='lines+markers',
                                         name=traj.mutation))
        # If Germline is False only plot trajectories that are not germline
        else:
            for traj in self.trajectories:
                if traj.germline is False:
                    fig.add_trace(go.Scatter(x=traj.data.wave, y=traj.data.AF,
                                             mode='lines+markers',
                                             name=traj.mutation))
        fig.update_layout(title=f'Trajectories of participant {self.id}',
                          xaxis_title='Time (years since first wave)',
                          yaxis_title='VAF')
        return fig


class trajectory:
    '''Definition of the trajectory class object.
    Attributes:
    - mutation: string with the name of the mutation this trajectory follows
    - germline: boolean with the germline status of the mutation
    - data: pandas dataframe with the trajectory data
    - gradient: gradient between the first and the last point of the trajectory
    '''

    def __init__(self, mutation=None, data=None, germline=False,
                 gradient=None):
        self.mutation = mutation
        self.germline = germline
        self.data = data
        self.gradient = gradient


def load(df):
    # Transform a dataset into a list of participant class objects

    cohort = []               # initialize complete list of participants
    total_grad = []           # initialize complete list of gradients
    df.wave = 3*(df.wave-1)   # transform the time into years since first wave

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

            data = part.data[part.data['key'] == key].copy()
            if len(data) < 2:   # avoid trajectories withh only 1 point
                continue

            mutation_list.append(key)

            # germline condition
            germline = False
            if np.mean(data.AF) > 0.45:
                germline = True

            # compute overall gradient and update total_grad
            gradient = np.gradient(data.AF.iloc[[0, -1]].tolist(),
                                   data.wave.iloc[[0, -1]].tolist())[0]
            total_grad.append(gradient)

            # compute the relative gradients in data.gradient
            data['gradient'] = np.gradient(data.AF.tolist(),
                                           data.wave.tolist())

            # append trajectory to the list of trajectories
            traj.append(trajectory(mutation=key,
                                   data=data[['AF', 'wave', 'gradient']],
                                   germline=germline,
                                   gradient=gradient
                                   ))

        # Add all trajectories to the participant
        part.trajectories = traj

        # add the list of all mutations to each individual
        part.mutation_list = mutation_list

    return cohort, total_grad


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

        mutation_list.append(key)

        # germline condition
        germline = False
        if np.mean(data.AF) > 0.45:
            germline = True

        # compute overall gradient and update total_grad
        gradient = np.gradient(data.AF.iloc[[0, -1]].tolist(),
                               data.wave.iloc[[0, -1]].tolist())[0]

        # compute the relative gradients in data.gradient
        data['gradient'] = np.gradient(data.AF.tolist(),
                                       data.wave.tolist())

        # append trajectory to the list of trajectories
        traj.append(trajectory(mutation=key,
                    data=data[['AF', 'wave', 'gradient']],
                    germline=germline,
                    gradient=gradient))

    # Add all trajectories to the participant
    part.trajectories = traj

    # add the list of all mutations to each individual
    part.mutation_list = mutation_list
    return part


def melt(cohort, mutation=None):
    # melts a cohort into one pandas dataframe
    # similar to the original pd.dataframe but with relative gradients computed

    full = pd.DataFrame()

    for part in cohort:
        for traj in part.trajectories:
            new = traj.data
            new['mutation'] = traj.mutation
            new['germline'] = traj.germline
            new['fitness'] = new['gradient']/new['AF']
            full = full.append(new, ignore_index=True)

    # subset rows containing mutations in a particular gene
    if mutation is not None:
        full = full[full['mutation'].str.contains(mutation)]

    return full
