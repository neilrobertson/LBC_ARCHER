from ARCH import basic
from ARCH import simulation_functions as functions

import pickle
import pandas as pd
import multiprocessing as mp
import numpy as np
from functools import partial

from sklearn.linear_model import LinearRegression
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

pio.templates.default = "plotly_white"
colors = px.colors.qualitative.Plotly


# load fitted trajectories stored in a pickle file in 'model_list'
infile = open('pickles/cohort_neutral.pkl', 'rb')
cohort_neutral = pickle.load(infile)
infile.close()



# N and lambda from Dingli, Pacheco (Plosone, 2017)
# They estimate HSC divisions at a rate of 1 per year
# Using the variance of the gradients in the neutral cohort, we can estimate
# lambda ~ theta * N, yielding lamb ~ 0.25.
# This means that each HSC divides every 2 years symmetrically
N = cohort_neutral.HSC
lamb =  cohort_neutral.rate

N = 500
# Estimate mutation rate from synonymous mutations
df = pd.read_csv(r'Datasets/LBC_ARCHER.1PCT_VAF.Mar21.synonymous.tsv',
                 sep='\t')
syn = basic.load(df)

max_part = []
for part in syn:
    for traj in part.trajectories:
        max_part.append(traj.data['AF'].max())

max_syn = np.mean(max_part)
mu = functions.find_mutation_rate(lamb, N, np.mean(max_syn))

tmax = 80     # final time of simulation

timepoints = list(range(0, tmax+1))

bd = partial(functions.bd, lamb=lamb,
             mu=lamb, tmax=tmax, timepoints=timepoints)

# Number of simulated mutations
simulated_part = 1
n_mutations = simulated_part*mu*tmax
n_mutations  = 1000

# Create a list of clone immigration time
# immigration is uniformly distributed between age 0 and tmax
start = tmax*np.random.random_sample(int(n_mutations))

# compute the trajectory simulations in parallel for all mutations
pool = mp.Pool()
result = pool.map(bd, start)
result = list(result)

# As fitted in our data, we observe that the average fitness advantage lambda*s
# is about 0.2
# Append a clone with such a fitness advantage
fit_lamb = lamb + lamb*0.2
result.append(bd(50, lamb=fit_lamb))

# Create a filtered list of checkpoints dataframe with only non-empty
checkpoints = []
for traj in result:
    if traj[1].empty is False:
        checkpoints.append(traj[1])


checkpoint_counts = {i: N for i in timepoints}


# Compute the total counts at each timepoint
# checkpoint_counts = {70: 0, 73: 0, 76: 0, 79: 0, 82: 0, 85: 0, 88: 0}
for index, row in checkpoints[-1].iterrows():
    checkpoint_counts[row['time']] = row['count'] + N


# Compute the VAF for each time point in each trajectory
for traj in checkpoints:
    time = traj['time']
    count = traj['count']
    traj['VAF'] = [c / (2*(N + int(checkpoint_counts[t])))
                   for t, c in zip(time, count)]

VAF_fig = functions.plot_VAF(checkpoints[:-1])
VAF_fig.add_trace(
    go.Scatter(x=checkpoints[-1]['time'], y=checkpoints[-1]['VAF'],
               marker_color=colors[0],
               mode='lines',
               name='Fit mutation'))
# VAF_fig.update_yaxes(type='log',tickmode = 'linear',
#         dtick=1
# )
VAF_fig.write_image('Simulation/Fit simulation VAF log.svg')


time_range = range(70, 80, 3)
checkpoint_range = []
for traj in checkpoints:
    new = traj.loc[traj.time.isin(time_range)]
    if new.empty is False:
        checkpoint_range.append(new)

VAF_range = functions.plot_VAF(checkpoint_range)
# VAF_range.write_image('Simulation/Checkpoint_range.svg')

filtered_melt_syn, no_filter, filter, filter_bin = (
    functions.plot_vignette(checkpoint_range))


no_filter.write_image('Simulation/no_filter.svg')
filter.write_image('Simulation/filter.svg')
filter_bin.write_image('Simulation/filter_bin.svg')


#
# # Find the distribution of gradients
# melt_syn = pd.DataFrame(columns=['VAF', 'regularized_gradient'])
#
# # for all mutations but the fit (last one) compute local reg_gradients
# for traj in checkpoint_range[:-1]:
#     for combination in list(combinations(traj.index, 2)):
#         data = traj.loc[list(combination)]
#         gradient = np.diff(data.VAF) / np.sqrt(np.diff(data.time))
#         melt_syn = (
#             melt_syn.append({'VAF': traj.loc[combination[0]]['VAF'],
#                              'regularized_gradient': gradient[0]},
#                             ignore_index=True))
#
# melt_syn['bin'] = pd.cut(melt_syn['VAF'], 25)
#
# fig = go.Figure()
# fig.add_trace(
#     go.Scatter(x=melt_syn['VAF'], y=melt_syn['regularized_gradient'],
#                mode='markers', name='Trajectories'))
# fig.update_layout(title='Gradient distribution of neutral mutations',
#                    xaxis_title="VAF", yaxis_title='Gradient')
#
# fig = go.Figure()
# for bin in melt_syn.bin.unique():
#     data = melt_syn[melt_syn['bin']==bin]
#     fig.add_trace(
#         go.Scatter(x=data['VAF'], y=data['regularized_gradient'],
#                    mode='markers'))
# fig
# px.scatter(melt_syn, x='VAF', y='regularized_gradient')
#
# # Exclude all mutations with VAF < 0.01
# # lack of data below this threshold alters the distribution of gradients
# filtered_melt_syn = melt_syn[melt_syn['VAF'] > 0.01]

# Variance and mean with bootstrap for linear regression and plotting

bin_center = []     # track bin centers
mean_dist = []      # mean distribution in each bin using bootstrap
variance_dist = []  # variance distribution in each bin using bootstrap
bin_size = []       # Number of trajectories in each bin

for bin in filtered_melt_syn['bin'].unique():
    # filter dataset to obtain rows in bin
    VAF_bin = filtered_melt_syn[filtered_melt_syn['bin'] == bin]
    bin_size.append(len(VAF_bin))               # points in bin
    bin_center.append(VAF_bin['VAF'].mean())     # compute bin center

    # Bootstrap to compute gradient mean and variance distributions
    n = len(VAF_bin)
    reps = 200   # samples drawn
    subsets = np.random.choice(VAF_bin['regularized_gradient'], (n, reps))
    mean_dist.append(subsets.mean(axis=0))
    variance_dist.append(subsets.var(axis=0))

# compute mean and quantiles for each bin's mean
mean = [x.mean() for x in mean_dist]
full = [np.abs(np.quantile(x, [0.1, 0.9]) - np.mean(x)) for x in mean_dist]
mean_low = [x[0] for x in full]
mean_high = [x[1] for x in full]
# points in bin
# compute mean and quantiles for each bin's variance
variance = [x.mean() for x in variance_dist]
full = [np.abs(np.quantile(x, [0.1, 0.9]) - np.mean(x))
        for x in variance_dist]
var_low = [x[0] for x in full]
var_high = [x[1] for x in full]

# Weighted linear regression on mean and variance ~ VAF
for n, i in enumerate(bin_size):
    if i < 10:
        bin_size[n] = 0
size_weight = np.array(bin_size)


# Mean weighted model
mean_regr = LinearRegression(fit_intercept=True)
mean_regr.fit(np.c_[bin_center], mean, size_weight)
mean_score = mean_regr.score(np.c_[bin_center], mean, size_weight)

# Mean model figure
mean_fig = px.scatter(x=bin_center, y=mean, hover_data=[bin_size],
                      error_y=mean_high, error_y_minus=mean_low)
# trace mean linear fit
x = [0, max(bin_center)]
mean_fig.add_trace(
    go.Scatter(x=x, y=mean_regr.predict(np.c_[x]),
               name='linear regression fit'))
mean_fig.update_layout(title=('Linear regression gradient ~ VAF'))
mean_fig.update_xaxes(title='VAF')
mean_fig.update_yaxes(title='Gradient')

# Variance weighted model
var_regr = LinearRegression(fit_intercept=False)
var_regr.fit(np.c_[bin_center], variance, size_weight)
var_score = var_regr.score(np.c_[bin_center], variance, size_weight)

# Variance model figure
var_fig = px.scatter(x=bin_center, y=variance, hover_data=[bin_size],
                     error_y=var_high, error_y_minus=var_low)

x = [0, max(bin_center)]
var_fig.add_trace(
    go.Scatter(x=x, y=var_regr.predict(np.c_[x]),
               name='Linear regression fit'))
var_fig.update_layout(title=('Linear regression of variance ~ VAF'))
var_fig.update_xaxes(title='VAF')
var_fig.update_yaxes(title='Gradients Variance')

# Create a dataframe with the regularized gradient
# between first and last timepoint
total_gradient = pd.DataFrame(columns=['VAF', 'gradient'])
for traj in checkpoint_range:
    if len(traj) > 1:
        endpoints = traj.iloc[[0, -1]]
        gradient = np.diff(endpoints.VAF) / np.sqrt(np.diff(endpoints.time))
        total_gradient = (
            total_gradient.append({'VAF': traj.iloc[0]['VAF'],
                                   'gradient': gradient[0]},
                                  ignore_index=True))


# Compute the percentage of succesfully filtered trajectories:
n_std = 2
x = total_gradient['VAF']
total_gradient['filter'] = total_gradient['gradient'] \
    - (mean_regr.predict(np.c_[x])
       + n_std*np.sqrt(var_regr.predict(np.c_[x])))
percentage = 100 * len(total_gradient[total_gradient['filter'] < 0]) \
    / (len(total_gradient)-1)


# Plot filter
fig = go.Figure()

fig.add_trace(
    go.Scatter(x=total_gradient['VAF'],
               y=total_gradient['gradient'],
               name='Synonymous mutations',
               mode='markers'))
# Confidence line
x = np.linspace(0, max(total_gradient['VAF']), 1000)
y_std = mean_regr.predict(np.c_[x]) \
    + n_std*np.sqrt(var_regr.predict(np.c_[x]))
y_mean = mean_regr.predict(np.c_[x])

fig.add_trace(go.Scatter(x=x, y=y_mean,
                         name='Mean linear regression',
                         marker=dict(color='Lightseagreen')))

fig.add_trace(go.Scatter(x=x, y=y_std,
                         name='Neutral growth filter',
                         marker=dict(color='darkorange')))

fig.update_layout(title=(f'Neutral growth filter <br>Filters '
                         f'{round(percentage,1)}% of all '
                         f'synonymous mutations'))
fig.update_xaxes(title='VAF')
fig.update_yaxes(title='Regularized gradient')



mean_fig.write_image('Simulation/mean_fit.svg', scale=10)
var_fig.write_image('Simulation/var_fit.svg', scale=10)
fig.write_image('Simulation/trajectory_filter.svg', scale=10)
