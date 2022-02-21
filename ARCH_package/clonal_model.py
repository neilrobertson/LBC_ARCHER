# =============================================================================
# Created By  : Eric Latorre
# Created Date: 2022-01
# =============================================================================
import distributions

import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from numpy import concatenate

colors = px.colors.qualitative.D3

from scipy.optimize import minimize
from scipy.stats import beta
from scipy.stats import nbinom
from scipy.stats import binom
from scipy.stats import expon

from operator import attrgetter

# =============================================================================
#region Creating classes
# =============================================================================

class clonal_fitness_model():
    def __init__(self, clonal_map, mutation_map, data,
                 trajectories_map, time_points,
                 common_origin_map=None, fit_list=None, optimal_fit=None):
        self.clonal_map = clonal_map
        self.mutation_map = mutation_map
        self.data = data
        self.trajectories_map = trajectories_map
        self.time_points = time_points
        self.fit_list = fit_list
        self.optimal_fit = optimal_fit
        self.common_origin_map = common_origin_map


#endregion
# =============================================================================
#region Create cohort for clonal modelling
# =============================================================================

def cohort_create_clonal_models(cohort):
    for part in cohort:
        part.trajectories = [traj for traj in part.trajectories if traj.filter is True]
        part.mutation_list = [traj.mutation for traj in part.trajectories ]
    
    cohort = [part for part in cohort if len(part.trajectories)>0]


    # Create clonal models for each participant
    for part in cohort:
        traj_count = len(part.trajectories)
        traj_index_list = list(range(traj_count))

        clonal_partition_map_list = [partition_map for partition_map 
                                    in partition(traj_index_list)]

        # extract time_points
        time_points = set()
        for traj in part.trajectories:
            time_points = time_points.union(set(traj.data.age))
        time_points = list(time_points)
        time_points.sort()

        # Extract participant VAF data
        data = np.array([traj.data.AF for traj in part.trajectories])

        #extract mutation list in participant
        mutation_list = np.array(part.mutation_list)

        # Create clonal model for each partition map
        part.clonal_models = []
        for clonal_map in clonal_partition_map_list:
            # create mutation map for each clonal map
            mutation_map = []
            for index_map in clonal_map:
                mutation_map.append(
                    list(mutation_list[index_map]))

            trajectories_map = []
            for clone_mutations in mutation_map:
                trajectories_map.append(
                    [traj.data for traj in part.trajectories
                     if traj.mutation in clone_mutations]
                )

            part.clonal_models.append(
                clonal_fitness_model(
                    clonal_map=clonal_map,
                    mutation_map=mutation_map,
                    trajectories_map=trajectories_map,
                    data=data,
                    time_points=time_points))

    return cohort

def partition(collection):
    if len(collection) == 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset 
        yield [ [ first ] ] + smaller

#endregion
# =============================================================================
#region Deterministic clonal modelling approach
# =============================================================================

def part_deterministic_fit (part):
    """Plot the optimal deterministic fit using the optimal model
    and fitness as suggested by the Bayesian """

    model = part.optimal_model
    # Create a list of deterministic models 
    model_fit = []
    for  i in range(100):
        # Each mutation has their own origin. 
        # Upper bound depends on first non-zero observation

        # First non-zero VAF ocurrence for each trajectory
        index_first_nonzero = (model.data!=0).argmax(axis=1)

        # assuming 100_000 stem cells, a lower bound of
        # -4 makes translates into a maximum init vaf of 0.4
        # uppder bound is determined by first non-zero vaf point
        traj_origin_bounds = [(-100, model.time_points[index])
                        for index in index_first_nonzero]

        bounds = traj_origin_bounds

        # initialise parameters uniformly inside bounds
        x0 = [np.random.uniform(low=param[0], high=param[1])
            for param in bounds]

        # find optimal origins of deterministic trajectories
        fit = minimize(clonal_residual_optimal_fit, x0=x0,
                    bounds=bounds,
                    args=part,
                    method='Nelder-Mead')
        model_fit.append(fit)

    # find deterministic model wth minimum squared error.
    part.optimal_model.deterministic_fit = min(model_fit, key=attrgetter('fun'))

    # Predict deterministic trajectories using optimal origins
    time_points, vaf_points = clonal_model_deterministic_prediction(part)

    # add deterministic trajectories attribute
    for i, traj in enumerate(part.trajectories):
        traj.deterministic_fit = np.vstack((time_points, vaf_points[:,i]))
        
    # Create deterministic plot
    fitness  = part.optimal_model.fitness_map

    fig = go.Figure()
    for i, map in enumerate(part.optimal_model.clonal_map, 1):
        for index in map:
            fig.add_trace(
                go.Scatter(
                    x=part.trajectories[index].data.age,
                    y=part.trajectories[index].data.AF,
                    legendgroup=f'clone {i}',  # this can be any string, not just "group"
                    legendgrouptitle_text=f'Clone {i}',
                    name=part.trajectories[index].mutation,
                    marker_color=colors[index],
                    mode='markers'
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=time_points,
                    y=vaf_points[:,index],
                    legendgroup=f'clone {i}',  # this can be any string, not just "group"
                    marker_color=colors[index],
                    showlegend=False,
                    # hovertext=(f"Fitness: {parameters[i]:.2f} <br>"
                    #             f"Origin:{parameters[1+len(part.optimal_model.clonal_map)+index]}")
                )
            )

    fig.update_layout(title=part.id,
                    xaxis_title='Age',
                    yaxis_title='AF')

    part.deterministic_plot = fig


def clonal_residual_optimal_fit (params, part):
    N_w = 100_000
    model = part.optimal_model
    # extract clonal map
    clonal_map = model.clonal_map

    # Extract parameters
    traj_origin = np.array(params)

    clonal_fitness = model.fitness_map
    traj_fitness = np.array([traj.fitness for traj in part.trajectories])    

    clonal_origin = clonal_fitness.copy()
    for i, index_list in enumerate(clonal_map): 
        clonal_origin[i] = np.min(traj_origin[index_list])

    clonal_origin = np.array(clonal_origin)
    # Compute residual for each time point
    total_residual = 0

    for i, t in enumerate(model.time_points):
        clonal_cells = exponential_cell_growth_alt(t, clonal_fitness, clonal_origin)
        total_cells = N_w + np.sum(clonal_cells)

        traj_cells = exponential_cell_growth_alt(t,traj_fitness, traj_origin)
        traj_vafs = traj_cells/(2*total_cells) 

        observed_vafs = model.data[:, i]
        time_point_residual = np.linalg.norm(traj_vafs-observed_vafs)**2 
        total_residual +=  time_point_residual  

    return total_residual


def clonal_model_deterministic_prediction (part):

    model = part.optimal_model
    params = model.deterministic_fit.x
    N_w = 100_000

    # extract clonal map
    clonal_map = model.clonal_map

    # Extract parameters
    traj_origin = np.array(params)

    clonal_fitness = model.fitness_map
    traj_fitness = np.array([traj.fitness for traj in part.trajectories])    

    clonal_origin = clonal_fitness.copy()
    for i, index_list in enumerate(clonal_map): 
        clonal_origin[i] = np.min(traj_origin[index_list])

    clonal_origin = np.array(clonal_origin)
    for i, index_list in enumerate(clonal_map): 
        clonal_origin[i] = np.min(traj_origin[index_list])


    time_points = np.linspace(model.time_points[0]-2,
                            model.time_points[-1]+2,
                            100)
    vaf_points = []
    for t in time_points:
        clonal_cells = exponential_cell_growth_alt(t, clonal_fitness, clonal_origin)
        total_cells = N_w + np.sum(clonal_cells)

        traj_cells = exponential_cell_growth_alt(t,traj_fitness, traj_origin)
        traj_vafs = traj_cells/(2*total_cells) 
        vaf_points.append(traj_vafs)

    return time_points, np.array(vaf_points)

def exponential_cell_growth_alt (t, fitness, origin):
    origin_positive = np.where(origin<0, 0 ,origin)

    origin_negative = np.where(origin<0, np.abs(origin) ,0)
    origin_negative

    cells = np.exp(fitness*(t-origin_positive)) + origin_negative*100_000  
    return cells

#endregion
# =============================================================================
#region Bayesian clonal model
# =============================================================================
def clone_conditional_likelihood(clonal_model, clone_index, s, N_w=100_000,
                                 alpha=0.01, resolution=25):
    if clone_index > len(clonal_model.clonal_map):
        print(f"Clone index error. Index must be < {len(clonal_model.clonal_map)}.")

    # Allocate memory for clonal vaf
    clonal_vaf = np.empty((len(clonal_model.clonal_map),
                                clonal_model.data.shape[1]))

    # Determine clonal vaf using the maximum avg mutation 
    # in each clone.
    leading_mutation_map = []
    for i, index_map in enumerate(clonal_model.clonal_map):
        #compute average vaf for each mutation in clone
        avg_mutation_vaf = clonal_model.data[index_map].mean(axis=1)
        
        # locate maximum average vaf
        max_id = np.argmax(avg_mutation_vaf)

        # leading_mutation map
        leading_mutation_map.append(index_map[max_id])

        # extract maximum clonal vaf
        clonal_vaf[i] = clonal_model.data[index_map][max_id]

        
    # Compute total clonal stem cells in participant
    total_clonal_cells = (
        distributions.vaf_to_clone_size(clonal_vaf.sum(axis=0), N_w))

    # Compute total stem cells in participant
    total_cells = N_w + total_clonal_cells

    # Compute cells with each mutation
    mutation_cells = distributions.deterministic_vaf_to_clone_size(
        clonal_model.data, total_cells)

    # Extract trajectories in clone and append cell counts
    # to each trajectory
    clone = clonal_model.trajectories_map[clone_index]
    clone_mutation_cells = mutation_cells[clonal_model.clonal_map[clone_index]]

    for traj, cells in zip(clone, clone_mutation_cells):
        traj['cell_count'] = cells

    # compute likelihood for each trajectory in clone
    clone_likelihood = []
    for traj in clone:
    # Step 1: compute range of reasonable initial clone sizes based on data alone
    # Beta function
        beta_p_conf_int = beta.ppf(q=[alpha, 1-alpha],
                                a=traj.iloc[0].AO+1,
                                b=(traj.iloc[0].DP
                                    - traj.iloc[0].AO
                                    + 1))

        # List of binomial p ranges and clone_size_ranges
        beta_p_range = np.linspace(beta_p_conf_int[0],
                            min(beta_p_conf_int[1], 0.5),
                            resolution)

        # Transform p values to clone sizes
        init_range = distributions.vaf_to_clone_size(
                                beta_p_range,
                                total_cells[0]-traj.iloc[0].cell_count)

        # initialise clone size ranges of starting birth and death processes
        bd_init_ranges = [(init_range[0], init_range[-1])]


        # Compute first term of recursion
        recursive_prob = y_1_cond(init_range, traj.iloc[0], total_cells[0]-traj.iloc[0].cell_count)

        previous_x_range = init_range

        for i in range(1,len(traj)):
            
            next_x_range = find_range (traj, i, bd_init_ranges, s, total_cells[i]-traj.iloc[i].cell_count, alpha, resolution)

            recursive_prob = np.array([y_k_cond(x, traj.iloc[i], recursive_prob, previous_x_range,
                                                s, total_cells[i]-traj.iloc[i].cell_count) for x in next_x_range])
            previous_x_range = next_x_range

        probability = np.trapz(x=previous_x_range, y=recursive_prob)
        clone_likelihood.append(probability)

    # Compute the product of the likelihood of each trajectory
    total_likelihood = np.product(clone_likelihood)
    
    return total_likelihood


# def non_markov_clone_conditional_likelihood(clonal_model, clone_index, s, N_w=100_000,
#                                  alpha=0.01, resolution=25):
#     """clone_index -> use as index in trajectories_map"""

#     if clone_index > len(clonal_model.clonal_map):
#         print(f"Clone index error. Index must be < {len(clonal_model.clonal_map)}.")
#         return

#     # Allocate memory for clonal vaf
#     clonal_vaf = np.empty((len(clonal_model.clonal_map),
#                                clonal_model.data.shape[1]))

#     # Determine clonal vaf using the maximum avg mutation 
#     # in each clone.
#     leading_mutation_map = []
#     for i, index_map in enumerate(clonal_model.clonal_map):
#         #compute average vaf for each mutation in clone
#         avg_mutation_vaf = clonal_model.data[index_map].mean(axis=1)
        
#         # locate maximum average vaf
#         max_id = np.argmax(avg_mutation_vaf)

#         # leading_mutation map
#         leading_mutation_map.append(index_map[max_id])

#         # extract maximum clonal vaf
#         clonal_vaf[i] = clonal_model.data[index_map][max_id]

        
#     # Compute total clonal stem cells in participant
#     total_clonal_cells = (
#         distributions.vaf_to_clone_size(clonal_vaf.sum(axis=0), N_w))

#     # Compute total stem cells in participant
#     total_cells = N_w + total_clonal_cells

#     # Compute cells with each mutation
#     mutation_cells = distributions.deterministic_vaf_to_clone_size(
#         clonal_model.data, total_cells)
    
#     # Extract trajectories in clone and append cell counts
#     # to each trajectory
#     clone = clonal_model.trajectories_map[clone_index]
#     clone_mutation_cells = mutation_cells[clonal_model.clonal_map[clone_index]]

#     for traj, cells in zip(clone, clone_mutation_cells):
#         traj['cell_count'] = cells
    
#     # compute likelihood for each trajectory in clone
#     clone_likelihood = []
#     for trajectory in clone:
        
#         time_point_likelihood = []
#         for i in range(len(trajectory[:-1])):
#             init_time_point = trajectory.iloc[i]
#             next_time_point = trajectory.iloc[i+1]

#             # Beta function
#             beta_p_conf_int = beta.ppf(q=[alpha, 1-alpha],
#                                 a=init_time_point.AO+1,
#                                 b=(init_time_point.DP
#                                     - init_time_point.AO
#                                     + 1))

#             # List of binomial p ranges and clone_size_ranges
#             beta_p_range = np.linspace(beta_p_conf_int[0],
#                                     beta_p_conf_int[1],
#                                     resolution)

#             # We only want to integrate over beta_p_range < 0.5
#             beta_p_range = beta_p_range[beta_p_range < 0.5]

#             # Compute the pdf associated to each p value
#             beta_p_prob = beta.pdf(x=beta_p_range,
#                                 a=init_time_point.AO+1,
#                                 b=(init_time_point.DP
#                                     - init_time_point.AO
#                                     + 1))

#             # Transform p values to clone sizes
#             init_cs_range = distributions.vaf_to_clone_size(
#                                 beta_p_range,
#                                 total_cells[i]-init_time_point.cell_count)

#             # Initialise probability of observing AO alternate 
#             # observations at next time point
#             total_AO_prob = []

#             for init_cs in init_cs_range:
#                 # Compute next time point clone size range.
#                 # Notice that fitness is the same for trajectories
#                 # in the same clone.

#                 # Special case of initial clone size 0
#                 if init_cs == 0:
#                     if next_time_point.AO == 0:
#                         total_AO_prob.append(1)
#                     else:
#                         total_AO_prob.append(0)
                
#                     # continue to next initital clone size 
#                     continue

#                 # Extract proportion and trials NB approximating BD process
#                 nb_p, nb_n = distributions.BD_parametrization(
#                 init_size=init_cs, s=s, delta_t=init_time_point.delta_t)

#                 # Compute next clone size range interval
#                 next_cs_int = nbinom.ppf(q=[alpha, 1-alpha],
#                                         n=nb_n, p=nb_p)

#                 # Create next time point clone size range
#                 next_cs_range = np.linspace(next_cs_int[0],
#                                             next_cs_int[1],
#                                             resolution,
#                                             dtype='int')
#                 # Compute associated BD probabilities
#                 next_cs_prob = nbinom.pmf(k=next_cs_range, n=nb_n, p=nb_p)

#                 # Transform clone sizes to VAFs                
#                 next_binom_p_range = distributions.clone_size_to_vaf(
#                                         next_cs_range,
#                                         total_cells[i+1]-next_time_point.cell_count)
                                        
#                 next_cs_AO_prob = binom.pmf(k=next_time_point.AO,
#                                             n=next_time_point.DP,
#                                             p=next_binom_p_range)                
#                 total_AO_prob.append(
#                     np.trapz(x=next_cs_range, y=next_cs_prob*next_cs_AO_prob))

#             time_point_likelihood.append(
#                 np.trapz(x=beta_p_range, y=total_AO_prob*beta_p_prob))

#         clone_likelihood.append(np.product(time_point_likelihood))
    
#     total_likelihood = np.product(clone_likelihood)

#     return np.array(total_likelihood)


def create_uniform_prior(min_range, max_range, resolution, dtype='float'):
    range = np.linspace(min_range, max_range,
                        resolution, dtype=dtype)
    prob = [1 / (max_range-min_range)]*len(range)

    return np.vstack((range, prob))


def create_exponential_fitness_prior(min_fitness, max_fitness, resolution):
    """Approximation of the distribution of fitness
    reported in Blundell by an exponential function."""

    # Retrieve lambda from the equation
    # Exponential CDF valued at 0.04 (fitness) = 0.9
    lamb = -np.log(0.1)/0.04

    # Create array of fitness values
    fitness_range = np.linspace(min_fitness, max_fitness, resolution)
    # Retrieve associate pdf values
    fitness_pdf = expon.pdf(fitness_range, scale=1/lamb)
    fitness_prob = fitness_pdf/np.trapz(x=fitness_range, y=fitness_pdf)

    return np.vstack((fitness_range, fitness_prob))

flat_fitness_prior = create_uniform_prior(0, 1, 51)
flat_N_w_prior = create_uniform_prior(10_000, 200_000, 10)
exponential_fitness_prior = create_exponential_fitness_prior(min_fitness=0, max_fitness=1, resolution=21)


def clone_fitness_range_probability(clonal_model,
    fitness_prior=flat_fitness_prior, N_w=100_000):
    """Compute the posterior distribution of a clone conditional on a fixed number of stem cells"""

    # Create a list of posterior distributions
    fitness_range_probability = []
    
    # For each clone in clonal map compute the posterior distribution
    # of the fitness associated to that clone.
    for i in range(len(clonal_model.clonal_map)):
        fitness_posterior = [clone_conditional_likelihood(clonal_model, i, fitness, N_w)
                             for fitness in fitness_prior[0]]

        # Append a 2d array with the posterior fitness distribution of the clone
        fitness_range_probability.append(
            np.vstack((fitness_prior[0], fitness_posterior)))

    return fitness_range_probability


def clone_conditional_N_likelihood(clonal_model, N_w,
    fitness_prior=flat_fitness_prior):

    conditional_fitness_posterior = clone_fitness_range_probability(
        clonal_model, fitness_prior, N_w)

    # marginalise over fitness prior
    likelihood_fitness_marginal = []
    for fitness_posterior in conditional_fitness_posterior:
        likelihood_fitness_marginal.append(
            np.trapz(x=fitness_prior[0],
            y=fitness_posterior[1]*fitness_prior[1]))

    return likelihood_fitness_marginal


def clonal_model_probability(model,
                             fitness_prior=flat_fitness_prior,
                             N_w_prior=flat_N_w_prior,
                             return_model=False):

    likelihood_N_w_marginal = []
    for N_w in N_w_prior[0]:
        # For each N_w append the product over clones of
        # conditional probabilities.
        likelihood_N_w_marginal.append(
            np.product(
            clone_conditional_N_likelihood(model, N_w, fitness_prior)))

    likelihood_N_w_marginal = np.array(likelihood_N_w_marginal)

    model_probability = np.trapz(x=N_w_prior[0], y=likelihood_N_w_marginal*N_w_prior[1])
    
    model.model_probability = model_probability

    if return_model is True:
        return model


def bayesian_clonal_model_comparison (part,
                                      fitness_prior_clonal_structure=flat_fitness_prior,
                                      fitness_prior_plots = flat_fitness_prior,
                                      N_w_prior=flat_N_w_prior,
                                      model_evidence=0.5,
                                      deterministic_fit=True):

    if len(part.clonal_models) > 1:
        # Set a variable recording maximal model probability
        for model in part.clonal_models:
            # Compute model probability
            clonal_model_probability(model,
                                    fitness_prior=fitness_prior_clonal_structure,
                                    N_w_prior=N_w_prior)
        
        determine_optimal_model(part, model_evidence=model_evidence)

    else:
        part.optimal_model = part.clonal_models[0]
        
    part = optimal_model_post_processing(part,
                                  fitness_prior=fitness_prior_plots,
                                  N_w_prior=N_w_prior,
                                  deterministic_fit=deterministic_fit)

    return part


def determine_optimal_model(part, model_evidence=0.5, return_part=False):

    weighted_model_probability = np.array(
        [model.model_probability*(1 + model_evidence*(len(model.clonal_map)-1))
        for model in part.clonal_models])

    part.optimal_model = (
        part.clonal_models[np.nanargmax(weighted_model_probability)])
    
    if return_part is True:
        return_part
    

def optimal_model_post_processing(part,
                                  fitness_prior=flat_fitness_prior,
                                  N_w_prior=flat_N_w_prior,
                                  deterministic_fit=True):
        # Compute fitnes posterior distribution
    fitness_posterior_distribution(part.optimal_model,
                                   fitness_prior=fitness_prior,
                                   N_w_prior=N_w_prior)

    # Create conditional fitness distribution plot
    part.optimal_model.posterior_fitness_plot = (
        part_fitness_distribution(part.optimal_model, part.id))
    
    # Create clonal fitness map
    fitness_map = []
    for posterior in part.optimal_model.clonal_fitness_posterior:
        max_idx = np.argmax(posterior[1])
        max_likelihood_fitness = posterior[0, max_idx]
        fitness_map.append(max_likelihood_fitness)
    
    part.optimal_model.fitness_map = fitness_map
 
    # Create clonal fintess_quantiles map
    fitness_quantiles = []
    for fitness in part.optimal_model.clonal_fitness_posterior:
        posterior_random_sample = np.random.choice(fitness[0], 10_000, p=fitness[1]/fitness[1].sum())
        posterior_quantile = np.quantile(posterior_random_sample, [0.05, 0.95])
        fitness_quantiles.append(posterior_quantile)
    
    part.optimal_model.fitness_quantiles = fitness_quantiles


    # Set optimal fitness and quantiles to participant's trajectories
    for fitness, fitness_quantiles, mutation_list in zip(part.optimal_model.fitness_map,
                                  part.optimal_model.fitness_quantiles,
                                  part.optimal_model.mutation_map):
        for traj in part.trajectories:
            if traj.mutation in mutation_list:
                traj.fitness = fitness
                traj.fitness_quantiles = fitness_quantiles

    # Append posterior distribution to each trajectory
    for posterior,mutation_list in zip(part.optimal_model.clonal_fitness_posterior,part.optimal_model.mutation_map):
        for traj in part.trajectories:
            traj.gene=traj.mutation.split()[0]
            if traj.mutation in mutation_list:
                traj.fitness_posterior=posterior
                traj.clonal_population=len(mutation_list)

    if deterministic_fit is True:
        # Create deterministic plots conditional on optimal fitness
        part_deterministic_fit(part)

    return part

def clonal_fitness_map(model):
    """Extract fitness and quantiles from clone and create a clone map to mutations."""
    # Create clonal fitness map
    fitness_map = []
    for posterior in model.clonal_fitness_posterior:
        max_idx = np.argmax(posterior[1])
        max_likelihood_fitness = posterior[0, max_idx]
        fitness_map.append(max_likelihood_fitness)
    
    model.fitness_map = fitness_map
 
    # Create clonal fintess_quantiles map
    fitness_quantiles = []
    for fitness in model.clonal_fitness_posterior:
        posterior_random_sample = np.random.choice(fitness[0], 10_000, p=fitness[1]/fitness[1].sum())
        posterior_quantile = np.quantile(posterior_random_sample, [0.05, 0.95])
        fitness_quantiles.append(posterior_quantile)
    
    model.fitness_quantiles = fitness_quantiles


def marginalise_rest_of_clones(conditional_probability_fitness_range, fitness_prior):
    marginalised_clones_probability_product = []
    for i, clone_probability in enumerate(conditional_probability_fitness_range):
        marginalised_clones_probability = []
        for j, marginal_clone_probability in enumerate(conditional_probability_fitness_range):
            if j != i:
                marginalised_clones_probability.append(
                    np.trapz(x=fitness_prior[0],
                    y=marginal_clone_probability[1]*fitness_prior[1]))
        marginalised_clones_probability_product.append(
            np.product(marginalised_clones_probability))
    
    return marginalised_clones_probability_product

def fitness_posterior_distribution (model, fitness_prior, N_w_prior):
    """Compute the posterior distribution of fitness for every clone in model"""
    N_int = [] 

    # Initialise clonal contours fitness_N_w contour plots
    clonal_contour =  []
    for i in model.clonal_map:
        clonal_contour.append([])

    for N_w in N_w_prior[0]:
        # Compute conditional fitness_distributions for each clone
        conditional_probability_fitness_range = clone_fitness_range_probability(
            model, fitness_prior=fitness_prior, N_w=N_w)
        
        # Marginalised alternative clones
        alt_marg_clones = marginalise_rest_of_clones(conditional_probability_fitness_range,
                                                     fitness_prior)

        conditional_probability_fitness_range = np.array(conditional_probability_fitness_range).T
        conditional_probability_fitness_range [:,1] = conditional_probability_fitness_range [:,1]*alt_marg_clones

        for clone_idx, contour in enumerate(clonal_contour):   
            contour.append(list(conditional_probability_fitness_range[:,1, clone_idx]))
        # Compute conditional_N_s probability
        N_int.append(conditional_probability_fitness_range[:,1])

    
    contour_titles = tuple([', '.join(mutation_list) for mutation_list in model.mutation_map])
    fig = make_subplots(rows=len(clonal_contour), cols=1,
                        subplot_titles=contour_titles, shared_xaxes=True)
    
    for i, contour in enumerate(clonal_contour):
        fig.add_trace(
            go.Contour(
                z=contour,
                x=fitness_prior[0], # horizontal axis
                y=N_w_prior[0]),
            row=i+1, col=1)

    fig.update_xaxes(title='Fitness')
    fig.update_yaxes(title='Stem cells') 
        
    model.contour_plot = fig

    N_int = np.array(N_int)
    clonal_fitness_posterior = []
    for clone in range(N_int.shape[-1]):
        fitness_posterior = []
        for i in range(len(fitness_prior[1])):
            conditional_fitness_clone_marginal = (
                np.trapz(x=N_w_prior[0], y=N_w_prior[1]*N_int[:,i,clone])
            )
            fitness_posterior.append(conditional_fitness_clone_marginal)
        clonal_fitness_posterior.append(np.array(fitness_posterior))

    clonal_fitness_posterior = [np.vstack((fitness_prior[0], clonal_posterior)) for clonal_posterior in clonal_fitness_posterior]
    model.clonal_fitness_posterior = clonal_fitness_posterior


def part_fitness_distribution(model, part_id):
    """Creates a plot for the normalised posterior distribution of fitness,
    conditional on N_w = 100_000."""

    if hasattr(model, 'clonal_fitness_posterior') is not True:
        # Check that model has conditional_fitness_posterior attribute
        print('Error: Model not fitted')
        return

    fig = go.Figure()
    for posterior, mutation_list in (
        zip(model.clonal_fitness_posterior,
            model.mutation_map)):

        if len(mutation_list) == 1:
            name = mutation_list[0]
        else:
            name = '<br>'.join(mutation_list)

        fig.add_trace(
            go.Scatter(
                x=posterior[0],
                y=posterior[1] / np.max(posterior[1]),
                name=name,
                showlegend=True
            )
        )
    fig.update_layout(title=f"Posterior fitness distributions in {part_id}",
                      xaxis_title='Fitness',
                      yaxis_title="Normalised probability")
    return fig

#endregion
# =============================================================================
#region Hidden Markov
# =============================================================================
def y_1_cond(x_range, data, N_w):
    """Hidden Markov first recursion term"""
    prob_init = binom.pmf(k=data.AO,
                    n=data.DP,
                    p=distributions.clone_size_to_vaf(x_range, N_w=N_w))

    normalisation = np.trapz(x=x_range, y=prob_init)

    normalised_prob_init = prob_init/normalisation

    return normalised_prob_init

def y_k_cond(x, data, recursive_prob, previous_x_range, fitness, N_w):
    """Hiden Markov k-th recursion term"""
    data_observation_prob = binom.pmf(k=data.AO,
                                      n=data.DP,
                                      p=distributions.clone_size_to_vaf(x, N_w=N_w))

    # Compute the transition probability of the hidden process for all initial clone sizes
    # Extract proportion and trials NB approximating BD process
    nb_p, nb_n = distributions.BD_parametrization(
        init_size=previous_x_range, s=fitness, delta_t=3)
        
    # Compute associated BD probabilities
    next_cs_prob = nbinom.pmf(k=x, n=nb_n, p=nb_p)
    marginal = np.trapz(x=previous_x_range, y=recursive_prob*next_cs_prob)

    return data_observation_prob*marginal

def find_range (traj, k, bd_init_ranges, fitness, N_w, alpha = 0.01, resolution = 25):

    # BD RANGE
    # Find clone sizes with non-zero probability according to the hidden process
    # Extract proportion and trials NB approximating BD process
    hidden_range_total = [] 
    for i, cs_range in enumerate(bd_init_ranges):
        # compute time-steps since observation
        time_steps = k-i
        nb_p, nb_n = distributions.BD_parametrization(
            init_size=cs_range, s=fitness, delta_t=time_steps*3)

        # Compute minimum probable clone size
        min_clone_size = nbinom.ppf(q=alpha,
                                n=nb_n[0], p=nb_p[0])

        # Compute maximum probable clone size
        max_clone_size = nbinom.ppf(q=1-alpha,
                                n=nb_n[1], p=nb_p[1])

        # create a list of sampled clone sizes
        hidden_range = list(np.linspace(min_clone_size, max_clone_size, resolution, dtype='int'))
        hidden_range_total += hidden_range
    
    # Merge list of sampled clone sizes
    hidden_range_total = np.array(hidden_range_total)
    
    # OBSERVATION RANGE
    # Use Beta function to find range of clone sizes near observation
    beta_p_conf_int = beta.ppf(q=[alpha, 1-alpha],
                               a=traj.iloc[k].AO+1,
                               b=(traj.iloc[k].DP
                                  - traj.iloc[k].AO
                                  + 1))

    # List of binomial p ranges and clone_size_ranges
    beta_p_range = np.linspace(beta_p_conf_int[0],
                               min(beta_p_conf_int[1], 0.5),
                               resolution)

    # Transform p values to clone sizes
    observation_range = distributions.vaf_to_clone_size(v=beta_p_range,
                                                        N_w=N_w)
    
    bd_init_ranges.append((observation_range[0], observation_range[-1]))

    next_range = concatenate((observation_range, hidden_range_total))
    next_range.sort(kind='mergesort')

    return next_range     

# def hidden_markov_clone_conditional_likelihood(clonal_model, clone_index, fitness, N_w, alpha = 0.01, resolution=25):
#     if clone_index > len(clonal_model.clonal_map):
#         print(f"Clone index error. Index must be < {len(clonal_model.clonal_map)}.")

#     # Allocate memory for clonal vaf
#     clonal_vaf = np.empty((len(clonal_model.clonal_map),
#                                 clonal_model.data.shape[1]))

#     # Determine clonal vaf using the maximum avg mutation 
#     # in each clone.
#     leading_mutation_map = []
#     for i, index_map in enumerate(clonal_model.clonal_map):
#         #compute average vaf for each mutation in clone
#         avg_mutation_vaf = clonal_model.data[index_map].mean(axis=1)
        
#         # locate maximum average vaf
#         max_id = np.argmax(avg_mutation_vaf)

#         # leading_mutation map
#         leading_mutation_map.append(index_map[max_id])

#         # extract maximum clonal vaf
#         clonal_vaf[i] = clonal_model.data[index_map][max_id]

        
#     # Compute total clonal stem cells in participant
#     total_clonal_cells = (
#         distributions.vaf_to_clone_size(clonal_vaf.sum(axis=0), N_w))

#     # Compute total stem cells in participant
#     total_cells = N_w + total_clonal_cells

#     # Compute cells with each mutation
#     mutation_cells = distributions.deterministic_vaf_to_clone_size(
#         clonal_model.data, total_cells)

#     # Extract trajectories in clone and append cell counts
#     # to each trajectory
#     clone = clonal_model.trajectories_map[clone_index]
#     clone_mutation_cells = mutation_cells[clonal_model.clonal_map[clone_index]]

#     for traj, cells in zip(clone, clone_mutation_cells):
#         traj['cell_count'] = cells

#     # compute likelihood for each trajectory in clone
#     clone_likelihood = []
#     for traj in clone:
#     # Step 1: compute range of reasonable initial clone sizes based on data alone
#     # Beta function
#         beta_p_conf_int = beta.ppf(q=[alpha, 1-alpha],
#                                 a=traj.iloc[0].AO+1,
#                                 b=(traj.iloc[0].DP
#                                     - traj.iloc[0].AO
#                                     + 1))

#         # List of binomial p ranges and clone_size_ranges
#         beta_p_range = np.linspace(beta_p_conf_int[0],
#                             min(beta_p_conf_int[1], 0.5),
#                             resolution)

#         # Transform p values to clone sizes
#         init_range = distributions.vaf_to_clone_size(
#                                 beta_p_range,
#                                 total_cells[0]-traj.iloc[0].cell_count)

#         # Compute first term of recursion
#         recursive_prob = y_1_cond(init_range, traj.iloc[0], total_cells[0]-traj.iloc[0].cell_count)

#         previous_x_range = init_range

#         for i in range(1,len(traj)):
            
#             next_x_range = find_range (traj, i, init_range, fitness, N_w, alpha, resolution)

#             recursive_prob = np.array([y_k_cond(x, traj.iloc[i], recursive_prob, previous_x_range,
#                                                 fitness, N_w-traj.iloc[i].cell_count) for x in next_x_range])
#             previous_x_range = next_x_range

#         probability = np.trapz(x=previous_x_range, y=recursive_prob)
#         clone_likelihood.append(probability)

#         total_likelihood = np.product(clone_likelihood)
#         return total_likelihood
