import numpy as np
import math
from scipy.stats import nbinom
from scipy.special import comb

# =============================================================================
# Clone size distributions
# =============================================================================


def exact_clone_size_prob(size: int, init_size: int,
                          delta_t: float, s: float) -> float:
    """Exact computation of P(x(t) = size | x(t_0)=init_size) for a
       birth and death process x(t).

    Parameters
    ----------
    size: int. Clone size where the probability is evaluated.
    init_size: int. Clone size during initial observation.
    delta_t: float. Interval of time between observations, t-t_0.
    s: float. Self-renewal advantage or fitness parameter.

    Returns
    ----------
    float. Probability of observing a given clone size at time t
           given init_size at time t_0.
    """
    # Self-replication rate / birth and death rates.
    lamb = 1.3

    # Case of fit clone
    if s > 0:
        # Compute alpha and beta values
        alpha_nom = lamb*(np.exp(s*delta_t)-1)
        beta_nom = (lamb+s)*(np.exp(s*delta_t)-1)
        denom = (lamb + s)*np.exp(s*delta_t)-lamb

        alpha = alpha_nom / denom
        beta = beta_nom / denom

        # For
        if size > 0:
            # Probability of observing a clone size
            total_sum = 0
            for j in range(min(size, init_size)+1):
                combinatorial_terms = (comb(init_size, j, exact=True)
                                       * comb(init_size + size - j - 1,
                                              init_size - 1,
                                              exact=True))

                binomial_sucess_alpha = np.power(alpha, init_size-j)
                binomial_sucess_beta = np.power(beta, size-j)
                binomial_failure = np.power(1-alpha-beta, j)

                total_sum += (combinatorial_terms
                              * binomial_sucess_alpha
                              * binomial_sucess_beta
                              * binomial_failure)
            return total_sum

        elif size == 0:
            return np.power(alpha, init_size)

        else:
            return print('Negative size input')
    if s == 0:
        # Notice that in the limit s-> 0 alpha = beta
        alpha = lamb*delta_t / (1+lamb*delta_t)

        if size > 0:
            total_sum = 0
            for j in range(min(size, init_size)+1):
                combinatorial_terms = (comb(init_size, j, exact=True)
                                       * comb(init_size + size - j - 1,
                                              init_size - 1, exact=True))
                binomial_success = np.power(alpha, init_size+size-2*j)
                binomial_failure = np.power(1-2*alpha, j)

                total_sum += (combinatorial_terms
                              * binomial_success
                              * binomial_failure)
            return total_sum
        elif size == 0:
            return np.power(alpha, init_size)
        else:
            return print('Negative size input')


def exact_clone_size_pmf(sizes: int, init_size: int,
                         delta_t: float, s: float) -> float:
    """Exact computation of P(x(t) = size | x(t_0)=init_size) for a birth and
    death process x(t).

    Parameters
    ----------
    size: Array. Clone size where the probabilities are evaluated.
    init_size: int. Clone size during initial observation.
    delta_t: float. Interval of time between observations, t-t_0.
    s: float. Self-renewal advantage or fitness parameter.

    Returns
    ----------
    float. Probability of observing a given clone size at time t given init_size at time t_0.
    """
    lamb = 1.3

    distribution = []

    # Case of fit clone
    if s > 0:
        # Compute alpha and beta values
        alpha_nom = lamb*(np.exp(s*delta_t)-1)
        beta_nom = (lamb+s)*(np.exp(s*delta_t)-1)
        denom = (lamb + s)*np.exp(s*delta_t)-lamb

        alpha = alpha_nom / denom
        beta = beta_nom / denom
        for x in sizes:
            if x > 0:
                # Probability of observing a clone size
                total_sum = 0
                for j in range(min(x, init_size)+1):
                    combinatorial_terms = (comb(init_size, j, exact=True)
                                           * comb(init_size + x - j - 1,
                                                  init_size - 1, exact=True))
                    binomial_sucess_alpha = np.power(alpha, init_size - j)
                    binomial_sucess_beta = np.power(beta, x - j)
                    binomial_failure = np.power(1 - alpha - beta, j)

                    total_sum += (combinatorial_terms
                                  * binomial_sucess_alpha
                                  * binomial_sucess_beta
                                  * binomial_failure)

                return distribution.append(total_sum)

            elif x == 0:
                distribution.append(np.power(alpha, init_size))

            else:
                return print('Negative size input')
    if s == 0:
        # Notice that in the limit s-> 0 alpha = beta
        alpha = lamb*delta_t / (1+lamb*delta_t)
        for x in sizes:
            if x > 0:
                total_sum = 0
                for j in range(min(x, init_size)+1):
                    combinatorial_terms = (comb(init_size, j, exact=True)
                                           * comb(init_size + x - j - 1,
                                                  init_size - 1,
                                                  exact=True))
                    binomial_success = np.power(alpha, init_size + x - 2*j)
                    binomial_failure = np.power(1 - 2*alpha, j)

                    total_sum += (combinatorial_terms
                                  * binomial_success
                                  * binomial_failure)
                distribution.append(total_sum)

            elif x == 0:
                distribution.append(np.power(alpha, init_size))

            else:
                return print('Negative size input')
    return distribution


def BD_stats(init_size: int, s: float, delta_t: float):
    """Extract the mean and variance of a birth and death process.

    Given an initial size, a time interval and a fitness, we can exatract the
    exact mean and variance of  P(x(t) = size | x(t_0)=init_size).

    Parameters
    ----------
    init_size: int. Clone size during initial observation.
    delta_t: float. Interval of time between observations, t-t_0.
    s: float. Self-renewal advantage or fitness parameter.

    Returns
    ----------
    mean, variance: tuple.
    """
    lamb = 1.3
    delta_t = np.array(delta_t)
    init_size = np.array(init_size, dtype=int)
    # make sure init_size is not 0
    init_size = np.where(init_size==0, 1, init_size)
    exp_term = np.exp(delta_t*s)

    mean = init_size*exp_term
    if s != 0:
        variance = init_size*(2*lamb + s)*exp_term*(exp_term-1)/s
    if s == 0:
        variance = 2*lamb*init_size*delta_t

    return mean, variance


def BD_parametrization(init_size: int, s: float, delta_t: float):
    """Extract the parameters p and n of a negative binomial
    to approximate a birth and death process.

    Given an initial size, a time interval and a fitness, we can exatract the
    exact mean and variance of  P(x(t) = size | x(t_0)=init_size).

    Parameters
    ----------
    init_size: int. Clone size during initial observation.
    delta_t: float. Interval of time between observations, t-t_0.
    s: float. Self-renewal advantage or fitness parameter.

    Returns
    ----------
    p, n: tuple.
    """
     # extract mean and variance of BD process
    mean, variance = BD_stats(init_size=init_size, s=s, delta_t=delta_t)
   

    # negative binomial probability and number of successes
    p = mean / variance
    n = np.power(mean, 2) / (variance-mean)  

    return p, n  



def nb_approx_pmf(size: int, init_size: int,
              s: float, delta_t: float) -> float:
    """
    Approximation of a birth and death process using the negative binomial.

    Parameters
    ----------
    size: int or Array. Clone size where the probability is evaluated.
    init_size: int. Clone size during initial observation.
    delta_t: float. Interval of time between observations, t-t_0.
    s: float. Self-renewal advantage or fitness parameter.

    Returns
    ----------
    prob: float. Negative binomial pmf with mean evaluated at size.
    """

    # extract mean and variance of BD process
    mean, variance = BD_stats(init_size=init_size, s=s, delta_t=delta_t)

    # negative binomial probability and number of successes
    p = mean / variance
    n = np.power(mean, 2) / (variance-mean)

    return nbinom.pmf(size, n, p)


def nb_approx_rvs(init_size: int,
              s: float, delta_t: float, size=1) -> float:
    """
    Generate expected clone sizes of a birth and death process using
    the negative binomial approxiamtion.

    Parameters
    ----------
    init_size: int. Clone size during initial observation.
    delta_t: float. Interval of time between observations, t-t_0.
    s: float. Self-renewal advantage or fitness parameter.
    size: int or Array. Number of generated samples.


    Returns
    ----------
    prob: float. Negative binomial pmf with mean evaluated at size.
    """

    # extract mean and variance of BD process
    mean, variance = BD_stats(init_size=init_size, s=s, delta_t=delta_t)

    # negative binomial probability and number of successes
    p = mean / variance
    n = np.power(mean, 2) / (variance-mean)

    return nbinom.rvs(n, p, size=size)


# =============================================================================
# VAF distributions
# =============================================================================
def vaf_nb_approx(vaf: int, init_vaf: int, N_w: int,
                  s: float, delta_t: float) -> float:
    """
    Approximation of the probability distribution of VAFs of a birth and death
    process using the negative binomial.

    Parameters
    ----------
    vaf: float or Array. VAF size where the probability is evaluated.
    init_vaf: float. VAF size during initial observation.
    N_w: int. Number of wild-type haematopoietic stem cells.
    delta_t: float. Interval of time between observations, t-t_0.
    s: float. Self-renewal advantage or fitness parameter.

    Returns
    ----------
    prob: float. Negative binomial pmf with mean evaluated at size.
    """
    size = vaf_to_clone_size(v=vaf, N_w=N_w)
    init_size = vaf_to_clone_size(v=init_vaf, N_w=N_w)

    # extract mean and variance of BD process
    mean, variance = BD_stats(init_size=init_size, s=s, delta_t=delta_t)

    # negative binomial probability and number of successes
    p = mean / variance
    n = np.power(mean, 2) / (variance-mean)

    return nbinom.pmf(size, n, p)


def vaf_gradient_nb_approx(gradient: float, init_vaf: float, N_w: int,
                           s: float, delta_t: float):
    """Probability of observing a given gradient between 2 VAF measurements
    steming from a birth and death process.

    Parameters
    ----------
    s: float. Self-renewal advantage or fitness parameter.
    N_w: int. Number of wild-type haematopoietic stem cells.
    gradient: float or Array. VAF gradient where the probability is evaluated.
    init_vaf: float. VAF size during initial observation.
    delta_t: float. Interval of time between observations, t-t_0.
    s: float. Self-renewal advantage or fitness parameter.

    Returns
    ----------
    prob: float. Negative binomial pmf with mean evaluated at size."""
    vaf = gradient_to_vaf(gradient, init_vaf, delta_t)

    return vaf_nb_approx(vaf=vaf, init_vaf=init_vaf,
                         N_w=N_w, s=s, delta_t=delta_t)


def gradient_to_vaf(gradient, init_vaf, delta_t):
    return gradient*delta_t + init_vaf


def vaf_to_clone_size(v: float, N_w: int):
    """Inverse function of v=x/(2N+2x) transforming VAF to clone sizes.

    Parameters
    ----------
    v: float. Variant allele frequency of a clone.
    N_w: int. Total number of wild type haematopoietic stem cells.

    Returns
    ----------
    int. Clone size corresponding to a given VAF.
    """
    clone_size = np.array(-N_w*v / (v-0.5), dtype='int')

    return clone_size


def clone_size_to_vaf(x: int, N_w: int):
    """VAF of a clone of size x given a total of N_w wild type stem cells.
    Function v(x) = x/(2x+2N)

    Parameters
    ----------
    x: int. Clone size.
    N_w: int. Total number of wild type haematopoietic stem cells.

    Returns
    ----------
    int. Clone size corresponding to a given VAF.
    """
    return x / (2*N_w+2*x)

def deterministic_vaf_to_clone_size(v, total_cells):

    clone_size = np.array(v*2*total_cells, dtype='int')
    return clone_size
    
def deterministic_clone_size_to_vaf(x, total_cells):
    return x/(2*total_cells)
