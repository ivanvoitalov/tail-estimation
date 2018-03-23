import sys
import numpy as np
from collections import defaultdict, Counter
from matplotlib import pyplot as plt

"""Script to estimate network degree distribution tail index.

This script uses several well-established techniques to give
estimate of a tail index for degree distribution of a given network.
Hill, moments, Pickands and kernel-type estimators are implemented
and their estimates are plotted on multiple panels along with
PDF/CCDF plot for degree distribution.
Network should be provided in a form of an edgelist.

Example:
    Usage example is here::

        $ python tail_estimation.py edgelist_file_path

"""

def get_ccdf(data_sequence):
    """Function to calculate complimentary CDF for a given dataset.

    Args:
        data_sequence (list): the list of data to be analyzed.

    Returns:
        x (list): the list of data in sorted order.
        y (list): CCDF values for corresponding entries in x.

    """
    x = sorted(data_sequence, reverse = True)
    N = len(x)
    y = 1. * np.arange(N) / (N - 1.)
    return x, y

def get_distribution(data_sequence, number_of_bins = 30):
    """Function to get log-binned PDF of a dataset.

    Args:
        data_sequence (list): the list of data to be analyzed.
        number_of_bins (int): number of logarithmic bins to use.

    Returns:
        x (list): mid points of logarithmic bins.
        y (list): PDF values of corresponding bin points from x.

    """
    # define the support of the distribution
    lower_bound = min(data_sequence)
    upper_bound = max(data_sequence)
    # define bin edges
    log = np.log10
    lower_bound = log(lower_bound) if lower_bound > 0 else -1
    upper_bound = log(upper_bound)
    bin_edges = np.logspace(lower_bound, upper_bound, number_of_bins)
    
    # compute the histogram using numpy
    y, _ = np.histogram(data_sequence, bins = bin_edges, density = True)
    # compute for each y the value of x
    x = bin_edges[1:] - np.diff(bin_edges) / 2.0
    # if bin is empty, drop it from the resulting list
    drop_indices = [i for i,k in enumerate(y) if k == 0.0]
    x = [k for i, k in enumerate(x) if i not in drop_indices]
    y = [k for i, k in enumerate(y) if i not in drop_indices]
    return x, y

def get_frequencies(data_sequence):
    """Function to compute appearance frequencies for values in a dataset.
    
    Args:
        data_sequence (list): the list of data to be analyzed.

    Returns:
        x (list): list of unique values appearing in the data.
        y (list): frequencies of corresponding data values from x.

    """
    N = float(len(data_sequence))
    data_counts = Counter(data_sequence)
    x, y = [], []
    for k, count in data_counts.iteritems():
        x.append(k)
        y.append(count / N)
    return x, y

def hill_estimator(data_sequence, r_smooth = 2):
    """Function implementing Hill tail index estimator.

    Calculates original Hill estimator, smoothed Hill estimator
    and adjusted Hill estimator for all order statistics of the dataset.
    Also calculates :math:`k^{*}` -- a threshold value for order
    statistic based on KS-distance test proposed in Danielsson et al (2016).
    This value is calculated based on the adjusted Hill estimator.

    Args:
        data_sequence (list): the list of data to be analyzed.
        r_smooth (int, optional): smoothing parameter for smooth Hill
            estimator, usually set to 2 or 3. Cannot exceed (n/k).

    Returns:
        results (list): list of tuples of the form :math:`(k,
            \hat{\gamma_{Hill}}, \hat{\gamma_{smooHill}}, 
            \hat{\gamma_{adjHill}}, KS_{k})`.

    """
    pass

def moments_estimator(data_sequence, n_bootstrap = 0.1, r_bootstrap = 1.0):
    """Function implementing moments tail index estimator.

    Calculates moment estimator for a given dataset and uses
    double-bootstrap approach proposed in Draisma et al. (2000)
    to estimate an optimal order statistic :math:`k^{*}`.

    Args:
        data_sequence (list): the list of data to be analyzed.
        n_bootstrap (float): fraction of the original data to be
            bootstrap-sampled, should be in (0,1) range.
        r_boostrap (float): parameter controlling number of
            bootstrap repetitions, total number of repetitions is
            given by :math:`N_{rep.} = r N`, where :math:`N` is
            the number of data points in the original sequence.

    Returns:
        results (list): list of tuples of the form :math:`(k,
            \hat{\gamma_{moments}}).
        k_star (int): optimal order statistic calculated via
            double-bootstrap approach.

    """
    pass

def pickands_estimator(data_sequence):
    """Function implementing Pickands tail index estimator.

    Calculates Pickands estimator for a given dataset.

    Args:
        data_sequence (list): the list of data to be analyzed.

    Returns:
        results (list): list of tuples of the form :math:`(k,
            \hat{\gamma_{Pickands}}).

    """
    pass

def kernel_type_estimator(data_sequence, lower_h = 0.01, upper_h = 1.0):
    """Function implementing kernel-type tail estimator.

    Calculates kernel-type tail index estimator for a given dataset
    based on the biweight kernel proposed in Groeneboom et al. (2003)
    for a given order statistics bandwidth given by :math:`h_{l}, h_{u}`
    bounds.

    Args:
        data_sequence (list): the list of data to be analyzed.
        lower_h (float): lower bound for order statistics to consider
            given by :math:`k_{l}=\floor{h_{l}N}`, where :math:`N` is the
            length of the original data sequence.
        upper_h (float): upper bound for order statistics to consider
            given by :math:`k_{u}=\floor{h_{u}N}`, where :math:`N` is the
            length of the original data sequence.

    Returns:
        results (list): list of tuples of the form :math:`(k,
            \hat{\gamma_{kernel}}).

    """
    pass

def estimate_tail_index(edgelist_file, mode = 'undirected'):
    """Function to estimate tail index of degree distribution.

    Estimates tail index of degree distribution of a network provided
    in the form of an edgelist. Builds a collection of plots for
    implemented estimators and PDF/CCDF of degree distribution.

    Args:
        edgelist_file (str): path to the file containing network
            edgelist to be analyzed.
        mode (str): either undirected or directed, type of the
            network provided.

    Returns:
        None

    """
    pass

def main():
    pass

if __name__ == '__main__':
    main()