'''
stat_helpers.py
wrapper functions for statistics
'''
import numpy as np
from scipy.stats import fisher_exact
import math

def calc_fisher_exact(l1, l2, size):
    '''
    Calculate Fisher's exact result for getting this enrichment or greater in a group.
    l1 = list1 genes
    l2 = list2 genes
    size = number of genes in population
    See DAVID documentation for example:
    https://david.ncifcrf.gov/helps/functional_annotation.html
    The contingency matrix should be: 
    [x, n-x] = hits, pop_hits - hits
    [N-x, M-N-n+x] = sample-hits, genes in population and not in other groups
    Sum of matrix => all genes
    x = successes in sample
    M = size of population
    n = successes in population
    N = size of sample
    '''
    M = size
    N = len(l1)
    n = len(l2)
    x = len(set(l1).intersection(set(l2)))

    table = [[x, n-x], [N-x, M-N-n+x]]

    odds_r, p = fisher_exact(table, alternative='greater')
    try:
        log_p = -math.log(p, 10)
    except ValueError:
        log_p = np.nan

    a = table[0][0]
    b = table[0][1]
    c = table[1][0]
    d = table[1][1]
    try:
        # https://select-statistics.co.uk/calculators/confidence-interval-calculator-odds-ratio/
        oddr_recalc = (a*d)/(b*c)
        lower = np.exp(np.log(oddr_recalc) - 1.96*np.sqrt(np.sum([1/a, 1/b, 1/c, 1/d])))
        upper = np.exp(np.log(oddr_recalc) + 1.96*np.sqrt(np.sum([1/a, 1/b, 1/c, 1/d])))
    except ZeroDivisionError:
        odds_r = np.nan
        lower = np.nan
        upper = np.nan
        log_p = np.nan
    return odds_r, log_p, lower, upper, table