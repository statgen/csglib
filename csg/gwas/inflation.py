from collections import namedtuple
import numpy as np
import scipy.stats

AssociationStats = namedtuple('AssociationStats', ['pvalue', 'effect', 'se'])

def inflation_from_effect_se(association_stats):
   pvalues = np.array([2 * scipy.stats.norm.cdf(-np.abs(x.effect / x.se)) for x in association_stats])
   return np.median(scipy.stats.chi2.ppf(1 - pvalues, 1)) / scipy.stats.chi2.ppf(0.5, 1)

def inflation_from_pvalue(association_stats):
   pvalues = np.array([x.pvalue for x in association_stats])
   return np.median(scipy.stats.chi2.ppf(1 - pvalues, 1)) / scipy.stats.chi2.ppf(0.5, 1)
