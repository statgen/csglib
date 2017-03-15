import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

def padjust_qvalue(pvalues):
  qvalue = importr("qvalue",on_conflict="warn")
  res = qvalue.qvalue(FloatVector(pvalues))
  return np.array(res.rx("qvalues")).flatten()

def padjust_bh(pvalues):
  import warnings
  warnings.warn("Use statsmodels.sandbox.stats.multicomp.multipletests instead",DeprecationWarning)

  stats = importr('stats',on_conflict="warn")

  adjusted = stats.p_adjust(FloatVector(pvalues),method='BH')

  series = pd.Series([x[1] for x in adjusted.items()])

  return series
