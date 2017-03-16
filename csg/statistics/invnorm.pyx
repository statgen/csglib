import scipy.stats
import numpy as np
cimport numpy as np
cimport cython
from cpython cimport bool

ctypedef fused long_or_double:
    cython.long
    cython.double

cpdef np.ndarray[long_or_double] rank_ties_random(np.ndarray[long_or_double] vec):
  """
  Sorts a vector, but breaks ties at random. 
  
  This is mainly only applicable for vectors with many ties
  that you may want to ultimately inverse normalize. 
  
  NaN values are ranked last. 

  This function only operates on numpy arrays of type long or double (np.int64 or np.float64)
  """
  
  # Sort vector
  cdef np.ndarray[long_or_double] s = np.sort(vec)
  
  # Find indices needed to "unsort" the rank vector
  # Equivalent to vec.argsort().argsort() but faster
  # on long vectors
  cdef np.ndarray[long] rev = np.empty(len(vec),dtype=np.int)
  rev[np.argsort(vec)] = np.arange(len(vec))
  
  # Rank vector (built below)
  cdef np.ndarray[long] ranks = np.zeros(shape=len(vec),dtype=np.int)
  
  cdef bool run = False
  cdef int ri = 0
  cdef float right = np.nan
  cdef np.ndarray[long] tr
  cdef size_t i
  
  for i in range(1,len(vec)+1):  
    right = np.nan if i >= len(vec) else s[i]
    
    if right == s[i-1]:
      if not run:
        run = True
        ri = i-1
      
    elif right != s[i-1]:
      if run:
        run = False
        
        # Ties broken at random
        tr = np.arange(ri+1,i+1,dtype=np.int)
        np.random.shuffle(tr)
        ranks[ri:i] = tr
        
      else:
        ranks[i-1] = i
        
  # Unsort ranks back to original array order
  return ranks[rev]

def inverse_normalize(np.ndarray[long_or_double] vec):
  """
  Inverse normalize an array of values. 

  The array should be a numpy array of type long or double (np.int64 or np.float64)
  """

  nonmiss = (~np.isnan(vec)).sum()
  qs = (rank_ties_random(vec)-0.5)/nonmiss
  normed = scipy.stats.distributions.norm.isf(1-qs)
  normed[np.isnan(vec)] = np.nan
  return normed

