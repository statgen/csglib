import cython
cimport numpy as np
import numpy as np

@cython.boundscheck(False)
def calleles(np.ndarray x, np.ndarray y):
  """
  Given two arrays, each of which contains an allele (e.g. an effect allele column, and non-effect allele column),
  return an array of the sorted alleles from each row. 
  
  For example: 
  
    allele1 allele2 calleles
          T       A      A,T
          G       A      A,G
          A       C      A,C
  
  This is mainly useful to test if sets of alleles are identical in two datasets. 
  
  It is cythonized to be used on very large datasets (eQTLs.)
  """
  
  cdef int i
  cdef np.ndarray result = np.empty(len(x),dtype=np.object)
  result[:] = ".,."
  
  for i in range(x.shape[0]):
    a = ",".join(sorted([x[i],y[i]]))
    result[i] = a
    
  return result

