import pandas as pd
import concurrent.futures
import multiprocessing
from csg.iterhelp import chunk_slices

def hdf_get_nrows(hdf_file,key="df"):
  """
  Given an HDF file that was written by pandas, quickly determine the number
  of rows in the dataset
  """

  store = pd.HDFStore(hdf_file,mode="r")
  with store:
    n = store.get_storer(key).nrows
  
  return n

def hdf_head(hdf_file,key="df",n=10):
  """
  Given an HDF file that was written by pandas: 

  Retrieve the first n rows from the dataset
  """

  store = pd.HDFStore(hdf_file,mode="r")
  with store:
    it = store.select(key,chunksize=n,iterator=True)
    for df in it:
      break
    
    return df

def hdf_parallel_read(hdf_file,key="df",nworkers=None,chunksize=int(1E6),*args,**kwargs):
  """ 
  Read from a pandas HDF file in parallel

  Notes on chunksize: 
    * Changing it can paradoxically affect performance in both directions
    * Too large chunks make the master process take forever to pull data from worker processes,
      and worker processes do not terminate and move on to the next chunk until they transmit
      their data back to the master process
    * Too small chunks incur a fair amount of overhead
    * Due to an unfortunate "bug" in pandas, columns=[] does not restrict until AFTER load, so each process will 
      try to load the entire data frame (but only those rows within the chunk). If there are a lot of columns, 
      setting a large chunksize could exceed the memory on the machine
    * Eventually chunks will become so large as to exceed the ability of pickle protocol 3 to serialize them and send
      back to the master process. There is currently no way (without modifying multiprocessing.queue.SimpleQueue.put) to set
      the pickle protocol to version 4, which can handle much larger data

  Same general idea applies with nworkers. 

  The defaults here are set such that chunks are roughly 1 million rows, and the number of workers
  will try to default to 4, unless there aren't that many CPUs on the machine (and then it tries to use
  1 less than the max number of CPUs.) 
  """

  if nworkers is None:
    #      ncpus_used  ncpus_available
    #           1               -1
    #           1                0
    #           1                1
    #           1                2
    #           2                3
    #           3                4
    #           4                5
    #           4                6
    nworkers = min(4,max(1,multiprocessing.cpu_count()-1))

  nrows = hdf_get_nrows(hdf_file)

  with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = []
    for sl in chunk_slices(nrows,chunksize):
      future = executor.submit(pd.read_hdf,hdf_file,key=key,start=sl[0],stop=sl[1],mode="r",*args,**kwargs)
      futures.append(future)

  final = pd.concat([f.result() for f in futures])
  return final

def _parallel_apply_execfn(hdf_file,func,*args,**kwargs):
  """
  func should be a function that takes a data frame and returns one
  Returned data frames should have the same columns (not necessarily the same columns as the original)
  """

  df = pd.read_hdf(hdf_file,*args,**kwargs)
  result = func(df)
  return result

def hdf_parallel_apply(hdf_file,func,key="df",nworkers=None,chunksize=int(1E6),*args,**kwargs):
  if nworkers is None:
    nworkers = multiprocessing.cpu_count() - 2

  nrows = hdf_get_nrows(hdf_file)

  with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = []
    for sl in chunk_slices(nrows,chunksize):
      future = executor.submit(_parallel_apply_execfn,hdf_file,func,key=key,start=sl[0],stop=sl[1],mode="r")
      futures.append(future)

  final = pd.concat([f.result() for f in futures])
  return final

