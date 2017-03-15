import contextlib

@contextlib.contextmanager
def timer(s):
  """
  Example: 

  with timer("new algorithm"):
    run_algorithm()

  Prints: 

  [new algorithm]: (time)s
  """

  import time

  start = time.time()
  yield
  end = time.time()

  print("[{}] {:.03g}s".format(s,end-start))

