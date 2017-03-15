from itertools import chain
from toolz import sliding_window
from six.moves import zip_longest
from itertools import islice

def chunk_slices(length,size):
  """
  Given a length, and a chunk size, return an iterator over tuples containing start/end values

  For example: 

    Length - 20
    Chunksize - 4
  
    In [10]: for c in chunk_slices(20,4):
        ...:     print(c)
        ...:
    (0, 4)
    (4, 8)
    (8, 12)
    (12, 16)
    (16, 20)

  Each one can be used as a slice: 

  array[slice(*c)]

  """

  endpoints = range(0,length,size)
  return zip_longest(endpoints,islice(endpoints,1,None),fillvalue=length)

def chunk_word(word):
  """
  Returns an iterator over every possible ordered subset of a word. 
  
  For example:

  apple ->

    a
    p
    p
    l
    e
    ap
    pp
    pl
    le
    app
    ppl
    ple
    appl
    pple
    apple

  """

  return it.chain(*(tz.sliding_window(i,word) for i in xrange(1,len(word))))

