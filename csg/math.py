from numba import vectorize
import math

@vectorize("float64(float64)",nopython=True,target="parallel")
def fast_log10(x):
  return -math.log10(x)

def convert_to_log10(value):
  """
  Take the log10 of a **string** representing either a plain number,
  or number in scientific notation. This can handle arbitrary precision.

  This is much faster (~ 3 orders of magnitude) than using decimal.Decimal(value).log10(). 

  In [1]: convert_to_log10("1.93e-780")
  Out[1]: -779.7144426909922
  """

  import re
  from math import log10

  p = re.compile("([\d\.\-]+)([\sxeE]*)([0-9\-]*)")
  base, _, exponent = p.search(value).groups()
  base = float(base)

  if exponent != '':
    exponent = float(exponent)
  else:
    exponent = 0

  if base == 0:
    return float("-inf")

  lv = log10(float(base)) + float(exponent)
  return lv

