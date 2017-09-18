from __future__ import absolute_import
import math, re
import numpy as np

def convert_to_log10(value,na_strings=["NA","NaN",".",""]):
  """
  Take the log10 of a **string** representing either a plain number,
  or number in scientific notation. This can handle arbitrary precision.

  In [1]: convert_to_log10("1.93e-780")
  Out[1]: -779.7144426909922
  """

  if value in na_strings:
    return np.nan

  p = re.compile("([\d\.\-]+)([\sxeE]*)([0-9\-]*)")
  base, _, exponent = p.search(value).groups()
  base = float(base)

  if exponent != '':
    exponent = float(exponent)
  else:
    exponent = 0

  if base == 0:
    return float("-inf")

  lv = math.log10(float(base)) + float(exponent)
  return lv
