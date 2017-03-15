import pandas as pd

def chrom_int(c):
  """
  Convert a chromosome string (e.g. chr9 or chrX) to an integer
  
  The conversion ends up as: 
    X -> 23
    Y -> 24
    MT -> 25
    XY -> 26

  The chromosome string can be either chromX, chrX, or just X, for example. 

  """

  c = str(c).replace("chrom","").replace("chr","").upper()

  if c == "X":
    return 23
  elif c == "Y":
    return 24
  elif c == "MT":
    return 25
  elif c == "XY":
    return 26
  else:
    try:
      c = int(c)
    except:
      raise Exception, "Bad chromosome: %s" % c

    return c

def parse_bp(fstr):
  """
  Parse a unit of genomic length from a string. Possible strings: 

  1.25 MB or 1.25 mb = 1,250,000 bp
  300 kb = 300,000 bp
  20000 = 20,000 bp
  """
  
  # Normalize first by lower-casing the unit and removing commas (if they were used instead of dots.) 
  spl = fstr.lower().replace(",","")

  try:
    value, unit = re.search("([\d\.]+)\s*(mb|kb|bp)?",spl).groups()
  except:
    raise ValueError, "Invalid unit of genomic position: {}".format(fstr)

  if unit == 'mb':
    bp = int(float(value) * 1E6)
  elif unit == 'kb':
    bp = int(float(value) * 1E3)
  elif unit == 'bp' or unit is None:
    bp = int(value)

  return bp

def assign_locus(df,col_chrom,col_pos,locus_name="LOCUS",dist=250000):
  """
  Assigns a locus number to each variant in a data frame of variants with chromosome and position specified. 

  A column will be inserted containing the locus number assigned to each variant. 

  Args:
    df: data frame, each row is a variant.
    chrom: chromosome column name.
    pos: position column name.
    locus_name (str): What should we call the locus number column? 
    dist (int): base-pair distance between variants required to assign a new locus number

  Returns:
    DataFrame: with an added column specifying the locus number
  """

  if locus_name in df.columns:
    del df[locus_name]

  srt = df.sort([col_chrom,col_pos])
  srt.insert(0,locus_name,pd.np.nan)

  chunks = []
  cur_loc = 1
  for chrom, grpdf in srt.groupby(col_chrom):
    cur_pos = None
    for index, row in grpdf.iterrows():
      if cur_pos is None:
        cur_pos = row[col_pos]
      else:
        if abs(row[col_pos] - cur_pos) > dist:
          cur_loc += 1

        cur_pos = row[col_pos]

      srt.ix[index,locus_name] = cur_loc

    cur_loc += 1

  srt[locus_name] = srt[locus_name].astype("int")

  return srt

