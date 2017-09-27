import pysam

def find_possible_chunks(fpath):
  """
  Given a tabix-indexed file, figure chunks (chr:start-end) that cover data in the file. 
  Could be used to load results in a parallel fashion.
  """
  
  tb = pysam.TabixFile(fpath)
  chromosomes = list((str(s) for s in range(1,23)))
  chunks = dict()
  size = int(10E6)
  continues = 1000
  
  for chrom in chromosomes:
    cont_i = 0
    s = 1
    e = size
    while 1:
      try:
        it = tb.fetch(chrom,s,e)
      except ValueError:
        # Reached when start or end exceed 32-bit int max
        break
    
      try:
        it.next()
      except StopIteration:
        if cont_i > continues:
          break
        else:
          cont_i += 1
      else:
        chunks.setdefault(chrom,[]).append((s,e))
        cont_i = 0
    
      s += size
      e += size
        
  return chunks
