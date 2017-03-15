from bx.intervals import IntervalTree, Interval 

class ChromTree:
  def __init__(self):
    self.chrom = {}

  def add_interval(self,chrom,start,end,value=None):
    chrom = str(chrom).replace("chr","")
    self.chrom.setdefault(chrom,IntervalTree()).insert_interval(Interval(start,end,value))

  def find_overlap(self,chrom,start,end=None):
    chrom = str(chrom).replace("chr","")
    node = self.chrom.get(chrom)

    if node is None:
      return None

    else:
      if end is None:
        end = start

      match = node.find(start,end)

      return match
