from contextlib import closing
from collections import namedtuple
import pysam

Alleles = namedtuple('Alleles', ['ref', 'alt'])

class Panel:

   def __init__(self, vcf, left_range = 500000, right_range = 500000):
      self.vcf = vcf
      self.chrom = ''
      self.start = 0
      self.end = 0
      self.variants = None
      self.left_range = left_range
      self.right_range = right_range

   def get_alleles(self, chrom, position):
      if self.chrom != chrom or self.start > position or position > self.end:
         self.variants = dict()
         self.chrom = chrom
         self.start = position - self.left_range
         self.start = self.start - 1 if self.start > 0 else 0
         self.end = position + self.right_range
         with closing(pysam.Tabixfile(self.vcf)) as tabix:
            for row in tabix.fetch(self.chrom, self.start, self.end):
               fields = row.split('\t', 5)
               bp = long(fields[1])
               if bp not in self.variants:
                  self.variants[bp] = list()
               self.variants[bp].append(Alleles(fields[3], set(fields[4].split(','))))
      return self.variants.get(position, None)

