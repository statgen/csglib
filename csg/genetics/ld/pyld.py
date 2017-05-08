import pysam
import numpy as np
from collections import namedtuple

np.seterr(invalid = 'ignore') # ignore division by zero in matrix operations.

Haplotypes = namedtuple('Haplotypes', ['size', 'chrom', 'position', 'ref', 'alt', 'haplotypes'], verbose = False)

class LD:

   chunk_M = 1000 # initial number of variants in chunk
   incr_M = 1000 # chunk increment size if needed

   def __init__(self):
      self.N = None
      self.chrom2tabix = dict()

   def add_vcf(self, vcf_path):
      tabix = pysam.Tabixfile(vcf_path)
      n_samples = None
      for row in tabix.header:
         if row.startswith('#CHROM'):
            n_samples = len(row.split('\t')) - 9
      if not n_samples:
         raise Exception('Header was not found in the VCF: %s' % vcf_path)
      if self.N is None:
         self.N = 2 * n_samples
      elif self.N != 2 * n_samples:
         raise Exception('Number of samples in the opened VCF is not consistent with previous VCF\'s.')
      for chrom in tabix.contigs:
         if chrom not in self.chrom2tabix:
            self.chrom2tabix[chrom] = tabix

   def release_vcfs(self):
      for chrom, tabix in self.chrom2tabix.iteritems():
         tabix.close()
      self.chrom2tabix.clear()
      self.N = None

   def get_variant_haplotypes(self, chrom, position):
      tabix = self.chrom2tabix.get(chrom, None)
      if tabix is not None:
         for row in tabix.fetch(chrom, position - 1, position, parser = pysam.asTuple()):
            if ',' in row[4]: # ignore multi-allelic
               continue
            haplotypes = np.empty([1, self.N], dtype = np.int, order = 'C')
            for i in xrange(9, len(row)):
               haplotypes[0, 2 * i - 18] = np.int(row[i][0])
               haplotypes[0, 2 * i - 17] = np.int(row[i][2])
            return Haplotypes(size = 1, chrom = [row[0]], position = [long(row[1])], ref = [row[3]], alt = [row[4]], haplotypes = haplotypes)
      return None

   def get_variant_haplotypes_strict(self, chrom, position, ref, alt):
      tabix = self.chrom2tabix.get(chrom, None)
      if tabix is not None:
         for row in tabix.fetch(chrom, position - 1, position, parser = pysam.asTuple()):
            if ',' in row[4]: # ignore multi-allelic
               continue
            if row[3] != ref or row[4] != alt:
               continue
            haplotypes = np.empty([1, self.N], dtype = np.int, order = 'C')
            for i in xrange(9, len(row)):
               haplotypes[0, 2 * i - 18] = np.int(row[i][0])
               haplotypes[0, 2 * i - 17] = np.int(row[i][2])
            return Haplotypes(size = 1, chrom = [row[0]], position = [long(row[1])], ref = [row[3]], alt = [row[4]], haplotypes = haplotypes)
      return None

   def get_region_haplotypes(self, chrom, start_position, end_position):
      tabix = self.chrom2tabix.get(chrom, None)
      if tabix:
         max_M = self.chunk_M
         variant_chrom = []
         variant_position = []
         variant_ref = []
         variant_alt= []
         haplotypes = np.empty([max_M, self.N], dtype = np.int, order = 'C')
         M = 0
         for row in tabix.fetch(chrom, start_position - 1, end_position, parser = pysam.asTuple()):
            if ',' in row[4]: # ignore multi-allelic
               continue
            variant_chrom.append(row[0])
            variant_position.append(long(row[1]))
            variant_ref.append(row[3])
            variant_alt.append(row[4])
            if M >= max_M:
               max_M += self.incr_M
               haplotypes.resize(max_M, self.N)
            for i in xrange(9, len(row)):
               haplotypes[M, 2 * i - 18] = np.int(row[i][0])
               haplotypes[M, 2 * i - 17] = np.int(row[i][2])
            M += 1
         if M > 0:
            if M != max_M:
               haplotypes.resize(M, self.N)
            return Haplotypes(size = M, chrom = variant_chrom, position = variant_position, ref = variant_ref, alt = variant_alt, haplotypes = haplotypes)
      return None

   def compute_variant_freq(self, haplotypes):
      if haplotypes is None:
          return None
      return haplotypes.haplotypes.sum() / float(haplotypes.haplotypes.size)

   def compute_region_freq(self, haplotypes):
      if haplotypes is None:
         return None
      m, n = haplotypes.haplotypes.shape
      ones = np.ones((n, 1), dtype = np.int)
      freqs = np.dot(haplotypes.haplotypes, ones) / float(n)
      freqs.resize(m)
      return freqs

   def compute_r(self, haplotypes1, haplotypes2):
      if haplotypes1 is None:
         return None
      if haplotypes2 is None:
         return None
      n = haplotypes1.haplotypes.size
      a_counts = haplotypes1.haplotypes.sum()
      b_counts = haplotypes2.haplotypes.sum()
      return ((n * np.dot(haplotypes1.haplotypes, haplotypes2.haplotypes.transpose()) - a_counts * b_counts) / np.sqrt( a_counts * b_counts * (n - a_counts) * (n - b_counts)))[0, 0]

   def compute_r_matrix(self, haplotypes):
      if haplotypes is None:
         return None
      m, n = haplotypes.haplotypes.shape
      ones = np.ones((n, 1), dtype = np.int)
      a_counts = np.dot(haplotypes.haplotypes, ones)
      b_counts = n - a_counts
      c1 = np.dot(a_counts, a_counts.transpose())
      return (n * np.dot(haplotypes.haplotypes, haplotypes.haplotypes.transpose()) - c1) / np.sqrt(c1 * np.dot(b_counts, b_counts.transpose()))

