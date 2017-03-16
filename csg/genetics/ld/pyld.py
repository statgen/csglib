import pysam
import numpy as np

chunk_M = 1000 # initial number of variants in chunk
incr_M = 1000 # chunk increment size if needed

N = None # number of haplotypes
haplotypes_all = None
haplotypes_per_chrom = dict()

def open_vcf(vcf_path, chrom = None):
   global N
   global haplotypes_all

   close_vcf(chrom)
   if chrom is None or chrom.upper() == 'ALL':
      haplotypes_all = pysam.Tabixfile(vcf_path)
      header = haplotypes_all.header 
   else:
      haplotypes_per_chrom[chrom] = pysam.Tabixfile(vcf_path)
      header = haplotypes_per_chrom[chrom].header
   n_samples = None
   for row in header:
      if row.startswith('#CHROM'):
         n_samples = len(row.split('\t')) - 9
   if not n_samples:
      raise Exception('Header was not found in the VCF: %s' % vcf_path)
   if N is None:
      N = 2 * n_samples
   elif N != 2 * n_samples:
      raise Exception('Number of samples in the opened VCF is not consistent with previous VCF\'s.')

def close_vcf(chrom = None):
   global N
   global haplotypes_all
   if chrom is None or chrom.upper() == 'ALL':
      if not haplotypes_all is None:
         haplotypes_all.close()
         haplotypes_all = None
      for chrom, haplotypes in haplotypes_per_chrom.iteritems():
         haplotypes.close()
      haplotypes_per_chrom.clear()
   elif chrom in haplotypes_per_chrom:
      haplotypes_per_chrom[chrom].close()
      del haplotypes_per_chrom[chrom]
   if haplotypes_all is None and len(haplotypes_per_chrom) == 0:
      N = None

def get_variant_haplotype(chrom, position):
   tabix = haplotypes_per_chrom.get(chrom, haplotypes_all)
   if tabix:
      for row in tabix.fetch(chrom, position - 1, position, parser = pysam.asTuple()):
         # ignore multi-allelic
         if ',' in row[4]:
            continue
         haplotype = np.empty([1, N], dtype = np.int, order = 'C')
         for i in xrange(9, len(row)):
            haplotype[0, 2 * i - 18] = np.int(row[i][0])
            haplotype[0, 2 * i - 17] = np.int(row[i][2])
         return ([row[0] + '_' + row[1] + '_' + row[3] + '_' + row[4]], haplotype)
   return None

def get_region_haplotypes(chrom, start_position, end_position):
   tabix = haplotypes_per_chrom.get(chrom, haplotypes_all)
   if tabix:
      max_M = chunk_M
      markers = []
      haplotypes = np.empty([max_M, N], dtype = np.int, order = 'C')
      M = 0
      for row in tabix.fetch(chrom, start_position - 1, end_position, parser = pysam.asTuple()):
         # ignore multi-allelic
         if ',' in row[4]:
            continue
         markers.append(row[0] + '_' + row[1] + '_' + row[3] + '_' + row[4])
         if M >= max_M:
            max_M += incr_M
            haplotypes.resize(max_M, N)
         for i in xrange(9, len(row)):
            haplotypes[M, 2 * i - 18] = np.int(row[i][0])
            haplotypes[M, 2 * i - 17] = np.int(row[i][2])
         M += 1
      if M > 0:
         if M != max_M:
            haplotypes.resize(M, N)
         return (markers, haplotypes)
   return None

def get_freq(haplo):
   return haplo[1].sum() / float(haplo[1].size)

def get_freq_array(haplos):
   m, n = haplos[1].shape
   ones = np.ones((n, 1), dtype = np.int)
   freqs = np.dot(haplos[1], ones) / float(n)
   freqs.resize(m)
   return freqs

def compute_r(haplo1, haplo2):
   n = haplo1[1].size
   a_counts = haplo1[1].sum()
   b_counts = haplo2[1].sum()
   return ((n * np.dot(haplo1[1], haplo2[1].transpose()) - a_counts * b_counts) / np.sqrt( a_counts * b_counts * (n - a_counts) * (n - b_counts)))[0, 0]

def compute_r_matrix(haplos):
   m, n = haplos[1].shape
   ones = np.ones((n, 1), dtype = np.int)
   a_counts = np.dot(haplos[1], ones)
   b_counts = n - a_counts
   c1 = np.dot(a_counts, a_counts.transpose())
   return (n * np.dot(haplos[1], haplos[1].transpose()) - c1) / np.sqrt(c1 * np.dot(b_counts, b_counts.transpose()))

