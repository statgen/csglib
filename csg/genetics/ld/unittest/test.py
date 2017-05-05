import unittest
import gzip
import math
from csg.genetics.ld import pyld

class TestLD(unittest.TestCase):

   @classmethod
   def setUpClass(self):
      self.precomputed_ld = dict()
      self.precomputed_freq = dict()
      with gzip.GzipFile('1000G_phase3.EUR.chr20.LD_vcftools.tbl.gz', 'r') as ifile:
         for line in ifile:
            header = line.rstrip().split('\t')
            break
         for line in ifile:
            fields = dict(zip(header, line.rstrip().split('\t')))
            id1 = fields['CHR'] + '_' + fields['POS1']
            id2 = fields['CHR'] + '_' + fields['POS2']
            rsq = float(fields['R^2'])
            d = float(fields['D'])
            dprime = float(fields['Dprime'])
            if id1 not in self.precomputed_ld:
                self.precomputed_ld[id1] = dict()
            if id2 not in self.precomputed_ld:
                self.precomputed_ld[id2] = dict()
            self.precomputed_ld[id1][id2] = rsq
            self.precomputed_ld[id2][id1] = rsq
      with gzip.GzipFile('1000G_phase3.EUR.chr20.freq_vcftools.tbl.gz', 'r') as ifile:
         for line in ifile:
            header = line.rstrip().split('\t')
            break
         for line in ifile:
            fields = dict(zip(header, line.rstrip().split('\t', len(header) - 1)))
            id1 = fields['CHROM'] + '_' + fields['POS']
            alt_freq = float(fields['{FREQ}'].split('\t')[1])
            self.precomputed_freq[id1] = alt_freq


   def test_freq_array(self):
      pyld.open_vcf('1000G_phase3.EUR.chr20.vcf.gz', chrom = '20')
      haplotypes = pyld.get_region_haplotypes('20', 11650214, 60759931)
      self.assertEqual(len(haplotypes[0]), 9)
      self.assertTupleEqual(haplotypes[1].shape, (9, 1006))
      freq = pyld.get_freq_array(haplotypes)
      for i in xrange(0, len(haplotypes[0])):
         id1 = haplotypes[0][i].rsplit('_', 2)[0]
         self.assertAlmostEqual(freq[i], self.precomputed_freq[id1], places = 6)
      pyld.close_vcf()


   def test_r_matrix(self):
      pyld.open_vcf('1000G_phase3.EUR.chr20.vcf.gz', chrom = '20')
      haplotypes = pyld.get_region_haplotypes('20', 11650214, 60759931)
      self.assertEqual(len(haplotypes[0]), 9)
      self.assertTupleEqual(haplotypes[1].shape, (9, 1006))
      r = pyld.compute_r_matrix(haplotypes)
      for i in xrange(0, len(haplotypes[0])):
         for j in xrange(i, len(haplotypes[0])):
            id1 = haplotypes[0][i].rsplit('_', 2)[0]
            id2 = haplotypes[0][j].rsplit('_', 2)[0]
            if id1 != id2:
               expected_r = self.precomputed_ld[id1][id2]
               if math.isnan(r[i, j]):
                  self.assertTrue(math.isnan(expected_r))
               else:
                  self.assertAlmostEqual(r[i, j] ** 2, expected_r)
            else:
               pass
               #print id1, id2, r[i, j]
               #self.assertEqual(r[i, j] ** 2, 1.0)
      pyld.close_vcf()


if __name__ == '__main__':
   suite = unittest.TestSuite()
   test_loader = unittest.TestLoader()
   suite.addTest(test_loader.loadTestsFromTestCase(TestLD))
   unittest.TextTestRunner(verbosity = 2).run(suite)
