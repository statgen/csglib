import unittest
import gzip
import math
from csg.genetics.ld.pyld import LD

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

   def test_get_variant_haplotypes(self):
      ld = LD()
      ld.add_vcf('1000G_phase3.EUR.chr20.vcf.gz')

      haplotypes = ld.get_variant_haplotypes('21', 11650214)
      self.assertEqual(haplotypes.size, 0)

      haplotypes = ld.get_variant_haplotypes('20', 11111111)
      self.assertEqual(haplotypes.size, 0)

      haplotypes = ld.get_variant_haplotypes('20', 11650214)
      self.assertEqual(haplotypes.size, 1)
      self.assertListEqual(haplotypes.chrom, ['20'])
      self.assertListEqual(haplotypes.position, [11650214])
      self.assertListEqual(haplotypes.ref, ['G'])
      self.assertListEqual(haplotypes.alt, ['A'])
      ld.release_vcfs()

   def test_get_variant_haplotypes_strict(self):
      ld = LD()
      ld.add_vcf('1000G_phase3.EUR.chr20.vcf.gz')

      haplotypes = ld.get_variant_haplotypes_strict('21', 11650214, 'G', 'A')
      self.assertEqual(haplotypes.size, 0)

      haplotypes = ld.get_variant_haplotypes_strict('20', 11111111, 'G', 'A')
      self.assertEqual(haplotypes.size, 0)

      haplotypes = ld.get_variant_haplotypes_strict('20', 11650214, 'G', 'T')
      self.assertEqual(haplotypes.size, 0)

      haplotypes = ld.get_variant_haplotypes_strict('20', 11650214, 'G', 'A')
      self.assertEqual(haplotypes.size, 1)
      self.assertListEqual(haplotypes.chrom, ['20'])
      self.assertListEqual(haplotypes.position, [11650214])
      self.assertListEqual(haplotypes.ref, ['G'])
      self.assertListEqual(haplotypes.alt, ['A'])
      ld.release_vcfs()

   def test_get_region_haplotypes(self):
      ld = LD()
      ld.add_vcf('1000G_phase3.EUR.chr20.vcf.gz')

      haplotypes = ld.get_region_haplotypes('21', 14403183, 19485821)
      self.assertEqual(haplotypes.size, 0)

      haplotypes = ld.get_region_haplotypes('20', 11111111, 11111112)
      self.assertEqual(haplotypes.size, 0)

      haplotypes = ld.get_region_haplotypes('20', 14403183, 19485821)
      self.assertEqual(haplotypes.size, 3)
      self.assertListEqual(haplotypes.chrom, ['20', '20', '20'])
      self.assertListEqual(haplotypes.position, [14403183, 16655993, 19485821])
      self.assertListEqual(haplotypes.ref, ['C', 'T', 'A'])
      self.assertListEqual(haplotypes.alt, ['A', 'G', 'G'])
      ld.release_vcfs()

   def test_compute_variant_freq(self):
       ld = LD()
       ld.add_vcf('1000G_phase3.EUR.chr20.vcf.gz')

       haplotypes = ld.get_variant_haplotypes('20', 11111111)
       freq = ld.compute_freq(haplotypes)
       self.assertListEqual(list(freq), [])

       haplotypes = ld.get_variant_haplotypes('20', 14403183)
       freq = ld.compute_freq(haplotypes)
       id1 = haplotypes.chrom[0] + '_' + str(haplotypes.position[0])
       self.assertAlmostEqual(freq, self.precomputed_freq[id1])

       haplotypes = ld.get_variant_haplotypes('20', 11650214)
       freq = ld.compute_freq(haplotypes)
       id1 = haplotypes.chrom[0] + '_' + str(haplotypes.position[0])
       self.assertAlmostEqual(freq, self.precomputed_freq[id1])
       ld.release_vcfs()

   def test_compute_region_freq(self):
      ld = LD()
      ld.add_vcf('1000G_phase3.EUR.chr20.vcf.gz')

      haplotypes = ld.get_region_haplotypes('20', 11111111, 11111112)
      freq = ld.compute_freq(haplotypes)
      self.assertListEqual(list(freq), [])

      haplotypes = ld.get_region_haplotypes('20', 11650214, 60759931)
      freq = ld.compute_freq(haplotypes)
      for i in xrange(0, haplotypes.size):
         id1 = haplotypes.chrom[i] + '_' + str(haplotypes.position[i])
         self.assertAlmostEqual(freq[i], self.precomputed_freq[id1], places = 6)
      ld.release_vcfs()

   def test_compute_r_cross(self):
       ld = LD()
       ld.add_vcf('1000G_phase3.EUR.chr20.vcf.gz')

       haplotypes1 = ld.get_variant_haplotypes('20', 1111112)
       haplotypes2 = ld.get_variant_haplotypes('21', 1111112)
       r = ld.compute_r_cross(haplotypes1, haplotypes2)
       self.assertListEqual(list(r), [])

       haplotypes1 = ld.get_variant_haplotypes('20', 1111112)
       haplotypes2 = ld.get_variant_haplotypes('20', 16655993)
       r = ld.compute_r_cross(haplotypes1, haplotypes2)
       self.assertListEqual(list(r), [])

       haplotypes1 = ld.get_variant_haplotypes('20', 11650214)
       haplotypes2 = ld.get_variant_haplotypes('20', 16655993)
       id1 = haplotypes1.chrom[0] + '_' + str(haplotypes1.position[0])
       id2 = haplotypes2.chrom[0] + '_' + str(haplotypes2.position[0])
       expected_rsquare = self.precomputed_ld[id1][id2]
       r = ld.compute_r_cross(haplotypes1, haplotypes2)
       self.assertTupleEqual(r.shape, (1, 1))
       self.assertAlmostEqual(r[0, 0] ** 2, expected_rsquare)

       haplotypes1 = ld.get_region_haplotypes('20', 11650214, 14403183)
       haplotypes2 = ld.get_region_haplotypes('20', 16655993, 19485821)
       r = ld.compute_r_cross(haplotypes1, haplotypes2)
       freq1 = ld.compute_freq(haplotypes1)
       freq2 = ld.compute_freq(haplotypes2)
       self.assertTupleEqual(r.shape, (2, 2))
       for i in xrange(0, haplotypes1.size):
          for j in xrange(0, haplotypes2.size):
             id1 = haplotypes1.chrom[i] + '_' + str(haplotypes1.position[i])
             id2 = haplotypes2.chrom[j] + '_' + str(haplotypes2.position[j])
             if id1 != id2:
                expected_rsquare = self.precomputed_ld[id1][id2]
                if math.isnan(r[i, j]):
                   self.assertTrue(math.isnan(expected_rsquare))
                else:
                   self.assertAlmostEqual(r[i, j] ** 2, expected_rsquare)
             else:
                if freq1[i] == 0.0 or freq1[i] == 1.0 or freq2[j] == 0.0 or freq2[j] == 0.0:
                   self.assertTrue(math.isnan(r[i, j]))
                else:
                   self.assertAlmostEqual(r[i, j] ** 2, 1.0)
       ld.release_vcfs()

   def test_compute_r_pairwise(self):
      ld = LD()
      ld.add_vcf('1000G_phase3.EUR.chr20.vcf.gz')

      haplotypes = ld.get_region_haplotypes('20', 11111111, 11111112)
      r = ld.compute_r_pairwise(haplotypes)
      self.assertEqual(haplotypes.size, 0)

      haplotypes = ld.get_region_haplotypes('20', 11650214, 60759931)
      r = ld.compute_r_pairwise(haplotypes)
      freq = ld.compute_freq(haplotypes)
      for i in xrange(0, haplotypes.size):
         for j in xrange(i, haplotypes.size):
            id1 = haplotypes.chrom[i] + '_' + str(haplotypes.position[i])
            id2 = haplotypes.chrom[j] + '_' + str(haplotypes.position[j])
            if id1 != id2:
               expected_rsquare = self.precomputed_ld[id1][id2]
               if math.isnan(r[i, j]):
                  self.assertTrue(math.isnan(expected_rsquare))
               else:
                  self.assertAlmostEqual(r[i, j] ** 2, expected_rsquare)
            else:
               if freq[i] == 0.0 or freq[i] == 1.0:
                  self.assertTrue(math.isnan(r[i, j]))
               else:
                  self.assertAlmostEqual(r[i, j] ** 2, 1.0)
      ld.release_vcfs()


if __name__ == '__main__':
   suite = unittest.TestSuite()
   test_loader = unittest.TestLoader()
   suite.addTest(test_loader.loadTestsFromTestCase(TestLD))
   unittest.TextTestRunner(verbosity = 2).run(suite)
