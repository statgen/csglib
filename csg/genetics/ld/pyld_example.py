import sys
import yaml
import argparse
import timeit
import pyld

argparser = argparse.ArgumentParser('')
argparser.add_argument('--config', metavar = 'file', dest = 'config_yaml_file', required = True, help = 'Configuration file in YAML format.')

if __name__ == '__main__':
   args = argparser.parse_args()

   config = yaml.load(file(args.config_yaml_file, 'r'))
   if 'vcf-file-path' not in config:
      raise Exception('Configuration file is incorrect! \'vcf-file-path\' is missing.')
   if 'ALL' not in config['vcf-file-path'] and any([str(i) not in config['vcf-file-path'] for  i in xrange(1, 23)]):
      raise Exception('Configuration file is incorrect! Path of one or more VCF files is missing.')

   for key, file_path in config['vcf-file-path'].iteritems():
      pyld.open_vcf(file_path, key)

   # Pairwise LD between 2 bi-allelic SNPs
   chr1 = '22'
   pos1 = 16051249

   chr2 = '22'
   pos2 = 16052167

   haplo1 = pyld.get_variant_haplotype(chr1, pos1)
   haplo2 = pyld.get_variant_haplotype(chr2, pos2)

   r = pyld.compute_r(haplo1, haplo2)

   print '# Pairwise LD (r, r^2) between 2 bi-allelic SNPs'
   print haplo1[0][0], haplo2[0][0], r, r ** 2

   # Pairwise LD (r, r^2) between all 2 bi-allelic SNPs in region
   haplos = pyld.get_region_haplotypes(chr1, pos1, pos2)
   r_matrix = pyld.compute_r_matrix(haplos)

   print '# Pairwise LD (r, r^2) between all 2 bi-allelic SNPs in ', chr1, pos1, '-', pos2
   for i in xrange(0, len(haplos[0])):
      for j in xrange(i, len(haplos[0])):
         print haplos[0][i], haplos[0][j], r_matrix[i, j], r_matrix[i, j] ** 2
   
   pyld.close_vcf()
