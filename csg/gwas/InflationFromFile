import argparse
import gzip
import os
from csg.gwas.inflation import *
from csg.intervaltree.IntervalTree import *

argparser = argparse.ArgumentParser('Computes inflation factor.')
argparser.add_argument('--file', metavar = 'file', dest = 'in_file', required = True, help = 'Input file compressed using gzip/bgzip.')
argparser.add_argument('--pvalue', metavar = 'name', dest = 'pvalue_field', required = False, help = 'Field that stores p-values.')
argparser.add_argument('--effect', metavar = 'name', dest = 'effect_field', required = False, help = 'Field that stores effect sizes.')
argparser.add_argument('--se', metavar = 'name', dest = 'se_field', required = False, help = 'Field that stores standard error.')
argparser.add_argument('--chrom', metavar = 'name', dest = 'chrom_field', required = False, help = 'Field that stores chromosome name.')
argparser.add_argument('--position', metavar = 'number', dest = 'position_field', required = False, help = 'Field that stores chromosomal position in base-pairs.')
argparser.add_argument('--sep', metavar = 'separator', dest = 'sep', required = True, choices = ['tab', 'comma', 'semicolon'], help = 'Field separator: tab, comma, semicolon.')
argparser.add_argument('--known-regions-file', metavar = 'file', dest = 'known_regions_file', required = False, help = 'File with regions that have known association hits. File has three columns')

def read_known_regions(known_regions_file):
   regions_by_chrom = dict()
   with open(known_regions_file, 'r') as ifile:
      for line in ifile:
         fields = line.rstrip().split('\t')
         chrom = fields[0]
         position_start = long(fields[1])
         position_end = long(fields[2])
         if chrom not in regions_by_chrom:
            regions_by_chrom[chrom] = IntervalTree()
         regions_by_chrom[chrom].add(position_start, position_end)
   return regions_by_chrom

def read_pvalue(in_file, pvalue_field):
   with gzip.GzipFile(in_file, 'r') as ifile:
      header = ifile.readline().strip().split('\t')
      try:
         pvalue_idx = header.index(pvalue_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % pvalue_field)
      for line in ifile:
         fields = line.rstrip().split('\t')
         yield AssociationStats(pvalue = np.float(fields[pvalue_idx]), effect = None, se = None)

def read_pvalue_filter_regions(in_file, pvalue_field, regions_by_chrom, chrom_field, position_field):
   with gzip.GzipFile(in_file, 'r') as ifile:
      header = ifile.readline().strip().split('\t')
      try:
         pvalue_idx = header.index(pvalue_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % pvalue_field)
      try:
         chrom_idx = header.index(chrom_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % chrom_field)
      try:
         position_idx = header.index(position_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % position_field)
      for line in ifile:
         fields = line.rstrip().split('\t')
         position = long(fields[position_idx])
         chrom = fields[chrom_idx]
         regions = regions_by_chrom.get(chrom, None)
         if regions is not None and any(regions.point_intersect(position)):
            continue
         yield AssociationStats(pvalue = np.float(fields[pvalue_idx]), effect = None, se = None)

def read_effect_se(in_file, effect_field, se_field):
   with gzip.GzipFile(in_file, 'r') as ifile:
      header = ifile.readline().strip().split('\t')
      try:
         effect_idx = header.index(effect_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % effect_field)
      try:
         se_idx = header.index(se_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % se_field)
      for line in ifile:
         fields = line.rstrip().split('\t')
         yield AssociationStats(pvalue = None, effect = np.float(fields[effect_idx]), se = np.float(fields[se_idx]))

def read_effect_se_filter_regions(in_file, effect_field, se_field, regions_by_chrom, chrom_field, position_field):
   with gzip.GzipFile(in_file, 'r') as ifile:
      header = ifile.readline().strip().split('\t')
      try:
         effect_idx = header.index(effect_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % effect_field)
      try:
         se_idx = header.index(se_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % se_field)
      try:
         chrom_idx = header.index(chrom_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % chrom_field)
      try:
         position_idx = header.index(position_field)
      except ValueError:
         raise Exception('Field \'%s\' was not found in input file!' % position_field)
      for line in ifile:
         fields = line.rstrip().split('\t')
         position = long(fields[position_idx])
         chrom = fields[chrom_idx]
         regions = regions_by_chrom.get(chrom, None)
         if regions is not None and any(regions.point_intersect(position)):
            continue
         yield AssociationStats(pvalue = None, effect = np.float(fields[effect_idx]), se = np.float(fields[se_idx]))

if __name__ == '__main__':
   args = argparser.parse_args()

   regions_by_chrom = None

   if args.known_regions_file is not None or args.chrom_field is not None or args.position_field is not None:
      if args.known_regions_file is not None and args.chrom_field is not None and args.position_field is not None:
         regions_by_chrom = read_known_regions(args.known_regions_file)
      else:
         print 'Known regions file must be provided together with chromosome and position field names.'

   if args.pvalue_field is not None:
      if args.effect_field is not None:
         print 'Specified effect column name will be ignored.'
      if args.se_field is not None:
         print 'Specified se column name wil be ignored.'
      if regions_by_chrom:
         stats_generator = read_pvalue_filter_regions(args.in_file, args.pvalue_field, regions_by_chrom, args.chrom_field, args.position_field)
      else:
         stats_generator = read_pvalue(args.in_file, args.pvalue_field)
      inflation = inflation_from_pvalue(stats_generator)
      print os.path.basename(args.in_file), inflation
   elif args.effect_field is not None and args.se_field is not None:
      if regions_by_chrom:
         stats_generator = read_effect_se_filter_regions(args.in_file, args.effect_field, args.se_field, regions_by_chrom, args.chrom_field, args.position_field)
      else:
         stats_generator = read_effect_se(args.in_file, args.effect_field, args.se_field)
      inflation = inflation_from_effect_se(stats_generator)
      print os.path.basename(args.in_file), inflation
   else:
      if args.effect_field is None and args.se_field is None:
         print 'Specify column name for p-values or column names for effects and standard errors.'
      elif args.effect_field is None:
         print 'Specify column name for effect values.'
      elif args.se_field is None:
         print 'Specify column name for se values.'

