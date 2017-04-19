import sys
import argparse
import gzip
from collections import namedtuple
import operator

argparser = argparse.ArgumentParser(description = 'Finds nuclear families from PC-relate output.')
argparser.add_argument('--kinship', metavar = 'file', dest = 'in_kinship_file', required = True, help = 'File with estimated kinship coefficients (KING or PC-relate format).')
argparser.add_argument('--format', metavar = 'name', dest = 'in_format', required = True, choices = ['pc-relate', 'king'], help = 'File format: pc-relate, king')
argparser.add_argument('--sex', metavar = 'file', dest = 'in_sex_file', required = True, help = 'File with two tab-separated columns: id, sex (M/F).')
argparser.add_argument('--out', metavar = 'file', dest = 'out_file', required = True, help = 'Output file.')

pcrelate_header = ['ID1', 'ID2', 'nsnp',  'kin',  'k0',  'k1',  'k2']
king_header = ['FID1', 'ID1', 'FID2', 'ID2', 'N_SNP', 'HetHet', 'IBS0', 'Kinship']

relatives_by_individual = dict()
sex = dict()
trios = list()

Entry = namedtuple('Entry', ['id1', 'id2', 'kinship', 'prob_ibd0', 'prop_ibs0'])
Degree = namedtuple('Degree', ['degree', 'low', 'high'])
Relation = namedtuple('Relation', ['kinship', 'degree', 'prob_ibd0', 'prop_ibs0'])
Trio = namedtuple('Trio', ['individual', 'father', 'mother', 'sex'])

# Based on Manichaikul et al., i degree relationship if kinship coefficient is in the interval (2^-(d + 3/2), 2^-(d + 1/2))
# d = 0 -- monozygotic twins
kinship_thresholds = []
for d in xrange(0, 4):
   kinship_thresholds.append(Degree( degree = d, low = 2**(-d - 1.5), high = 2**(-d - 0.5) ))


def pcrelate_reader(in_kinship_file):
   with open(in_kinship_file, 'r') as ifile:
      line = ifile.readline()
      header = line.rstrip().split('\t')
      if len(header) != len(pcrelate_header) or any([obs != exp for obs, exp in zip(header, pcrelate_header)]):
         raise Exception('Kinship file header does not match the expected PC-relate format (tab-separated): %s!' % (', '.join(pcrelate_header)))  
      id1_idx = header.index('ID1')
      id2_idx = header.index('ID2')
      kinship_idx = header.index('kin')
      prob_ibd0_idx = header.index('k0') # probability of sharing zero alleles IBD
      for line in ifile:
         row = line.rstrip().split('\t')
         yield Entry(row[id1_idx], row[id2_idx], float(row[kinship_idx]), float(row[prob_ibd0_idx]), None)


def king_reader(in_kinship_file): 
   with open(in_kinship_file, 'r') as ifile:
      line = ifile.readline()
      header = line.rstrip().split('\t')
      if len(header) != len(king_header) or any([obs != exp for obs, exp in zip(header, king_header)]):
         raise Exception('Kinship file header does not match the expected KING format (tab-separated): %s!' % (', '.join(king_header)))
      id1_idx = header.index('ID1')
      id2_idx = header.index('ID2')
      kinship_idx = header.index('Kinship')
      prop_ibs0_idx = header.index('IBS0') # proportion of SNPs with zero IBS

      for line in ifile:
         row = line.rstrip().split('\t')
         yield Entry(row[id1_idx], row[id2_idx], float(row[kinship_idx]), None, float(row[prop_ibs0_idx]))


def read_sex(in_sex_file):
   with open(in_sex_file, 'r') as ifile:
      for line in ifile:
         row = line.rstrip().split('\t')
         if row[1] != 'M' and row[1] != 'F':
            continue
         sex[row[0]] = row[1]
      

def read_kinship(in_kinship_file, in_format):
   if in_format == 'pc-relate':
      reader = pcrelate_reader(in_kinship_file)
   elif in_format == 'king':
      reader = king_reader(in_kinship_file)

   for entry in reader:
      degree = None
      for threshold in kinship_thresholds:
         if entry.kinship > threshold.low and entry.kinship <= threshold.high:
            degree = threshold.degree
            break

      if degree is None:
         continue

      id1_relatives = relatives_by_individual.get(entry.id1, None)
      if id1_relatives is None:
         id1_relatives = dict()
         relatives_by_individual[entry.id1] = id1_relatives

      id2_relatives = relatives_by_individual.get(entry.id2, None)
      if id2_relatives is None:
         id2_relatives = dict()
         relatives_by_individual[entry.id2] = id2_relatives

      if not entry.id2 in id1_relatives:
         id1_relatives[entry.id2] = Relation(entry.kinship, degree, entry.prob_ibd0, entry.prop_ibs0)

      if not entry.id1 in id2_relatives:
         id2_relatives[entry.id1] = Relation(entry.kinship, degree, entry.prob_ibd0, entry.prop_ibs0)


def create_trios(individuals):
   possible_parents_by_individual = dict()

   # 1. for every individual keep all 1st degree relatives that may be parents
   for id1, relatives in relatives_by_individual.iteritems():
      possible_parents_by_individual[id1] = list()   
      for id2, relation in relatives.iteritems():
         # Parent-offspring thresholds are from Conomos et. al., 2016
         if relation.degree == 1 and (relation.prob_ibd0 is not None and relation.prob_ibd0 < 2**(-4.5) or relation.prop_ibs0 is not None and relation.prop_ibs0 < 0.005):
             possible_parents_by_individual[id1].append(id2)

   # 2. Drop all individuals that can't form trio family
   for id1 in relatives_by_individual.iterkeys():
      if len(possible_parents_by_individual[id1]) < 2:
          del possible_parents_by_individual[id1]

   # 3. Two individuals are biological parents if (1) they are not related (degree > 3), and (2) they have opposite sex. Trio is created if there is only one possible pair of biological parents. 
   for id1, possible_parents in possible_parents_by_individual.iteritems():
      unrelated_parents = set()
      for i in xrange(0, len(possible_parents)):
         for j in xrange(i + 1, len(possible_parents)):
            id_i = possible_parents[i]
            id_j = possible_parents[j]
            if id_i not in sex or id_j not in sex:
               continue
            if sex[id_i] == sex[id_j]:
               continue
            relation = relatives_by_individual[id_i].get(id_j, None)
            if relation is None or relation.degree > 3:
               unrelated_parents.add(id_i)
               unrelated_parents.add(id_j)
      if len(unrelated_parents) == 2:
         father = None
         mother = None
         for parent in unrelated_parents:
            if sex[parent] == 'M':
               father = parent
            else:
               mother = parent
         trio = Trio(id1, father, mother, sex.get(id1, 'NA'))
         trios.append(trio)


def write_trios(out_file):
   with open(out_file, 'w') as ofile:
      ofile.write('IID\tFID\tMID\tSEX\n')
      for trio in trios:
         ofile.write('%s\t%s\t%s\t%s\n' % (trio.individual, trio.father, trio.mother, trio.sex))


if __name__ == '__main__':
   args = argparser.parse_args()

   read_kinship(args.in_kinship_file, args.in_format)
   read_sex(args.in_sex_file)
   create_trios(relatives_by_individual)
   write_trios(args.out_file)         

