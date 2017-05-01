from collections import namedtuple

KinshipEntry = namedtuple('KinshipEntry', ['id1', 'id2', 'kinship', 'prob_ibd0', 'prop_ibs0'])
SexEntry = namedtuple('SexEntry', ['id', 'sex'])
Relation = namedtuple('Relation', ['kinship', 'degree', 'prob_ibd0', 'prop_ibs0'])
Trio = namedtuple('Trio', ['individual', 'father', 'mother', 'sex'])

def kinship2degree(kinship):
   """Returns relation degree derived from genetic kinship value.

   Returns
   -------
   degree: int or None

   Raises
   ------
   Exception
      Raises exception if kinship value is greater than 2^-0.5.

   """

   # Based on Manichaikul et al., i degree relationship if kinship coefficient is in the interval (2^-(d + 3/2), 2^-(d + 1/2))
   # d = 0 -- monozygotic twins
   if kinship > 2**(-0.5):
      raise Exception('Unexpectedly high kiship value %g. Check for duplicated samples.' % kinship)
   for d in xrange(0, 32):
      if kinship > 2**(-d - 1.5) and kinship <= 2**(-d - 0.5):
         return d
   return None


def pcrelate_reader(in_kinship_file):
   pcrelate_header = ['ID1', 'ID2', 'nsnp',  'kin',  'k0',  'k1',  'k2']
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
         yield KinshipEntry(row[id1_idx], row[id2_idx], float(row[kinship_idx]), float(row[prob_ibd0_idx]), None)


def king_reader(in_kinship_file):
   king_header = ['FID1', 'ID1', 'FID2', 'ID2', 'N_SNP', 'HetHet', 'IBS0', 'Kinship']
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
         yield KinshipEntry(row[id1_idx], row[id2_idx], float(row[kinship_idx]), None, float(row[prop_ibs0_idx]))


def sex_reader(in_sex_file):
   with open(in_sex_file, 'r') as ifile:
      for line in ifile:
         row = line.rstrip().split('\t')
         if row[1] != 'M' and row[1] != 'F':
            continue
         yield SexEntry(row[0], row[1])


def get_trios(kinship_entries, sex_entries, min_degree):
   individual_relatives = dict()
   individual_sex = dict()

   for kinship_entry in kinship_entries:
      degree = kinship2degree(kinship_entry.kinship)
      if degree is None or degree > min_degree:
         continue

      id1_relatives = individual_relatives.get(kinship_entry.id1, None)
      if id1_relatives is None:
         id1_relatives = dict()
         individual_relatives[kinship_entry.id1] = id1_relatives

      id2_relatives = individual_relatives.get(kinship_entry.id2, None)
      if id2_relatives is None:
         id2_relatives = dict()
         individual_relatives[kinship_entry.id2] = id2_relatives

      if not kinship_entry.id2 in id1_relatives:
         id1_relatives[kinship_entry.id2] = Relation(kinship_entry.kinship, degree, kinship_entry.prob_ibd0, kinship_entry.prop_ibs0)

      if not kinship_entry.id1 in id2_relatives:
         id2_relatives[kinship_entry.id1] = Relation(kinship_entry.kinship, degree, kinship_entry.prob_ibd0, kinship_entry.prop_ibs0)

   for sex_entry in sex_entries:
      individual_sex[sex_entry.id] = sex_entry.sex

   individual_possible_parents = dict()

   # 1. for every individual keep all 1st degree relatives that may be parents
   for id1, relatives in individual_relatives.iteritems():
      individual_possible_parents[id1] = list()
      for id2, relation in relatives.iteritems():
         # Parent-offspring thresholds are from Conomos et. al., 2016
         if relation.degree == 1 and \
                 (relation.prob_ibd0 is not None and relation.prob_ibd0 < 2**(-4.5) or relation.prop_ibs0 is not None and relation.prop_ibs0 < 0.005):
            individual_possible_parents[id1].append(id2)

   # 2. Drop all individuals that can't form trio family
   for id1 in list(individual_possible_parents):
      if len(individual_possible_parents[id1]) < 2:
         del individual_possible_parents[id1]

   # 3. Two individuals are most likely to be biological parents if:
   #   (1) they are not related (degree > 3), and
   #   (2) they have opposite sex. Trio is created if there is only one pair of most likely biological parents.
   for id1, possible_parents in individual_possible_parents.iteritems():
      unrelated_possible_parents = set()
      for i in xrange(0, len(possible_parents)):
         for j in xrange(i + 1, len(possible_parents)):
            id_i = possible_parents[i]
            id_j = possible_parents[j]
            if id_i not in individual_sex or id_j not in individual_sex:
               continue
            if individual_sex[id_i] == individual_sex[id_j]:
               continue
            relation = individual_relatives[id_i].get(id_j, None)
            if relation is None or relation.degree > min_degree:
               unrelated_possible_parents.add(id_i)
               unrelated_possible_parents.add(id_j)
      if len(unrelated_possible_parents) == 2:
         father = None
         mother = None
         for parent in unrelated_possible_parents:
            if individual_sex[parent] == 'M':
               father = parent
            else:
               mother = parent
         yield Trio(id1, father, mother, individual_sex.get(id1, 'NA'))


def write_trios(trio_entries, out_file):
   with open(out_file, 'w') as ofile:
      ofile.write('IID\tFID\tMID\tSEX\n')
      for trio in trio_entries:
         ofile.write('%s\t%s\t%s\t%s\n' % (trio.individual, trio.father, trio.mother, trio.sex))

