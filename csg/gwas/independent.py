from collections import namedtuple
from csg.genetics.ld.pyld import LD
import math


Association = namedtuple('Association', ['chrom', 'position', 'ref', 'alt', 'pvalue'], verbose = False)


def get_independent_associations(associations, pvalue_threshold, rsq_threshold, max_locus_bp, in_VCFs):
   associations_by_chrom = dict()
   for association in associations:
      if association.pvalue > pvalue_threshold:
         continue
      if association.chrom not in associations_by_chrom:
         associations_by_chrom[association.chrom] = list()
      associations_by_chrom[association.chrom].append(association)

   ld = LD()
   for in_VCF in in_VCFs:
      ld.add_vcf(in_VCF)

   for chrom, associations in associations_by_chrom.iteritems():
      associations.sort(key = lambda association: (association.pvalue, association.position))
      while associations:
         top_association = associations.pop(0)
         yield top_association
         top_haplotype = ld.get_variant_haplotypes_strict(top_association.chrom, top_association.position, top_association.ref, top_association.alt)
         if top_haplotype is None:
            continue
         to_remove = set()
         for i, association in enumerate(associations):
            if abs(top_association.position - association.position) >= max_locus_bp / 2.0:
               continue
            haplotype = ld.get_variant_haplotypes_strict(association.chrom, association.position, association.ref, association.alt)
            if haplotype is None:
               continue
            r = ld.compute_r(top_haplotype, haplotype)
            if math.isnan(r):
               continue
            if r ** 2 >= rsq_threshold:
               to_remove.add(i)
         associations[:] = (x for i, x in enumerate(associations) if i not in to_remove)

   ld.release_vcfs()

