import gzip
import ImputationQuality as imp
import itertools
import Panel
from collections import OrderedDict


headers = {
        'frequentist1': ['alternate_ids', 'rsid', 'chromosome', 'position', 'alleleA', 'alleleB', 'index', 'average_maximum_posterior_call', 'info',
            'cohort_1_AA', 'cohort_1_AB', 'cohort_1_BB', 'cohort_1_NULL', 'all_AA', 'all_AB', 'all_BB', 'all_NULL', 'all_total', 'cases_AA', 'cases_AB', 'cases_BB', 'cases_NULL', 'cases_total',
            'controls_AA', 'controls_AB', 'controls_BB', 'controls_NULL', 'controls_total', 'all_maf', 'cases_maf', 'controls_maf', 'missing_data_proportion',
            'het_OR', 'het_OR_lower', 'het_OR_upper', 'hom_OR', 'hom_OR_lower', 'hom_OR_upper', 'all_OR', 'all_OR_lower', 'all_OR_upper',
            'frequentist_add_pvalue', 'frequentist_add_info', 'frequentist_add_beta_1', 'frequentist_add_se_1', 'comment']
        }


counter_labels = ['n_variants', 'n_removed_variants', 'n_remained_variants',
         'duplicated', 'not_in_panel', 'allele_mismatch', 'pvalue_is_na', 'pvalue_le_zero', 'pvalue_ge_one',
         'info_is_na', 'info_lt_min', 'se_is_na', 'se_le_zero', 'se_gt_max', 'mac_lt_min']


def duplicates(in_file):
   with gzip.GzipFile(in_file, 'r') as ifile:
      for line in ifile:
         if not line.startswith('#'): # skip comments
            break
      if not line:
         return
      header = line.rstrip().split(' ')
      if any(x not in header for x in headers['frequentist1']): # less strict header requirement than the one below
      #if len(snptest_headers['frequentist1']) != len(header) or any(x != y for x, y in zip(snptest_headers['frequentist1'], header)):
         raise Exception('File header does not match expected format!')

      chrom_idx = header.index('chromosome')
      alternateid_idx = header.index('alternate_ids')
      position_idx = header.index('position')
      alleleA_idx = header.index('alleleA')
      alleleB_idx = header.index('alleleB')

      variants = dict()
      n_variants = 0
      for line in ifile:
         if line.startswith('#'): # skip if comments
            continue
         if line.startswith(header[0]): # skip repeated header (left after merging)
            continue

         n_variants += 1

         fields = line.rstrip().split(' ')

         chrom = fields[chrom_idx]
         if chrom == 'NA' or chrom == 'na': # if chromosome fields is NA, then try to extract chromosome name from alternate_id fields assuming CHR:POS... or CHR_POS... format.
            chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[alternateid_idx])).lstrip('0')

         name = chrom + ':' + fields[position_idx] + ':' + fields[alleleA_idx] + ':' + fields[alleleB_idx]

         if name not in variants:
            variants[name] = [n_variants]
         else:
            variants[name].append(n_variants)

      duplicates = set()
      for variant, locations in variants.iteritems():
         if len(locations) > 1:
            for location in locations:
               duplicates.add(location)
   return duplicates


def harmonize_frequentist1(in_file, panel_vcf, imputation_file, min_info, min_mac, max_se, out_file):
   impq = None
   if imputation_file:
      impq = imp.get_imputation_quality(imputation_file)
   duplicated = duplicates(in_file)
   panel = Panel.Panel(panel_vcf, 100000, 1000000)
   with gzip.GzipFile(in_file, 'r') as ifile, gzip.GzipFile(out_file, 'w') as ofile:
      for line in ifile:
         if not line.startswith('#'): # skip comments
            break
      if not line:
         return
      header = line.rstrip().split(' ')
      if any(x not in header for x in headers['frequentist1']): # less strict header requirement than the one below
      #if len(snptest_headers['frequentist1']) != len(header) or any(x != y for x, y in zip(snptest_headers['frequentist1'], header)):
          raise Exception('File header does not match expected format!')

      chrom_idx = header.index('chromosome')
      alternateid_idx = header.index('alternate_ids')
      position_idx = header.index('position')
      alleleA_idx = header.index('alleleA')
      alleleB_idx = header.index('alleleB')
      info_idx = header.index('info')
      pvalue_idx = header.index('frequentist_add_pvalue')
      effect_idx = header.index('frequentist_add_beta_1')
      se_idx = header.index('frequentist_add_se_1')
      total_idx = header.index('all_total')
      cases_total_idx = header.index('cases_total')
      controls_total_idx = header.index('controls_total')
      maf_idx = header.index('all_maf')
      all_ab_idx = header.index('all_AB')
      all_bb_idx = header.index('all_BB')

      ofile.write('UNIQUE_ID\tCHR\tPOSITION\tREF_ALLELE\tALT_ALLELE\tCODED_ALLELE\tNONCODED_ALLELE\tTOTAL\tCASES\tCONTROLS\tAF\tMAF\tEFFECT\tSE\tPVALUE\tINFO\n')

      counters = OrderedDict(zip(counter_labels, [0] * len(counter_labels)))

      for line in ifile:
         if line.startswith('#'): # skip if comments
            continue
         if line.startswith(header[0]): # skip repeated header (left after merging)
            continue

         counters['n_variants'] += 1

         if counters['n_variants'] in duplicated:
            counters['duplicated'] += 1
            counters['n_removed_variants'] += 1
            continue

         fields = line.rstrip().split(' ')

         chrom = fields[chrom_idx].lstrip('0')
         if chrom == 'NA' or chrom == 'na': # if chromosome fields is NA, then try to extract chromosome name from alternate_id fields assuming CHR:POS... or CHR_POS... format.
            chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[alternateid_idx])).lstrip('0')
         position = long(fields[position_idx])

         reference_alleles = panel.get_alleles(chrom, position)
         if reference_alleles is None:
            counters['not_in_panel'] += 1
            counters['n_removed_variants'] += 1
            continue

         noncoded_allele = fields[alleleA_idx]
         coded_allele = fields[alleleB_idx]
         ref_allele = None
         alt_allele = None
         for alleles in reference_alleles:
            if noncoded_allele == alleles.ref and coded_allele in alleles.alt:
               ref_allele = noncoded_allele
               alt_allele = coded_allele
               break
            elif noncoded_allele in alleles.alt and coded_allele == alleles.ref:
               ref_allele = coded_allele
               alt_allele = noncoded_allele
               break
         if ref_allele is None or alt_allele is None:
            counters['allele_mismatch'] += 1
            counters['n_removed_variants'] += 1
            continue

         if impq:
            info_str = impq.get(chrom + ':' + fields[position_idx] + ':' + ref_allele + ':' + alt_allele, 'NA')
         else:
            info_str = fields[info_idx]

         if info_str == 'NA':
            counters['info_is_na'] += 1
            counters['n_removed_variants'] += 1
            continue

         info = float(info_str)
         if info < min_info:
            counters['info_lt_min'] += 1
            counters['n_removed_variants'] += 1
            continue

         total = int(fields[total_idx])
         maf = float(fields[maf_idx])
         mac = 2 * total * maf
         if mac < min_mac:
            counters['mac_lt_min'] += 1
            counters['n_removed_variants'] += 1
            continue
         af = (float(fields[all_ab_idx]) + 2.0 * float(fields[all_bb_idx])) / (2.0 * total)

         if fields[pvalue_idx] == 'NA':
            counters['pvalue_is_na'] += 1
            counters['n_removed_variants'] += 1
            continue
         else:
            pvalue = float(fields[pvalue_idx])
            if pvalue <= 0.0:
               counters['pvalue_le_zero'] += 1
               counters['n_removed_variants'] += 1
               continue
            elif pvalue >= 1.0:
               counters['pvalue_ge_one'] += 1
               counters['n_removed_variants'] += 1
               continue

         if fields[se_idx] == 'NA':
            counters['se_is_na'] += 1
            counters['n_removed_variants'] += 1
            continue
         else:
            se = float(fields[se_idx])
            if se <= 0:
               counters['se_le_zero'] += 1
               counters['n_removed_variants'] += 1
               continue
            elif se > max_se:
               counters['se_gt_max'] += 1
               counters['n_removed_variants'] += 1
               continue

         counters['n_remained_variants'] += 1

         ofile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s\t%s\n' % (
            chrom + ':' + fields[position_idx] + ':' + ref_allele + ':' + alt_allele,
            chrom,
            fields[position_idx],
            ref_allele,
            alt_allele,
            coded_allele,
            noncoded_allele,
            fields[total_idx],
            fields[cases_total_idx],
            fields[controls_total_idx],
            af,
            fields[maf_idx],
            fields[effect_idx],
            fields[se_idx],
            fields[pvalue_idx],
            info_str
          ))

      return counters

