import gzip
import ImputationQuality as imp
from collections import OrderedDict
import itertools
import Panel

headers = {
           'default': ['chr', 'rsid', 'pos', 'allele_A', 'allele_B', 'all_AA', 'all_AB', 'all_BB', 'cases_AA', 'cases_AB', 'cases_BB',
           'controls_AA', 'controls_AB', 'controls_BB', 'all_maf', 'cases_maf',	'controls_maf',	'all_hwe', 'cases_hwe',	'controls_hwe',
           'Imputed', 'info', 'beta', 'SE', 'P']
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
      header = line.rstrip().split('\t')
      if any(x not in header for x in headers['default']): # less strict header requirement than the one below
      #if len(quicktest_headers['default']) != len(header) or any(x != y for x, y in zip(quicktest_headers['default'], header)):
         raise Exception('File header does not match expected format!')

      chrom_idx = header.index('chr')
      rsid_idx = header.index('rsid')
      position_idx = header.index('pos')
      alleleA_idx = header.index('allele_A')
      alleleB_idx = header.index('allele_B')

      variants = dict()
      n_variants = 0
      for line in ifile:
         if line.startswith('#'): # skip if comments
            continue
         if line.startswith(header[0]): # skip repeated header (left after merging)
            continue

         n_variants += 1

         fields = line.rstrip().split('\t')

         chrom = fields[chrom_idx]
         if chrom == 'NA' or chrom == 'na': # if chromosome fields is NA, then try to extract chromosome name from rsid fields assuming CHR:POS... or CHR_POS... format.
            chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[rsid_idx])).lstrip('0')

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


def harmonize_default(in_file, panel_vcf, imputation_file, min_info, min_mac, max_se, out_file):
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
      header = line.rstrip().split('\t')
      if any(x not in header for x in headers['default']): # less strict header requirement than the one below
      #if len(quicktest_headers['default']) != len(header) or any(x != y for x, y in zip(quicktest_headers['default'], header)):
          raise Exception('File header does not match expected format!')

      chrom_idx = header.index('chr')
      rsid_idx = header.index('rsid')
      position_idx = header.index('pos')
      alleleA_idx = header.index('allele_A')
      alleleB_idx = header.index('allele_B')
      info_idx = header.index('info')
      pvalue_idx = header.index('P')
      effect_idx = header.index('beta')
      se_idx = header.index('SE')
      all_AA_idx = header.index('all_AA')
      all_AB_idx = header.index('all_AB')
      all_BB_idx = header.index('all_BB')
      cases_AA_idx = header.index('cases_AA')
      cases_AB_idx = header.index('cases_AB')
      cases_BB_idx = header.index('cases_BB')
      controls_AA_idx = header.index('controls_AA')
      controls_AB_idx = header.index('controls_AB')
      controls_BB_idx = header.index('controls_BB')
      maf_idx = header.index('all_maf')

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

         fields = line.rstrip().split('\t')

         chrom = fields[chrom_idx].lstrip('0')
         if chrom == 'NA' or chrom == 'na': # if chromosome fields is NA, then try to extract chromosome name from rsid fields assuming CHR:POS... or CHR_POS... format.
            chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[rsid_idx])).lstrip('0')

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

         total = float(fields[all_AA_idx]) + float(fields[all_AB_idx]) + float(fields[all_BB_idx])
         cases = float(fields[cases_AA_idx]) + float(fields[cases_AB_idx]) + float(fields[cases_BB_idx])
         controls = float(fields[controls_AA_idx]) + float(fields[controls_AB_idx]) + float(fields[controls_BB_idx])
         maf = float(fields[maf_idx])
         mac = 2.0 * total * maf
         if mac < min_mac:
            counters['mac_lt_min'] += 1
            counters['n_removed_variants'] += 1
            continue
         af = (float(fields[all_AB_idx]) + 2.0 * float(fields[all_BB_idx])) / (2.0 * total)

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
            elif pvalue > 1.0:
               counters['pvalue_gt_one'] += 1
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

         ofile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%s\t%s\t%s\t%s\t%s\n' % (
            chrom + ':' + fields[position_idx] + ':' + ref_allele + ':' + alt_allele,
            chrom,
            fields[position_idx],
            ref_allele,
            alt_allele,
            coded_allele,
            noncoded_allele,
            total,
            cases,
            controls,
            af,
            fields[maf_idx],
            fields[effect_idx],
            fields[se_idx],
            fields[pvalue_idx],
            info_str
          ))

      return counters

