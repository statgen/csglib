import gzip
import ImputationQuality  as imp
import Panel
import itertools
from collections import OrderedDict


headers = {
        'firth': ['#CHROM', 'BEGIN', 'END', 'MARKER_ID', 'NS', 'AC', 'CALLRATE', 'MAF', 'PVALUE', 'BETA', 'SEBETA', 'CHISQ', 'NS.CASE', 'NS.CTRL', 'AF.CASE', 'AF.CTRL'],
        'emmax': ['#CHROM', 'BEG', 'END', 'MARKER_ID', 'NS', 'AC', 'CALLRATE', 'GENOCNT', 'MAF', 'STAT', 'PVALUE', 'BETA', 'SEBETA', 'R2', 'CTRLCNT', 'CASECNT']
        }


counter_labels = [
        'n_variants', 'n_removed_variants', 'n_remained_variants',
        'duplicated', 'not_in_panel', 'allele_mismatch', 'pvalue_is_na', 'pvalue_le_zero', 'pvalue_gt_one',
        'info_is_na', 'info_lt_min', 'se_is_na', 'se_le_zero', 'se_gt_max', 'mac_lt_min'
        ]


def duplicates(in_file):
   with gzip.GzipFile(in_file, 'r') as ifile:
      line = ifile.readline()
      if not line:
         return
      header = line.rstrip().split('\t')

      if len(headers['firth']) == len(header) and all(x == y for x, y in zip(headers['firth'], header)):
         chrom_idx = header.index('#CHROM')
         position_idx = header.index('BEGIN')
         markerid_idx = header.index('MARKER_ID')
      elif len(headers['emmax']) == len(header) and all(x == y for x, y in zip(headers['emmax'], header)):
         chrom_idx = header.index('#CHROM')
         position_idx = header.index('BEG')
         markerid_idx = header.index('MARKER_ID')
      else:
         raise Exception('File header does not match expected format!')

      variants = dict()
      n_variants = 0
      for line in ifile:
         if line.startswith('#'): # skip if comment or repeated header (left after merging)
            continue
         n_variants += 1
         fields = line.rstrip().split('\t')
         chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[chrom_idx])).lstrip('0')
         markerid = fields[markerid_idx]
         noncoded_allele, coded_allele = markerid.split('_')[1].split('/')
         name = chrom + ':' + fields[position_idx] + ':' + noncoded_allele + ':' + coded_allele

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


def harmonize_firth(in_file, panel_vcf, imputation_file, min_info, min_mac, max_se, out_file):
   impq = imp.get_imputation_quality(imputation_file)
   duplicated = duplicates(in_file)
   panel = Panel.Panel(panel_vcf, 100000, 1000000)

   with gzip.GzipFile(in_file, 'r') as ifile, gzip.GzipFile(out_file, 'w') as ofile:
      line = ifile.readline()
      if not line:
         return
      header = line.rstrip().split('\t')

      if len(headers['firth']) != len(header) or any(x != y for x, y in zip(headers['firth'], header)):
         raise Exception('File header does not match expected format!')

      chrom_idx = header.index('#CHROM')
      position_idx = header.index('BEGIN')
      markerid_idx = header.index('MARKER_ID')
      pvalue_idx = header.index('PVALUE')
      effect_idx = header.index('BETA')
      se_idx = header.index('SEBETA')
      total_idx = header.index('NS')
      cases_total_idx = header.index('NS.CASE')
      controls_total_idx = header.index('NS.CTRL')
      maf_idx = header.index('MAF')
      ac_idx = header.index('AC')

      ofile.write('UNIQUE_ID\tCHR\tPOSITION\tREF_ALLELE\tALT_ALLELE\tCODED_ALLELE\tNONCODED_ALLELE\tTOTAL\tCASES\tCONTROLS\tAF\tMAF\tEFFECT\tSE\tPVALUE\tINFO\n')

      counters = OrderedDict(zip(counter_labels, [0] * len(counter_labels)))

      for line in ifile:
         if line.startswith('#'): # skip if comment or repeated header (left after merging)
            continue

         counters['n_variants'] += 1

         if counters['n_variants'] in duplicated:
            counters['duplicated'] += 1
            counters['n_removed_variants'] += 1
            continue

         fields = line.rstrip().split('\t')
         chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[chrom_idx])).lstrip('0')
         position = long(fields[position_idx])

         reference_alleles = panel.get_alleles(chrom, position)
         if reference_alleles is None:
            counters['not_in_panel'] += 1
            counters['n_removed_variants'] += 1
            continue

         markerid = fields[markerid_idx]
         noncoded_allele, coded_allele = markerid.split('_')[1].split('/')
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

         info_str = impq.get(chrom + ':' + fields[position_idx] + ':' + ref_allele + ':' + alt_allele, 'NA')

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
         af = float(fields[ac_idx]) / (2.0 * total)

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


def harmonize_emmax(in_file, panel_vcf, imputation_file, min_info, min_mac, max_se, out_file):
   impq = imp.get_imputation_quality(imputation_file)
   duplicated = duplicates(in_file)
   panel = Panel.Panel(panel_vcf, 100000, 1000000)

   with gzip.GzipFile(in_file, 'r') as ifile, gzip.GzipFile(out_file, 'w') as ofile:
      line = ifile.readline()
      if not line:
         return
      header = line.rstrip().split('\t')

      if len(headers['emmax']) != len(header) or any(x != y for x, y in zip(headers['emmax'], header)):
         raise Exception('File header does not match expected format!')

      chrom_idx = header.index('#CHROM')
      position_idx = header.index('BEG')
      markerid_idx = header.index('MARKER_ID')
      pvalue_idx = header.index('PVALUE')
      effect_idx = header.index('BETA')
      se_idx = header.index('SEBETA')
      total_idx = header.index('NS')
      cases_total_idx = header.index('CASECNT')
      controls_total_idx = header.index('CTRLCNT')
      maf_idx = header.index('MAF')
      ac_idx = header.index('AC')

      ofile.write('UNIQUE_ID\tCHR\tPOSITION\tREF_ALLELE\tALT_ALLELE\tCODED_ALLELE\tNONCODED_ALLELE\tTOTAL\tCASES\tCONTROLS\tAF\tMAF\tEFFECT\tSE\tPVALUE\tINFO\n')

      counters = OrderedDict(zip(counter_labels, [0] * len(counter_labels)))

      for line in ifile:
         if line.startswith('#'): # skip if comment or repeated header (left after merging)
            continue

         counters['n_variants'] += 1

         if counters['n_variants'] in duplicated:
            counters['duplicated'] += 1
            counters['n_removed_variants'] += 1
            continue

         fields = line.rstrip().split('\t')
         chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[chrom_idx])).lstrip('0')
         position = long(fields[position_idx])

         reference_alleles = panel.get_alleles(chrom, position)
         if reference_alleles is None:
            counters['not_in_panel'] += 1
            counters['n_removed_variants'] += 1
            continue

         markerid = fields[markerid_idx]
         noncoded_allele, coded_allele = markerid.split('_')[1].split('/')
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

         info_str = impq.get(chrom + ':' + fields[position_idx] + ':' + ref_allele + ':' + alt_allele, 'NA')

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
         af = float(fields[ac_idx]) / (2.0 * total)

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

         ofile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%g\t%s\t%s\t%s\t%s\t%s\n' % (
            chrom + ':' + fields[position_idx] + ':' + ref_allele + ':' + alt_allele,
            chrom,
            fields[position_idx],
            ref_allele,
            alt_allele,
            coded_allele,
            noncoded_allele,
            fields[total_idx],
            sum([int(x) for x in fields[cases_total_idx].split('/')]),
            sum([int(x) for x in fields[controls_total_idx].split('/')]),
            af,
            fields[maf_idx],
            fields[effect_idx],
            fields[se_idx],
            fields[pvalue_idx],
            info_str
          ))

      return counters
