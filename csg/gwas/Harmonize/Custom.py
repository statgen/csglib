import gzip
import itertools
from collections import OrderedDict
import Panel


counter_labels = [
        'n_variants', 'n_removed_variants', 'n_remained_variants',
        'duplicated', 'not_in_panel', 'allele_mismatch', 'pvalue_is_na', 'pvalue_le_zero', 'pvalue_ge_one',
        'info_is_na', 'info_lt_min', 'se_is_na', 'se_le_zero', 'se_gt_max', 'mac_lt_min'
        ]


def duplicates(in_file, sep, comment, chrom_field, position_field, alleleA_field, alleleB_field):
   with gzip.GzipFile(in_file, 'r') as ifile:
      for line in ifile:
         if comment is not None and line.startswith(comment):
            continue
         header = line.rstrip().split(sep)
         break

      try:
         chrom_idx = header.index(chrom_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % chrom_field)

      try:
         position_idx = header.index(position_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % position_field)

      try:
         alleleA_idx = header.index(alleleA_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % alleleA_field)

      try:
         alleleB_idx = header.index(alleleB_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % alleleB_field)

      variants = dict()
      n_variants = 0
      for line in ifile:
         if comment is not None and line.startswith(comment):
            continue
         if line.startswith(header[0]): # Repeated header in the middle of a file. Probably felt after merging.
            continue

         n_variants += 1

         fields = line.rstrip().split(sep)
         chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[chrom_idx])).lstrip('0')
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


def harmonize(in_file, sep, comment, chrom_field, pos_field, coded_allele_field, noncoded_allele_field, coded_allele_freq_field, effect_field, se_field, pvalue_field, cases, controls, min_mac, max_se, panel_vcf, out_file):
   if sep == 'tab':
      sep_char = '\t'
   elif sep == 'whitespace':
      sep_char = ' '
   elif sep == 'comma':
      sep_char = ','
   elif sep == 'semicolon':
      sep_char = ';'
   else:
      raise Exception('Field separator %s is not supported!' % sep)

   duplicated = duplicates(in_file, sep_char, comment, chrom_field, pos_field, coded_allele_field, noncoded_allele_field)
   panel = Panel.Panel(panel_vcf, 100000, 1000000)
   with gzip.GzipFile(in_file, 'r') as ifile, gzip.GzipFile(out_file, 'w') as ofile:
      for line in ifile:
         if comment is not None and line.startswith(comment):
            continue
         header = line.rstrip().split(sep_char)
         break

      try:
         chrom_idx = header.index(chrom_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % chrom_field)

      try:
         position_idx = header.index(pos_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % pos_field)

      try:
         coded_allele_idx = header.index(coded_allele_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % coded_allele_field)

      try:
         noncoded_allele_idx = header.index(noncoded_allele_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % noncoded_allele_field)

      try:
         effect_idx = header.index(effect_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % effect_field)

      try:
         se_idx = header.index(se_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % se_field)

      try:
         pvalue_idx = header.index(pvalue_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % pvalue_field)

      try:
         eaf_idx = header.index(coded_allele_freq_field)
      except ValueError:
         raise Exception('\'%s\' column was not found!' % coded_allele_freq_field)

      total = cases + controls

      ofile.write('UNIQUE_ID\tCHR\tPOSITION\tREF_ALLELE\tALT_ALLELE\tCODED_ALLELE\tNONCODED_ALLELE\tTOTAL\tCASES\tCONTROLS\tAF\tMAF\tEFFECT\tSE\tPVALUE\tINFO\n')

      counters = OrderedDict(zip(counter_labels, [0] * len(counter_labels)))

      for line in ifile:
         if comment is not None and line.startswith(comment):
            continue
         if line.startswith(header[0]):  # Repeated header in the middle of a file. Probably felt after merging.
            continue

         counters['n_variants'] += 1

         if counters['n_variants'] in duplicated:
            counters['duplicated'] += 1
            counters['n_removed_variants'] += 1
            continue

         fields = filter(bool, line.rstrip().split(sep_char))
         chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[chrom_idx])).lstrip('0')
         position = long(fields[position_idx])

         panel_alleles = panel.get_alleles(chrom, position)
         if panel_alleles is None:
            counters['not_in_panel'] += 1
            counters['n_removed_variants'] += 1
            continue

         coded_allele = fields[coded_allele_idx]
         noncoded_allele = fields[noncoded_allele_idx]

         ref_allele = None
         alt_allele = None
         for alleles in panel_alleles:
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

         eaf = float(fields[eaf_idx])
         maf = 1.0 - eaf if eaf > 0.5 else eaf
         mac = 2.0 * total * maf
         if mac < min_mac:
            counters['mac_lt_min'] += 1
            counters['n_removed_variants'] += 1
            continue

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

         ofile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%s\t%g\t%s\t%s\t%s\t%s\n' % (
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
               fields[eaf_idx],
               maf,
               fields[effect_idx],
               fields[se_idx],
               fields[pvalue_idx],
               'NA'
            ))

      return counters
