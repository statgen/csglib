import gzip
import itertools


def vcf_entries_it(entries):
   for entry in entries:
      if entry.startswith('#'):
         continue
      fields = entry.rstrip().split('\t', 8)
      chrom = fields[0]
      position = fields[1]
      ref = fields[3]
      alt = fields[4]
      name = chrom + ':' + position + ':' + ref + ':' + alt
      for key_value in fields[7].split(';'):
         if key_value.startswith('R2='):
            yield (name, key_value.split('=')[1])
            break


def snpstat_entries_it(entries, variantid_idx, chrom_idx, position_idx, ref_idx, alt_idx, quality_idx):
   for entry in entries:
      fields = entry.rstrip().split('\t')
      chrom = fields[chrom_idx]
      # if chromosome field is NA, then try to extract chromosome name from SNPID fields assuming CHR:POS... or CHR_POS... format.
      if chrom == 'NA' or chrom == 'na' or chrom == '.':
         chrom = "".join(itertools.takewhile(lambda x: x != ':' and x != '_', fields[variantid_idx])).lstrip('0')
      name = chrom + ':' + fields[position_idx] + ':' + fields[ref_idx] + ':' + fields[alt_idx]
      yield (name, fields[quality_idx])


def get_imputation_quality(in_file):
   vcf_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
   snpstat_header = ['SNPID', 'RSID', 'chromosome', 'position', 'A_allele', 'B_allele', 'minor_allele', 'major_allele',
           'AA', 'AB', 'BB', 'AA_calls', 'AB_calls', 'BB_calls', 'MAF', 'HWE', 'missing', 'missing_calls', 'information']
   unique_variants = dict()
   with gzip.GzipFile(in_file, 'r') as ifile:
      for line in ifile:
         if not line.startswith('##'): #skip all possible VCF-like comments
            break
      if line.startswith('\t'.join(vcf_header)): # VCF
         variants = vcf_entries_it(ifile)
      elif line.startswith('\t'.join(snpstat_header)): # snpstat
         variants = snpstat_entries_it(ifile,
                 snpstat_header.index('SNPID'), snpstat_header.index('chromosome'), snpstat_header.index('position'),
                 snpstat_header.index('A_allele'), snpstat_header.index('B_allele'), snpstat_header.index('information'))
      else:
         raise Exception('Provided file format with imputation qualities is not supported!')
      for variant in variants:
         if variant[0] in unique_variants:
            raise Exception('Duplicated variant in imputation file: %s!' % variant[0])
         unique_variants[variant[0]] = variant[1]
   return unique_variants

