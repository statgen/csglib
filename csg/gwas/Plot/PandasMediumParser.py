import re
import gzip

class PandasMediumParser(object):

   def __init__(self, file_name, sep, chrom_field_regex, position_field_regex, pvalue_field):
      self.file_name = file_name
      self.sep = sep
      self.chrom_field_regex = chrom_field_regex
      self.position_field_regex = position_field_regex
      self.pvalue_field = pvalue_field
      self.it = iter(self)


   def initialize(self):
      self.sep_char = self.get_separator_character(self.sep)
      self.chrom_field, self.chrom_regex = self.parse_regex_argument(self.chrom_field_regex)
      self.position_field, self.position_regex = self.parse_regex_argument(self.position_field_regex)
      self.ifile = gzip.GzipFile(self.file_name, 'r')
      line = self.ifile.readline()
      if line is None:
         return
      header = line.rstrip().split(self.sep_char)
      indices = self.get_fields_indices(header, self.chrom_field, self.position_field, self.pvalue_field)
      self.chrom_idx = indices[self.chrom_field]
      self.position_idx = indices[self.position_field]
      self.pvalue_idx = indices[self.pvalue_field]


   def __iter__(self):
      self.initialize()
      return self


   def next(self):
      line = self.ifile.readline()
      if not line:
         self.ifile.close()
         raise StopIteration()

      fields = line.rstrip().split(self.sep_char)

      if self.chrom_regex is not None:
         match = self.chrom_regex.search(fields[self.chrom_idx])
         if match is None or not match.group(1):
            raise Exception('Error while extracting chromosome name form %s!' % fields[self.chrom_idx])
         chrom = match.group(1)
      else:
         chrom = fields[self.chrom_idx]

      if self.position_regex is not None:
         match = self.position_regex.search(fields[self.position_idx])
         if match is None or not match.group(1):
            raise Exception('Error while extracting chromosome name form %s!' % fields[self.position_idx])
         position = match.group(1)
      else:
         position = fields[self.position_idx]

      pvalue = fields[self.pvalue_idx]
      return chrom + '\t' + position + '\t' + pvalue + '\n'


   def read(self, n = 0):
      try:
         return next(self.it)
      except StopIteration:
         return ''


   def get_separator_character(self, name):
      if name == 'tab':
         return '\t'
      elif name == 'whitespace':
         return ' '
      elif name == 'comma':
         return ','
      elif name == 'semicolon':
         return ';'
      else:
         raise Exception('Field separator %s is not supported!' % name)


   def parse_regex_argument(self, field_regex):
      field_regex = field_regex.split(',')
      field = field_regex[0]
      regex = None
      if len(field_regex) > 1:
         regex = re.compile(field_regex[1])
      return (field, regex)


   def get_fields_indices(self, header, *field_names):
      field_indices = dict()
      for field_name in field_names:
         try:
            index = header.index(field_name)
         except ValueError:
            raise Exception('Field \'%s\' was not found in input file!' % field_name)
         field_indices[field_name] = index
      return field_indices

