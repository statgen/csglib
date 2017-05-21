import rdflib
from rdflib import URIRef
from rdflib.plugins import sparql
from collections import namedtuple
import gzip
import re

PublishedSNVAssociation = namedtuple('PublishedSNVAssociation', ['chrom', 'position', 'rsid', 'trait', 'pvalue', 'pubmed'], verbose = False)

#graph = rdflib.Graph()
#graph.parse('gwas-kb.owl')
#input: graph
def get_known_associations_owl(graph, trait):
   query = sparql.prepareQuery("""
      SELECT ?chrom_name ?bp ?rsid ?trait_name ?pvalue ?pubmed_id
      WHERE {
         ?assoc rdf:type gwas:TraitAssociation .
         ?assoc gwas:has_gwas_trait_name ?trait_name .
         ?assoc ro:part_of ?study .
         ?study gwas:has_pubmed_id ?pubmed_id .
         ?assoc oban:has_subject ?snp .
         ?assoc gwas:has_p_value ?pvalue .
         ?snp gwas:has_snp_reference_id ?rsid .
         ?snp gwas:has_basepair_position ?bp .
         ?snp ro:located_in ?cyto_region .
         ?cyto_region ro:part_of ?chrom .
         ?chrom gwas:has_name ?chrom_name .
         FILTER (?trait_name = ?trait) }
      """,
      initNs = {
         'gwas' : "http://rdf.ebi.ac.uk/terms/gwas/",
         'ro': "http://www.obofoundry.org/ro/ro.owl#",
         'oban': "http://purl.org/oban/"
      })

   result = graph.query(query, initBindings = {'trait' : rdflib.term.Literal(trait)})

   for row in result:
      yield PublishedSNVAssociation(chrom = row[0].toPython(), position = row[1].toPython(), rsid = row[2].toPython(), trait = row[3].toPython(), pvalue = row[4].toPython(), pubmed = row[5].toPython())


def get_known_snv_associations_tsv(catalog_associations_file, trait = None):
   chr_pos_pattern = re.compile('^(?:chr)?([XY0-9]+):([0-9]+)$', re.IGNORECASE)
   with open(catalog_associations_file, 'r') as ifile:
      line = ifile.readline()
      if not line:
         return
      header = line.rstrip().split('\t')
      for line in ifile:
         fields = dict(zip(header, line.rstrip().split('\t')))
         pvalue = None
         pubmed = None
         study_trait = None

         if not fields['DISEASE/TRAIT']:
            continue
         study_trait = fields['DISEASE/TRAIT']

         if trait and study_trait != trait:
            continue

         if not fields['P-VALUE']:
            continue
         pvalue = fields['P-VALUE']

         if not fields['PUBMEDID']:
            continue
         pubmed = fields['PUBMEDID']

         if fields['SNP_ID_CURRENT']:
            if fields['SNP_ID_CURRENT'].startswith('rs'):
               rsid = fields['SNP_ID_CURRENT']
            else:
               rsid = 'rs' + fields['SNP_ID_CURRENT']
            chrom = None
            position = None
            if fields['CHR_ID']:
               chrom = fields['CHR_ID']
            if fields['CHR_POS']:
               position = long(fields['CHR_POS'])
            if chrom is None or position is None:
               yield PublishedSNVAssociation(chrom = None, position = None, rsid = rsid, trait = study_trait,  pvalue = pvalue, pubmed = pubmed)
               pass
            else:
               yield PublishedSNVAssociation(chrom = chrom, position = position, rsid = rsid, trait = study_trait,  pvalue = pvalue, pubmed = pubmed)
         else:
            m = chr_pos_pattern.match(fields['SNPS'])
            if not m:
               continue
            chrom =  m.group(1)
            position = long(m.group(2))
            yield PublishedSNVAssociation(chrom = chrom, position = position, rsid = None, trait = study_trait,  pvalue = pvalue, pubmed = pubmed)


def get_known_haplotype_associations_tsv(catalog_associations_file, trait = None):
   pass

# dbsnp_file -- VCF file from dbSNP
# rsIds -- set() of rs identifiers that needs to be mapped to chromosome and position
def map2dbsnp(rsids, dbsnp_file):
   rsid2pos = dict()
   with gzip.GzipFile(dbsnp_file, 'r') as ifile:
      for line in ifile:
         if line.startswith('#'):
            if line.startswith('##'):
               continue
            else:
               break
      if any(x != y for (x, y) in zip(line.rstrip().split('\t'), ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])):
         raise Exception('Error while parsing VCF header!')
      for line in ifile:
         fields = line.rstrip().split('\t')
         chrom = fields[0]
         position = long(fields[1])
         rsid = fields[2]

         if rsid in rsids:
            rsid2pos[rsid] = (chrom, position)
   return rsid2pos

# ucsc_snp_file -- from UCSC Table Browser: group = Variaion, track = All SNPs (XYZ), table = snpXYZ
# rsIds --  set() of rs identifiers that needs to be mapped to chromosome and position
def map2ucsc(rsIds, ucsc_snp_file):
   required_header = [
      '#bin',
      'chrom',
      'chromStart',
      'chromEnd',
      'name',
      'score',
      'strand',
      'refNCBI',
      'refUCSC',
      'observed',
      'molType',
      'class',
      'valid',
      'avHet',
      'avHetSE',
      'func',
      'locType',
      'weight',
      'exceptions',
      'submitterCount',
      'submitters',
      'alleleFreqCount',
      'alleles',
      'alleleNs',
      'alleleFreqs',
      'bitfields'
   ]

   chrom_idx = required_header.index('chrom')
   pos_idx = required_header.index('chromStart')
   rsId_idx = required_header.index('name')
   strand_idx = required_header.index('strand')
   alleles_idx = required_header.index('observed')

   with gzip.GzipFile(ucsc_snp_file, 'r') as ifile:
      header = ifile.readline().rstrip().split('\t')

      if len(header) != len(required_header) or any([i != j for i, j in zip(header, required_header)]):
         raise Exception('Header is not valid!')

      for line in ifile:
         line = line.rstrip()
         if not line:
            continue

         fields = line.split('\t')
         try:
            rsId = fields[rsId_idx]
            pos = long(fields[pos_idx]) + 1
            chrom = fields[chrom_idx]
            strand = fields[strand_idx]
            alleles = fields[alleles_idx]
         except IndexError:
            print 'Error while parsing UCSC SNPs file line:', line
            continue

         if not rsId:
            continue

         if rsId in rsIds:
            yield {'CHROM': re.sub('^chr', '', chrom), 'POS': pos , 'STRAND': strand, 'ALLELES': alleles}


#Input: gwas-catalog-ancestry.tsv
def get_study_by_ancestry(catalog_ancestry_file, ancestry):
   with open(catalog_ancestry_file, 'r') as ifile:
      line = ifile.readline()
      if not line:
         return
      header = line.rstrip().split('\t')
      for line in ifile:
         fields = dict(zip(header, line.rstrip().split('\t')))
         if ancestry in fields['BROAD ANCESTRAL CATEGORY']:
            print fields['PUBMEDID'], fields['FIRST AUTHOR'], fields['BROAD ANCESTRAL CATEGORY']


