import rdflib
from rdflib import URIRef
from rdflib.plugins import sparql
from collections import namedtuple
import gzip
import re

PublishedAssociation = namedtuple('PublishedAssociation', ['chrom', 'position', 'rsid', 'trait', 'pvalue', 'pubmed'], verbose = False)

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
      yield PublishedAssociation(chrom = row[0].toPython(), position = row[1].toPython(), rsid = row[2].toPython(), trait = row[3].toPython(), pvalue = row[4].toPython(), pubmed = row[5].toPython())


def get_known_associations_tsv(catalog_associations_file, trait):
   with open(catalog_associations_file, 'r') as ifile:
      line = ifile.readline()
      if not line:
         return
      header = line.rstrip().split('\t')
      for line in ifile:
         fields = dict(zip(header, line.rstrip().split('\t')))
         if not fields['SNP_ID_CURRENT']:
            #print fields['SNPS'], fields['MERGED']
            continue
         if fields['DISEASE/TRAIT'] != trait:
            continue
         yield PublishedAssociation(
            chrom = fields['CHR_ID'],
            position = fields['CHR_POS'],
            rsid = fields['SNP_ID_CURRENT'],
            trait = fields['DISEASE/TRAIT'],
            pvalue = fields['P-VALUE'],
            pubmed = fields['PUBMEDID']
         )

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


