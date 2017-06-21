csglib's modules
================

.. _interval_tree:

Interval tree
-------------

The *IntervalTree* module implements interval tree data structure for time efficient (log complexity) interval queries.
The implementation is based on red-black binary tree and supports the following queries:

* find all intervals overlapping given point;
* find all intervals overlapping given interval;
* find nearest interval to the left of the given point;
* find nearest interval to the right of the given point;
* find *K* nearest intervals to the left of the given point;
* find *K* nearest intervals to the right of the given point;
* get *K* first (leftmost) intervals;
* get *K* last (rightmost) intervals;
* traverse all intervals in ascending/descending order;
* merge overlapping intervals;
* construct complementary intervals (i.e. extracts all gaps between non-overlapping intervals).

Example:
::

   from csg.intervaltree.IntervalTree import IntervalTree

   intervals = IntervalTree() # create new interval tree

   for start, end in zip(range(1, 10), range(5, 14)): # add intervals
      intervals.add(start, end)

   print 'Number of intervals:', intervals.get_intervals_count()
   
   for interval in intervals.descending(): # list intervals in descending order
      print interval.start, interval.end
   
   for interval in intervals.point_intersect(10): # find all intervals that intersect position 10
      print interval.start, interval.end
   
   for interval in intervals.interval_overlap(2, 5): # find all intervals that overlap interval [2, 5]
      print interval.start, interval.end
   
   merged_intervals = intervals.merge() # merge all overlapping intervals
   print 'Number of intervals after merging:', merged_intervals.get_intervals_count()

Please, refer to :doc:`API documentation <csg.intervaltree>` for further details.

.. _pyld:

Linkage disequilibrium
----------------------

The *pyld* module implements computation of *r* coefficient of linkage disequilibrium (LD).
No special file format is required. Phased genotypes must be provided in VCF files compressed using bgzip and indexed using tabix. 
Supported operations:

* compute pairwise LD between pair of SNPs from any given chromosomal region;
* compute pairwise LD between a single SNP and all other SNPs from a given region;
* compute pairwise LD between all SNPs from a given region;
* compute alternate allele frequencies.

Example:
::

   from csg.genetics.ld.pyld import LD

   ld = LD()

   ld.add_vcf('genotypes.phased.vcf.gz') # open VCF with all chromosomes
   # Alternatively, you may load VCF files by chromosome:
   # ld.add_vcf('chr1.vcf.gz') 
   # ld.add_vcf('chr2.vcf.gz')
   # ...
   # ld.add_vcf('chr22.vcf.gz')

   haplotypes = ld.get_region_haplotypes('20', 11650214, 60759931) # read phased genotypes in 20:11650214-60759931 
   
   freqs = ld.compute_freq(haplotypes) # compute alternative allele frequencies
   for i in xrange(0, haploypes.size):
      print haplotypes.chrom[i], haplotypes.position[i], freqs[i]

   r = ld.compute_r_pairwise(haplotypes) # compute LD between all variants in 20:11650214-60759931
   for i in xrange(0, haplotypes.size):
      for j in xrange(i, haplotypes.size):
         print haplotypes.chrom[i], haplotypes.position[i], haplotypes.position[2], r[i, j] ** 2

   haplotypes1 = ld.get_variant_haplotypes('20', 11650214)
   r = ld.compute_r_cross(haplotypes1, haplotypes) # compute LD between variant 20:11650214 and all variants in 20:11650214-60759931
   for i in xrange(0, haplotypes1.size):
      for j in xrange(0, haplotypes.size):
         print haplotypes1.chrom[i], haplotypes1.position[i], haplotypes.position[j], r[i, j] ** 2

   haplotypes2 = ld.get_variant_haplotypes('20', 16655993)
   r = ld.compute_r_cross(haplotypes1, haplotypes2) # compute LD between vatiants 20:11650214 and 20:16655993
   print haplotypes1.chrom[0], haplotypes1.position[0], haplotypes2.position[0], r[0, 0] ** 2

   ld.release_vcfs() # close VCF files


Please, refer to :doc:`API documentation <csg.genetics.ld>` for further details.
