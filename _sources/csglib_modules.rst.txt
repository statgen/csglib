csglib's modules
================

.. _interval_tree:

IntervalTree
-------------

This module implements interval tree data structure for time efficient (log complexity) interval queries.
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

Please, refer to :doc:`API documentation <csg.intervaltree>` for further details.

.. _pyld:

LD
--

This module is for computing linkage-disequilibrium (LD) (r, r^2) from VCFs.
No special file format is needed - just provide phased VCFs (e.g. 1000G) compressed using bgzip and indexed using tabix. 
Supported operations:

* compute pairwise LD between pair of SNPs from any chromosomal region in VCF;
* compute pairwise LD between all SNPs from a specified region in VCF;
* provides additional routines to compute allele frequencies.

Please, refer to :doc:`API documentation <csg.genetics.ld>` for further details.
