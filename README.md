# Development

An easy way to install this while developing is: 

```
pip install --user -e .
```

# Structure

### csg / intervaltree /
Module intervaltree implements interval tree data structure for time efficient (log complexity) interval queries. The implementation is based on red-black binary tree and supports the following queries: 
1. find all intervals overlapping the given point; 
2. find all intervals overlapping the given interval; 
3. find nearest interval from the left of the given point; 
4. find nearest interval from the right of the given point; 
5. find K nearest intervals from the left of the given point; 
6. find K nearest intervals from the right of the given point; 
7. get first K intervals; 
8. get last K intervals.

### csg / genetics / ld /
Module pyld for computing linkage-disequilibrium (LD) (r, r^2) from VCFs. No special file format is needed - just provide phased VCFs (e.g. 1000G) compressed using bgzip and indexed using tabix. Supported operations: 
1. compute pairwise LD between pair of SNPs from any chromosomal region in VCF;
2. compute pairwise LD between all SNPs from a specified region in VCF;
3. provides additional routines to compute allele frequencies.

### csg / gwas /
Module gwas has miscellaneous functions for manipulating genome-wide associations results such as:
1. independent.py -- get independent most significant association hits based on linkage-disequilibrium

### csg / pedigree / trios / 
Module trios for constructing trio families (two parents and a child) based on genetic kinship estimates and sex. Includes utilities to read PC-relate and KING software outputs. 
