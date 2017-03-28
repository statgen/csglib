# Development

An easy way to install this while developing is: 

```
pip install --user -e .
```

# Structure

### csg / genetics / ld /
Module pyld for computing linkage-disequilibrium (LD) (r, r^2) from VCFs. No special file format is needed - just provide phased VCFs (e.g. 1000G) compressed using bgzip and indexed using tabix. Supported operations: (1) compute pairwise LD between pair of SNPs from any chromosomal region in VCF; (2) compute pairwise LD between all SNPs from a specified region in VCF; Provides additional routines to compute allele frequencies.
