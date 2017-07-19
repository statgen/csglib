csglib's tools
==============

.. _harmonize:

Harmonize
---------

The *Harmonize* provides a set of tools for filtering and formatting GWAS files before meta-analysis.
It processes single GWAS input file as follows: 

* applies filters on minor allele count, imputation quality, and standard error;
* renames and reorders field names;
* changes field separator to tabulation;
* adds imputation quality (if qualities were in different input file);
* check alleles against provided reference panel.

EPACTS, SNPTEST, and QuickTest input formats are supported. The custom input file formats maybe processed as well, by providing additional commandline arguments.


.. _transform_effect:

TransformEffect
---------------

The *TransformEffect* tool transforms effect from linear regression model to log-OR scale. See `PMC5237383 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5237383/>`_ for details.

.. _inflation_from_file:

InflationFromFile
-----------------

The *InflationFromFile* tool computes genomic inflation factor (lambda) from provided the GWAS files.

.. _filter:

Filter
------

The *Filter* tool provide the following filtering options for GWAS files:

* filter variants based on the number of missing direction effects after meta-analysis using METAL software;
* filter significant associations based on P-value;
* filters/remove associations falling into the specified region(s).

.. _plot:

Plot
----

The *Plot* provides a set of tools for creating plots from files with GWAS results.
Currently, it supports Manhattan and Q-Q plots.
