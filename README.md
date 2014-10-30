ducking-ironman
===============

simple fine mapping script in R

Given a number of index SNPs, and complete summary statistics for a dense GWAS (in this case, from ImmunoChip), we can calculate Bayes Factors through Wakefield's approximations [1] and from these we can generate posterior probabilities that a single SNP is causal in the region given two key assumptions:
- all SNPs have been genotyped
- there is only a single causal variant
Using these, we can generate credible sets for the causal variant - the smallest set of SNPs that capture a given cumulative posterior probability for a region.

The first assumption is rarely met, but the credible sets can still be taken to represent a "best guess" given existing data.

The second assumption may be met, perhaps often is, but when it is violated *and* multiple causal variants are in some LD, this simple method may produce completely spurious results.  Caveat emptor!

[1] Jon Wakefield (2008) Bayes factors for genome-wide association studies: comparison with P-values. Genet Epidemiol DOI: 10.1002/gepi.20359

Files
-----

- table2-mod.csv is Table 2 from Eyre et al Nat Genet. 2012 Dec;44(12):1336-40. doi: 10.1038/ng.2462.
- bf-functions.R is an R script processing these results

Missing file
------------

These functions only work given a complete set of p values from the RA ImmunoChip study.  These may be downloaded from http://www.immunobase.org
