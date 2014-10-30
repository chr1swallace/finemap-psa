ducking-ironman
===============

simple fine mapping script in R

Given a number of index SNPs, and complete summary statistics for a dense GWAS (in this case, from ImmunoChip), we can calculate Bayes Factors through Wakefield's approximations [1] and from these we can generate posterior probabilities that a single SNP is causal in the region given two key assumptions:
- all SNPs have been genotyped
- there is only a single causal variant
Using these, we can generate credible sets for the causal variant - the smallest set of SNPs that capture a given cumulative posterior probability for a region.

The first assumption is rarely met, but the credible sets can still be taken to represent a "best guess" given existing data.

The second assumption may be met, perhaps often is, but when it is violated *and* multiple causal variants are in some LD, this simple method may produce completely spurious results.  Caveat emptor!

[1] Jon Wakeifled (2008) Bayes factors for genome-wide association studies: comparison with P-values. Genet Epidemiol DOI: 10.1002/gepi.20359
