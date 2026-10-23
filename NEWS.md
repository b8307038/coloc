# coloc 8.0.0
* restructured package
* new `coloc_dataset()` and `as_coloc_dataset()` constructors
* new `wakefield_bf()`, Wakefield's approximate Bayes factor for each snp of a dataset, from beta and varbeta or from p-values; its `effect_priors` applies on both routes
* `coloc()` checks both datasets before fine-mapping either, and its messages about a dataset begin "trait 1: " or "trait 2: "
* `check_dataset()`, and every function that checks a dataset, gives its own message where R's error showed before: for a `beta`, `varbeta` or `pvalues` that isn't numeric, a negative `varbeta`, a `type` or `s` that isn't a single value, and an `LD` without two dimensions. A `coloc_dataset` changed to hold a negative `varbeta`, which gave NA posteriors, is refused
* a dataset is refused when its `N` isn't a finite number above 0, or, on a quantitative trait given as `beta` and `varbeta`, its `sdY` isn't a single finite number above 0. 6.0.3 accepted most of these, and an `sdY` of 0, or an `N` of 0 with `sdY` estimated from it, changed the answer with no sign
* on the SuSiE route, `LD` must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal. An asymmetric LD, which susieR made symmetric with a warning, is refused
* `check_dataset()` names the elements whose lengths differ, as in "snp has 500 values but beta has 499", where it said "lengths of inputs don't match"; a per-snp `N` gets "N must be a single number"
* a dataset may carry `z`, which `coloc_dataset()` takes and `runsusie()` fine-maps from in place of `beta / sqrt(varbeta)`, as it always has. On the SuSiE route a `z` with missing or infinite values, of the wrong length, or not numeric is refused, where it failed with R's errors or was fine-mapped as given
* `check_alignment()` counts each pair of snps once, matches `LD` to the snps by name, and uses `thr`, which it ignored for 0.2. It counted each snp paired with itself, which always agrees, and divided by `LD` by position, so an `LD` in another order than the snps gave the wrong share, and one with more snps failed with R's error. With `thr` at its default, its result for the packaged datasets moves in the fourth decimal place
* new `qc_dataset()` lists every check of a dataset without stopping: the problems `check_dataset()` stops at, and data that looks wrong without being refused, such as an `sdY` far from what `varbeta`, `MAF` and `N` suggest, a `z` that differs from `beta / sqrt(varbeta)`, or z-scores whose signs disagree with `LD`. Printed, it shows a count and only the rows to look at

# coloc 6.0.3
* fix issue #208, #209
# coloc 6.0.0
* Add new ability to calculate varying prior probabilities based on distance to the gene in eQTL datasets
* Add new plot_extended_data function to produce nice QC plots for coloc

# coloc 5.2.3
* minor change, for another package to be able to access one of the internal functions, it was exported

# coloc 5.2.2
* better checking of dataset inputs, and better documentation of what that data should look like

# coloc 5.2.0
* update to new susieR >= 0.12.06 which requires sample size
* general tidy up of susie-related code
* new functions to help find subregions for finemapping - findpeaks, findends
* updated FAQ.md
* tidied empty items in documentation to make CRAN happy
* slimmed down example dataset to make CRAN even happier

# coloc 5.1.2
* deprecate runsusie() argument p

# coloc 5.1.1
* update plot_dataset() and sensitivity() to fix bug for very small p values displaying as 0

# coloc 5.1.0
* release to CRAN of coloc.susie
* remove nref now not recommended

# coloc 5.0.0.9002
* in runsusie, allow nref to be overridden if user passes z_ld_weight

# coloc 5.0.0.9001
* deprecate coloc.signals for multiple causal variants
* introduce coloc.susie for multiple causal variants
* new vignette 02_data giving more details on structuring your data properly

# coloc 4.0-5
* update error checking

# coloc 4.0-3
* BUGFIX
- missing sdY for type="quant" would error in finemap.signals or coloc.signals, now it is estimated if missing.
- warning for factors used as snp names added (issue #29)

# coloc 4.0-2
* BUGFIX
- snps in the second dataset might not have been masked as intended. This would result in effect using the "single" option for trait 2 when masking was expected.

# coloc 4.0-0
*new functions: coloc.signals, finemap.abf*
- analogues of coloc.abf and finemap.abf that allow for multiple causal variants.  See vignette on conditioning/masking
*new function: sensitivity*
- post-hoc, determine the sensitivity of coloc results to changes in the prior.  See vignette on sensitivity

2019-06-25

# coloc 3.2-1
* BUGFIX: finemap.abf()
- in low power situations, the posterior for H0 was previously too low.  This will only affect datasets where the minimum p value was > 1e-7 - ie where the posterior for H0 would be expected to be much above 0.

# coloc 3.2 
- fix bug in process.datasets which suggested MAF was needed for cc data when beta/varbeta also present
- added pkgdown
	
2018-09-28

# coloc 3.1
*new function: finemap.abf*
- improved clarity of error messages in coloc.abf() sub functions
- added finemap.abf to fine map a single trait
- fixed warnings caused by CRAN disliking the BioConductor devel branch
	
2018-02-22
	
# coloc 2.4 
- added stratification to proportional testing approaches
	
2017-03-21 

# coloc 2.3-2 
- tidied code relating mainly to proportional
- colocalisation testing methods, making more methods confirm to S4.
- pcs.prepare now imputes missing genotypes by default
2013-05-22 

# coloc 2.3 
## coloc BUGFIX
- Introduced a function to estimate trait variance from supplied coefficients and standard errors.  This is used within the approach implemented in coloc.abf(), and replaces the earlier version which implicitly assumed that var(Y)=1 for quantitative traits, which could lead to incorrect inference when var(Y) was far from 1.

2013-09-25  

# coloc 2.2
- Merged coloc.abf and coloc.abf.imputed(), so that datasets for which beta, var(beta) are available can be matched to datasets with only p values and maf.2 This means the arguments to coloc.abf() have been changed!  Please
check ?coloc.abf for the new function.

2013-19-06  

# coloc 2.1
- Bug fix for coloc.abf() function, which used p12 instead of log(p12) to calculate L4.
- New function coloc.abf.imputed() to make better use of fuller information on imputed data.

2013-03-06  

# coloc 2.0
- New function, coloc.abf(), to implement the colocalisation approach described by Giambartolomei et. al.
- Changes in the coloc.test() and coloc.bma() functions to make them consistent with regards arguments and output.

2013-22-05  

# coloc 1.11
- added principal components functions pcs.prepare(), pcs.model().
- Restructed the coloc objects to separate Bayesian and non-Bayesian inference.

2012-12-10  

# coloc 1.10
- added Credible Interval calculation to coloc.test().

2012-11-28  

# coloc 1.09
- updated to return u and Var(u) in addition to chisq statistic.

2012-07-12  

# coloc 1.08
- fixed error in documentation, added MASS to Depends.

2012-04-12  

# coloc 1.06
- some tidying.

2012-01-05  

# coloc 1.05
- moved to S4 methods.

2012-01-03  

# coloc 1.04
- made the means of generating plots more flexible.

2011-12-27  

