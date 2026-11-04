# coloc 8.0.0
* restructured package
* new `coloc_dataset()` and `as_coloc_dataset()` constructors
* new `wakefield_bf()`, Wakefield's approximate Bayes factor for each snp of a dataset, from beta and varbeta or from p-values; its `effect_priors` applies on both routes
* `coloc()` checks both datasets before fine-mapping either, and the checks' messages about a dataset begin "dataset 1: " or "dataset 2: ", as in the older functions. The checks' errors and warnings show no call, so the older functions print "Error: dataset 1: ..." where 6.0.3 printed "Error in check_dataset(d = dataset1, 1) : dataset 1: ..."
* `check_dataset()` puts no prefix on its messages when called without a label, where 6.0.3 began them "dataset : ". With a label, as in `check_dataset(d, 2)`, every message begins "dataset 2: ", the minimum p-value warning and the messages about `N` and `LD` included
* `coloc.susie_bf()` refuses a `bf2` that isn't numeric, such as a dataset or a SuSiE fit, saying what it was given and what `bf2` must hold, where 6.0.3 failed inside `coloc.bf_bf()` with R's errors
* `check_dataset()`, and every function that checks a dataset, gives its own message where R's error showed before: for a `beta`, `varbeta` or `pvalues` that isn't numeric, or a negative `varbeta`, where each is read; for a `type` or `s` that isn't a single value; and for an `LD` with fewer than two dimensions. The exception is a logical `beta` or `varbeta`, which 6.0.3 computed with and which is now refused too. An `LD` with more than two dimensions is refused wherever it is given, though 6.0.3 took one where nothing reads `LD`. A `coloc_dataset` changed to hold a negative `varbeta`, which gave NA posteriors, is refused
* `N` must be a finite number above 0 wherever it is given, and above 1 on the SuSiE route, as susieR needs, so a dataset with a smaller one is refused before either trait is fine-mapped. `sdY` must be a single finite number above 0 where Wakefield's method uses it, for a quantitative trait given as `beta` and `varbeta`; SuSiE doesn't use it. 6.0.3 accepted most of these. An `sdY` of 0, or an `N` of 0 with `sdY` estimated from it, changed the answer with no sign; a negative `sdY`, or an `N` of 0 beside `sdY`, gave the right answer and is refused too
* on the SuSiE route, `LD` must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal. An asymmetric LD, which susieR made symmetric with a message beginning "WARNING:", is refused
* on the SuSiE route, an `LD` that isn't a numeric matrix, such as a data frame or a sparse `Matrix`, gets "LD must be a numeric matrix" from the checks, where `runsusie()` said "LD must be of class matrix"
* `check_dataset()` names the elements whose lengths differ, as in "snp has 500 values but beta has 499", where it said "lengths of inputs don't match"; a per-snp `N` gets "N must be a single number"
* a dataset may carry `z`, which `coloc_dataset()` takes but doesn't compute. When a dataset has a `z`, `runsusie()` fine-maps from it in place of `beta / sqrt(varbeta)`, as 6.0.3 did. On the SuSiE route a `z` with missing or infinite values, of the wrong length, or not numeric is refused, where it failed with R's errors or was fine-mapped as given
* `check_alignment()` counts each pair of snps once, matches `LD` to the snps by name, and uses `thr`, which it ignored for 0.2. It counted each snp paired with itself, which always agrees, and divided by `LD` by position, so an `LD` in another order than the snps gave the wrong share, and one with more snps failed with R's error. With `thr` at its default, `check_alignment()`'s result for the packaged datasets drops by less than 0.0004, as it no longer counts each snp paired with itself, which always agrees
* `check_alignment()` refuses a dataset without `LD`, `beta` or `varbeta`, or whose `LD` isn't a numeric correlation matrix, as on the SuSiE route, where 6.0.3 failed with R's error or, for an `LD` that isn't a correlation matrix, gave a share
* `check_dataset()` warns when a quantitative trait's `sdY` is more than twice or under half its estimate from `varbeta`, `MAF` and `N`, as when `varbeta` holds standard errors. It warns for a `coloc_dataset` too, each time it is checked, and so the old functions warn as well
* `coloc()` now converts a plain list to a `coloc_dataset`, as `wakefield_bf()` already did, and the first time either is given a plain list in a session, it says so in a warning; `check_dataset()`, `runsusie()` and the older functions take plain lists without the warning
* warnings a caller may want to silence on their own have a class of their own, such as `coloc_min_p` and `coloc_sdY`, for `suppressWarnings(classes = )` in R 4.0 and later. coloc now needs R 3.6.0 or later
* the checks' messages say what is wrong and what to do, naming the element and, where it helps, a snp at fault, as in "snp has duplicated ids, such as s1; each snp must appear once"
* `coloc.abf()`, `coloc.bf_bf()`, `coloc.susie()` and `coloc.susie_bf()` refuse a `p12` outside [p1 * p2, min(p1, p2)], as `coloc()` does, where 6.0.3 took it. `coloc.abf()` also refuses a prior weight of 0, as 6.0.3's `coloc.bf_bf()` did, and gives coloc's message for a `p1`, `p2` or `p12` that isn't a single number, where R's error showed
* `coloc.bf_bf()` stops on a missing or `+Inf` log Bayes factor in either matrix, as `coloc()` does, where 6.0.3 failed with R's error or, with `trim_by_posterior = FALSE`, reported a certainty; `-Inf` is still taken
* when every signal pair is left out under `overlap.min`, `PP.H0.abf` is NA, like `PP.H1.abf` to `PP.H4.abf`, as the hypotheses are weighed on the snps the traits share; 6.0.3 reported a number for it. Each row names each trait's own top snp, where 6.0.3 could name the null column
* `check_ld()` refuses an `LD` that isn't a numeric correlation matrix or has missing values, where 6.0.3 computed with it or failed with R's error. `finemap.signals()` and `coloc.signals()` run it on a dataset they fine-map by conditioning, with `method = "cond"`
* `sensitivity()` refuses an `obj` that isn't the result of one of the functions it reads, with coloc's error, where it failed with R's

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

