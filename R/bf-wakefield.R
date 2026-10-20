## Wakefield's approximate Bayes factors, one per snp

##' Wakefield's approximate Bayes factor for each snp
##'
##' `wakefield_bf()` gives Wakefield's approximate Bayes factor (ABF) for each
##' snp in a dataset, on the log scale: the evidence for an effect at that snp
##' against no effect. These are the per-snp Bayes factors [coloc.abf()] and
##' [finemap.abf()] use, here one row per snp with the working kept.
##'
##' The route depends on what the dataset has:
##'
##' - `beta` and `varbeta`: z is `beta / sqrt(varbeta)` and V is `varbeta`.
##'   For a quantitative trait the effect prior is in units of the trait's
##'   standard deviation, so it is multiplied by `sdY`; when the dataset has
##'   no `sdY`, it is estimated from `varbeta`, `MAF` and `N`, with a warning.
##' - `pvalues`, `MAF` and `N`, and `s` for a case-control trait: z comes from
##'   the p-value, and V from `MAF`, `N` and `s`, for a trait with standard
##'   deviation 1.
##'
##' Either way, with W the prior standard deviation squared, r = W / (W + V)
##' and the log ABF is (log(1 - r) + r z^2) / 2. With the default
##' `effect_priors`, the numbers are those of [process.dataset()].
##'
##' The dataset is checked with [check_dataset()] first.
##'
##' @param d a dataset made with [coloc_dataset()], or a list with the same
##'     elements
##' @param effect_priors the prior standard deviation of a causal snp's
##'     effect, named by the type of trait: `quant` in units of the trait's
##'     standard deviation, and `cc` on the log odds ratio scale
##' @return a data frame with one row per snp and the columns `snp`,
##'     `position` (`NA` when the dataset has none), `z`, `V` (the variance
##'     of the effect estimate), `r` (the shrinkage factor) and `lABF` (the
##'     log ABF)
##' @references Wakefield J (2009). Bayes factors for genome-wide association
##'     studies: comparison with P-values. Genetic Epidemiology 33(1),
##'     79-86. \doi{10.1002/gepi.20359}
##' @seealso [finemap.abf()], which turns these into posterior probabilities
##' @author Claudia Giambartolomei, Vincent Plagnol, Chris Wallace
##' @export
##' @examples
##' data(coloc_test_data)
##' head(wakefield_bf(coloc_test_data$D1))
##'
##' ## a wider prior for the effect of a causal snp
##' head(wakefield_bf(coloc_test_data$D1, effect_priors = c(quant = 0.3, cc = 0.2)))
wakefield_bf <- function(d, effect_priors = c(quant = 0.15, cc = 0.2)) {
    if (!is.numeric(effect_priors) ||
        !all(c("quant", "cc") %in% names(effect_priors)) ||
        !all(is.finite(effect_priors[c("quant", "cc")]) &
             effect_priors[c("quant", "cc")] > 0))
        stop("effect_priors should be positive numbers named quant and cc")
    check_dataset(d)
    bf <- .wakefield_bf(d, effect_priors)
    data.frame(snp = as.character(d$snp),
               position = if ("position" %in% names(d)) d$position
                          else rep(NA_integer_, nrow(bf)),
               bf, row.names = NULL)
}

## The arithmetic of wakefield_bf(), unchecked: z, V, r and lABF for each snp,
## from beta and varbeta when the dataset has both, else from p-values, MAF
## and N. process.dataset() calls it directly, as its callers check the
## dataset first, so a dataset is checked, and warned about, once.
.wakefield_bf <- function(d, effect_priors) {
    nd <- names(d)
    ## by the type's label: [[ would take a factor type by its code. A type
    ## other than quant gets cc's prior, as in 6.0.3
    sd_prior <- effect_priors[[if (d$type == "quant") "quant" else "cc"]]
    if ("beta" %in% nd && "varbeta" %in% nd) {
        z <- d$beta / sqrt(d$varbeta)
        V <- d$varbeta
        ## a quantitative trait's prior is in units of its sd, estimated, with
        ## a warning, when not given
        if (d$type == "quant")
            sd_prior <- sd_prior *
                (if ("sdY" %in% nd) d$sdY else sdY.est(d$varbeta, d$MAF, d$N))
    } else {
        z <- qnorm(0.5 * d$pvalues, lower.tail = FALSE)
        V <- if (d$type == "quant") Var.data(d$MAF, d$N)
             else Var.data.cc(d$MAF, d$N, d$s)
    }
    r <- sd_prior^2 / (sd_prior^2 + V)
    data.frame(z = z, V = V, r = r, lABF = 0.5 * (log(1 - r) + (r * z^2)))
}
