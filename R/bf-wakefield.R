## Wakefield's approximate Bayes factors, one per snp

##' Wakefield's approximate Bayes factor for each snp
##'
##' `wakefield_bf()` gives Wakefield's approximate Bayes factor (ABF) for each
##' snp in a dataset, on the log scale: the evidence for an effect at that snp
##' against no effect.
##'
##' Wakefield showed we could get ABF as
##'
##' log ABF = (log(1 - r) + r z^2) / 2
##'
##' We assume the true effect size ~ N(0,W), and set r = W / (W + V)
##' where V is the variance of the estimated effect size. z is the
##' usual Z score (beta/se(beta)). W is given by the `effect_priors`
##' parameter.
##'
##' Where the dataset includes beta and varbeta, this is
##' straightforward. Where it does not, we approximate the result from
##' pvalues, MAF and N. Please give beta and varbeta where available.
##'
##' For a quantitative trait the effect prior is in units of the
##' trait's standard deviation, so it is multiplied by `sdY`; when the
##' dataset has no `sdY`, it is estimated from `varbeta`, `MAF` and
##' `N`, with a warning. Again, please give sdY if known.
##'
##' @param d a dataset made with [coloc_dataset()], or a list, which it
##'     converts with [as_coloc_dataset()] and says so once per session
##' @param effect_priors the prior standard deviation of a causal snp's
##'     effect, named by the type of trait: `quant` in units of the trait's
##'     standard deviation, and `cc` on the log odds ratio scale. This option
##'     exists for special cases - do not change the standard priors unless
##'     you know what you are doing!
##' @return a data frame with one row per snp and the columns `snp`,
##'     `position` (`NA` when the dataset has none), `z`, `V` (the variance
##'     of the effect estimate), `r` (the shrinkage factor) and `lABF` (the
##'     log ABF)
##' @references Wakefield J (2009). Bayes factors for genome-wide association
##'     studies: comparison with P-values. Genetic Epidemiology 33(1),
##'     79-86. \doi{10.1002/gepi.20359}
##' @seealso [finemap.abf()], which turns these into posterior probabilities
##' [check_dataset()], which gives more detail on how to structure a coloc dataset
##' @author Wei-Yu Lin, Claudia Giambartolomei, Vincent Plagnol, Chris Wallace
##' @export
##' @examples
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##' head(wakefield_bf(D1))
wakefield_bf <- function(d, effect_priors = c(quant = 0.15, cc = 0.2)) {
    if (!is.numeric(effect_priors) ||
        !all(c("quant", "cc") %in% names(effect_priors)) ||
        !all(is.finite(effect_priors[c("quant", "cc")]) &
             effect_priors[c("quant", "cc")] > 0))
        stop("effect_priors should be positive numbers named quant and cc")
    if (inherits(d, "coloc_dataset")) {
        ## checked again, as coloc's functions do: its elements may have been
        ## changed since it was made
        check_dataset(d)
    } else {
        d <- as_coloc_dataset(d)
        .note_once("plain_list", .plain_list_note("wakefield_bf"))
    }
    nd <- names(d)
    ## by the type's label: [[ would take a factor type by its code
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
    data.frame(snp = as.character(d$snp),
               position = if ("position" %in% nd) d$position
                          else rep(NA_integer_, length(z)),
               z = z, V = V, r = r, lABF = 0.5 * (log(1 - r) + (r * z^2)),
               row.names = NULL)
}
