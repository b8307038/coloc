## ---------------------------------------------------------------------------
## The prior and control objects.
##
## coloc_prior() holds p1, p2, p12 and any prior weights, and is the one place
## their defaults are written. A prior is checked when it is made, and again
## by coloc() and coloc.bf_bf(), as they can be given a plain list, or a
## prior edited with $<-, which keeps its class.
## ---------------------------------------------------------------------------

new_coloc_prior <- function(p1, p2, p12, weights1, weights2) {
    x <- structure(list(p1 = p1, p2 = p2, p12 = p12,
                        weights1 = weights1, weights2 = weights2),
                   class = "coloc_prior")
    validate_coloc_prior(x)
}

validate_coloc_prior <- function(x) {
    ## first, as the bounds below, with any() and pmin(), would pass a vector
    for (nm in c("p1", "p2", "p12"))
        if (length(x[[nm]]) != 1)
            stop("p1, p2 and p12 must each be a scalar; ",
                 "to vary the prior across SNPs, use weights1/weights2")
    ## p12 may equal either bound, but p1 * p2 can round above it: 0.2 * 0.2 is
    ## 0.04000000000000001, so eps lets p12 = 0.04 pass
    eps <- sqrt(.Machine$double.eps)
    if (any(x$p12 < x$p1 * x$p2 * (1 - eps)) || any(x$p12 > pmin(x$p1, x$p2) * (1 + eps)))
        stop("p12 must lie in [p1*p2, min(p1, p2)]")
    for (w in list(x$weights1, x$weights2))
        if (!is.null(w) && any(is.na(w) | w <= 0))
            stop("prior weights must be non-NA and positive")
    x
}

##' Prior probabilities for a colocalisation analysis
##'
##' One object, validated once, in place of the same three probabilities being
##' re-declared in every entry point.
##'
##' `p1`, `p2` and `p12` are each the prior probability for any one SNP: that
##' it is causal for trait 1 only, trait 2 only, or both. Each is a single
##' number. To let the prior vary across SNPs, give weights, which scale it
##' SNP by SNP.
##'
##' @param p1 the prior probability for any one SNP that it is associated with
##'   trait 1 only
##' @param p2 the prior probability for any one SNP that it is associated with
##'   trait 2 only
##' @param p12 the prior probability for any one SNP that it is associated with
##'   both traits
##' @param weights1 optional: one weight per SNP of trait 1, positive, no `NA`
##' @param weights2 optional: one weight per SNP of trait 2, positive, no `NA`
##' @return an object of class coloc_prior
##' @export
##' @examples
##' coloc_prior()                    # the defaults
##' coloc_prior(p12 = 1e-5)          # coloc.abf()'s default p12
##' ## favour SNP 3 of 5 for trait 1, ten to one
##' coloc_prior(weights1 = c(1, 1, 10, 1, 1))
coloc_prior <- function(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                        weights1 = NULL, weights2 = NULL) {
    prior <- new_coloc_prior(p1, p2, p12, weights1, weights2)
    ## coloc.abf()'s advice on weights for both traits, so coloc() gives it too
    if (!is.null(weights1) && !is.null(weights2))
        warning("Prior weights specified for both traits.\n",
                "  The two weight vectors should be derived from independent sources of\n",
                "  information, i.e. not computed using the same method.")
    prior
}

##' Control parameters for a colocalisation analysis
##'
##' @param overlap.min minimum proportion of a signal's posterior that must
##'   fall on snps shared with the other trait before the pair is considered
##' @param trim_by_posterior whether to drop signal pairs failing overlap.min
##' @return an object of class coloc_control
##' @export
##' @examples
##' coloc_control()
##' ## compare pairs of signals even when only 10% of a signal's posterior
##' ## falls on SNPs the other trait has
##' coloc_control(overlap.min = 0.1)
coloc_control <- function(overlap.min = 0.5, trim_by_posterior = TRUE) {
    ## is.na() before the comparisons, so an NA gets this message rather than
    ## R's "missing value where TRUE/FALSE needed"
    if (!is.numeric(overlap.min) || length(overlap.min) != 1 || is.na(overlap.min) ||
        overlap.min <= 0 || overlap.min > 1)
        stop("overlap.min must be a single number above 0 and at most 1")
    if (!is.logical(trim_by_posterior) || length(trim_by_posterior) != 1 ||
        is.na(trim_by_posterior))
        stop("trim_by_posterior must be TRUE or FALSE")
    structure(list(overlap.min = overlap.min, trim_by_posterior = trim_by_posterior),
              class = "coloc_control")
}
