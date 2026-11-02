##' Check coloc dataset inputs for errors
##'
##' A coloc dataset is a list, containing a mixture of vectors
##' capturing quantities that vary between snps (these vectors must
##' all have equal length) and scalars capturing quantities that
##' describe the dataset.
##'
##' Coloc is flexible, requiring perhaps only p values, or z scores, or effect
##' estimates and standard errors, but with this flexibility, also comes
##' difficulties describing exactly the combinations of items required.
##'
##' Required vectors are some subset of
##'
##' \describe{
##'   \item{beta}{regression coefficient for each SNP from dataset 1}
##'   \item{varbeta}{variance of beta}
##'   \item{pvalues}{P-values for each SNP in dataset 1}
##'   \item{MAF}{minor allele frequency of the variants}
##'   \item{snp}{a character vector of snp ids, optional. It will be used to merge dataset1 and dataset2 and will be retained in the results.}
##' }
##'
##' Preferably, give \code{beta} and \code{varbeta}.  But if these are not available, sufficient statistics can be approximated from \code{pvalues} and \code{MAF}.
##'
##' Required scalars are some subset of
##'
##' \describe{
##'   \item{N}{Number of samples in dataset 1}
##'   \item{type}{the type of data in dataset 1 - either "quant" or "cc" to denote quantitative or case-control}
##'   \item{s}{for a case control dataset, the proportion of samples in dataset 1 that are cases}
##'   \item{sdY}{for a quantitative trait, the population standard deviation of the trait.  if not given, it can be estimated from the vectors of varbeta and MAF}
##' }
##'
##' You must always give {\code{type}}.  Then,
##' \describe{
##' \item{if \code{type}=="cc"}{\code{s}}
##' \item{if \code{type}=="quant" and \code{sdY} known}{\code{sdY}}
##' \item{if beta, varbeta not known}{\code{N}}
##' }
##' If \code{sdY} is unknown, it will be approximated, and this will require
##' \describe{
##' \item{summary data to estimate \code{sdY}}{\code{beta}, \code{varbeta}, \code{N}, \code{MAF}}
##' }
##'
##' Beyond what must be present, \code{check_dataset} refuses values coloc
##' can't use. Wherever they are given: a \code{type} or \code{s} that isn't a
##' single value, an \code{N} that isn't a finite number above 0, and an
##' \code{LD} without two dimensions. Where they are read: a \code{beta} or
##' \code{varbeta} that isn't numeric, or a negative \code{varbeta}, when both
##' are given; \code{pvalues} that aren't numeric, when \code{beta} and
##' \code{varbeta} aren't both given; and an \code{sdY} that isn't a single
##' finite number above 0, on a quantitative trait given as \code{beta} and
##' \code{varbeta}, unless \code{req} includes "LD", as SuSiE doesn't read
##' \code{sdY}. When \code{req} includes "LD", as for SuSiE, \code{LD} must be
##' a correlation matrix for the dataset's snps: numeric, with values from -1
##' to 1, symmetric, and 1 on the diagonal. A \code{z}, if given, must then
##' have a finite number for each snp.
##'
##' Optional vectors are
##'
##' \describe{
##'   \item{position}{a vector of snp positions, required for \code{plot_dataset}}
##'   \item{z}{z-scores, which \code{runsusie} fine-maps from in place of \code{beta / sqrt(varbeta)}}
##' }
##'
##' \code{check_dataset} calls stop() unless a series of expectations on dataset
##' input format are met
##'
##' This is a helper function for use by other coloc functions, but
##' you can use it directly to check the format of a dataset to be
##' supplied to coloc.abf(), coloc.signals(), finemap.abf(), or
##' finemap.signals().
##' @title check_dataset
##' @param d dataset to check
##' @param suffix a label for the dataset, such as 1 or 2: messages then begin
##'     "dataset 1: ". With "" they have no prefix
##' @param req names of elements that must be present. When it includes "LD",
##'     as \code{runsusie}'s does, \code{LD} must also be a correlation matrix,
##'     and a \code{z} must have a finite number for each snp
##' @param warn.minp print warning if no p value < warn.minp
##' @return NULL if no errors found
##' @export
##' @author Chris Wallace
##' @examples
##' data(coloc_test_data)
##' check_dataset(coloc_test_data$D1)     # NULL: no problems found
##' bad <- coloc_test_data$D1
##' bad$beta[1] <- NA
##' try(check_dataset(bad, 1))
check_dataset <- function(d,
                          suffix="",
                          req=c("type","snp"),
                          warn.minp=1e-6) {
    ## A label, such as 1 or 2, starts each message as "dataset 1: "; "", NULL,
    ## NA or more than one value is no label. No call is shown, as the message
    ## names the dataset.
    dataset <- if (length(suffix) == 1L && !is.na(suffix) && !identical(as.character(suffix), ""))
                   suffix
    prefix <- if (is.null(dataset)) "" else paste0("dataset ", dataset, ": ")
    fail <- function(...) stop(prefix, ..., call. = FALSE)

    if (!is.list(d))
        fail(if (is.null(dataset)) "the dataset is not a list" else "is not a list")
    nd <- intersect(names(d), .dataset_items)          # the elements coloc reads
    p_route <- !("beta" %in% nd && "varbeta" %in% nd)  # p-values unless both are given
    susie <- "LD" %in% req                             # SuSiE reads LD, N and z, not sdY

    ## The checks stop at the first problem, so their order decides which
    ## message a dataset with two problems gets.
    if (length(setdiff(req, nd)))
        fail("missing required element(s) ", paste(setdiff(req, nd), collapse = ", "))
    ## no missing values - make people clean their own data rather than make
    ## assumptions here for datasets I don't know. z is checked where SuSiE
    ## reads it, below
    for (v in setdiff(nd, "z"))
        if (any(is.na(d[[v]])))
            fail(v, " contains missing values")
    if ("type" %in% nd) {
        if (length(d$type) != 1L)
            fail("type must be a single value, quant or cc")
        if (!(d$type %in% c("quant", "cc")))
            fail("type must be quant or cc")
    }
    if ("snp" %in% nd && any(duplicated(d$snp)))
        fail("duplicated snps found")
    if ("snp" %in% nd && is.factor(d$snp))
        fail("snp should be a character vector but is a factor")
    if ("MAF" %in% nd && (!is.numeric(d$MAF) || any(d$MAF <= 0 | d$MAF >= 1)))
        fail("MAF should be a numeric, strictly >0 & <1")
    ## each per-snp element but z has as many values as snp, or without snp as
    ## the first of them, in .per_snp_items' order; z is compared with the
    ## same element below
    per_snp <- intersect(setdiff(.per_snp_items, "z"), nd)
    ref <- if ("snp" %in% per_snp) "snp" else per_snp[1]
    for (v in per_snp)
        if (length(d[[v]]) != length(d[[ref]]))
            fail(ref, " has ", length(d[[ref]]), if (length(d[[ref]]) == 1L) " value" else " values",
                 " but ", v, " has ", length(d[[v]]))
    if ("N" %in% nd && (!is.numeric(d$N) || length(d$N) != 1L))
        fail("N must be a single number; coloc does not take per-snp sample sizes")
    if (!("type" %in% nd))
        fail("variable type not set")
    if ("s" %in% nd) {
        if (length(d$s) != 1L)
            fail("s must be a single number between 0 and 1")
        if (!is.numeric(d$s) || d$s <= 0 || d$s >= 1)
            fail("s must be between 0 and 1")
    }

    ## without beta and varbeta, coloc works from p-values, MAF and N
    if (p_route) {
        if (!("pvalues" %in% nd && "MAF" %in% nd))
            fail("require p values and MAF if beta, varbeta are unavailable")
        ## p-values that aren't numbers would be compared as strings here;
        ## they are refused below
        if (is.numeric(d$pvalues) && any(d$pvalues <= 0 | d$pvalues > 1))
            fail("pvalues must be above 0 and at most 1")
        if (d$type == "cc" && !("s" %in% nd))
            fail("require, s, proportion of samples who are cases, if beta, varbeta are unavailable")
        if (!("N" %in% nd) || is.null(d$N) || any(d$N <= 0))
            fail("sample size N <=0 or not set")
        if (!is.numeric(d$pvalues))
            fail("pvalues must be numeric")
    } else {
        for (v in c("beta", "varbeta"))
            if (!is.numeric(d[[v]]))
                fail(v, " must be numeric")
    }
    if (any(is.infinite(d$beta) | is.infinite(d$varbeta)))
        fail("Infinite values in beta and/or varbeta")
    if (any(d$varbeta == 0))
        fail("Zero values in varbeta")
    if (!p_route && any(d$varbeta < 0))
        fail("varbeta must be positive")

    ## a warning when no p-value is below warn.minp; a coloc_dataset had it
    ## when it was made
    if (!inherits(d, "coloc_dataset")) {
        p <- min(if (p_route) d$pvalues else pnorm(-abs(d$beta / sqrt(d$varbeta))) * 2)
        if (p > warn.minp)
            warning(warningCondition(
                paste0(prefix, "minimum p value is: ", format.pval(p),
                       "\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument. If that's not the explanation, please check the 02_data vignette."),
                class = "coloc_min_p", call = NULL))
    }

    if (d$type == "quant" && !("sdY" %in% nd) && !("MAF" %in% nd && "N" %in% nd))
        fail("must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
    if ("LD" %in% nd) {
        if (length(dim(d$LD)) < 2L)
            fail("LD must be a matrix")
        if (nrow(d$LD) != ncol(d$LD))
            fail("LD not square")
        if (!identical(colnames(d$LD), rownames(d$LD)))
            fail("LD rownames != colnames")
        if (length(setdiff(d$snp, colnames(d$LD))))
            fail("colnames in LD do not contain all SNPs")
        if (length(dim(d$LD)) > 2L)
            fail("LD must be a matrix")
    }
    if ("N" %in% nd && (!is.finite(d$N) || d$N <= 0))
        fail("N must be a finite number above 0")
    ## sdY is read for a quantitative trait on the beta route, and not by SuSiE
    if ("sdY" %in% nd && !p_route && d$type == "quant" && !susie &&
        (!is.numeric(d$sdY) || length(d$sdY) != 1L || !is.finite(d$sdY) || d$sdY <= 0))
        fail("sdY must be a single finite number above 0")
    ## Where sdY is read and MAF and N are given, an sdY more than twice or
    ## under half its estimate from varbeta, MAF and N gets a warning, as when
    ## varbeta holds standard errors. It doesn't stop, as covariates or a mixed
    ## model can move a correct sdY away from the estimate. A coloc_dataset
    ## gets it too, as it may have been edited since it was made.
    if (all(c("sdY", "MAF", "N") %in% nd) && !p_route && d$type == "quant" && !susie &&
        length(d$varbeta)) {
        fit <- .sdY_fit(d$varbeta, d$MAF, d$N)
        if (is.finite(fit) && fit > 0 && max(d$sdY / sqrt(fit), sqrt(fit) / d$sdY) > 2)
            warning(warningCondition(
                paste0(prefix, sprintf("sdY is %.3g, but varbeta, MAF and N suggest %.3g: does varbeta hold variances?",
                                       d$sdY, sqrt(fit))),
                class = "coloc_sdY", call = NULL))
    }

    ## what SuSiE reads; LD is given, as req names it
    if (susie) {
        if (!(is.matrix(d$LD) && is.numeric(d$LD)))
            fail("LD must be a numeric matrix")
        if ("N" %in% nd && d$N <= 1)
            fail("N must be above 1 for SuSiE")
        if ("z" %in% nd) {
            if (any(is.na(d$z)))
                fail("z contains missing values")
            if (length(d$z) != length(d[[ref]]))
                fail(ref, " has ", length(d[[ref]]), if (length(d[[ref]]) == 1L) " value" else " values",
                     " but z has ", length(d$z))
            if (!is.numeric(d$z))
                fail("z must be numeric")
            if (any(is.infinite(d$z)))
                fail("Infinite values in z")
        }
        ## only the dataset's snps, which are all SuSiE uses, copied out only
        ## when LD covers more; with as many, they are the same snps
        LD <- if (ncol(d$LD) > length(d$snp)) d$LD[d$snp, d$snp] else d$LD
        if (!.is_correlation(LD))
            fail("LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal")
    }
    NULL
}

#'@rdname check_dataset
#' @param ... arguments passed to check_dataset()
#'@export
check.dataset=function(...) {
  warning("Deprecated, use check_dataset() in future")
  check_dataset(...)
}


## LD as a correlation matrix: values in [-1, 1], 1 on the diagonal, and equal
## to its transpose, each to within tol, with no missing values. Any LD worked
## out in 64-bit arithmetic passes with room to spare, and so does one stored
## as 32-bit floats. Only one worked out in 32-bit arithmetic by rescaling each
## snp to sd 1 first can miss by more, about 2e-5 at 1,000 people, and is
## refused. A covariance matrix or a value of 1.5 is far outside it. The
## transpose is compared in blocks of 500 columns, so a large LD isn't copied
## whole.
.is_correlation <- function(LD, tol = 1e-6) {
    if (anyNA(LD))
        return(FALSE)
    r <- range(LD)
    if (r[1] < -1 - tol || r[2] > 1 + tol)
        return(FALSE)
    if (any(abs(diag(LD) - 1) > tol))
        return(FALSE)
    for (j in split(seq_len(ncol(LD)), ceiling(seq_len(ncol(LD)) / 500)))
        if (max(abs(LD[, j, drop = FALSE] - t(LD[j, , drop = FALSE]))) > tol)
            return(FALSE)
    TRUE
}

check_ld <- function(D,LD) {
    if(is.null(LD))
      stop("LD required")
    if(!is.matrix(LD))
      stop("LD must be of class matrix")
    ## if(any(LD[upper.tri(LD)]==1))
    ##   stop("LD includes SNPs in perfect LD")
    if(nrow(LD)!=ncol(LD))
        stop("LD not square")
    if(!identical(colnames(LD),rownames(LD)))
      stop("LD rownames != colnames")
    if(length(setdiff(D$snp,colnames(LD))))
        stop("colnames in LD do not contain all SNPs")
    ## finemap.signals() and coloc.signals() compute with LD, so it must be a
    ## numeric correlation matrix, as for SuSiE. Only D's snps are checked, or
    ## all of LD when D names none
    if(!is.numeric(LD))
        stop("LD must be a numeric matrix")
    if(!.is_correlation(if(length(D$snp) && ncol(LD) > length(D$snp)) LD[D$snp,D$snp] else LD))
        stop("LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal")
}

##' check alignment between beta and LD
##'
##' @title check alignment
##' @param D a coloc dataset
##' @param thr count and plot the pairs of SNPs in absolute LD > thr, each
##'     pair once
##' @param do_plot if TRUE (default) plot the diagnostic
##' @export
##' @return proportion of those pairs whose product of z scores, divided by
##'     their LD, is positive
##' @author Chris Wallace
##' @examples
##' data(coloc_test_data)
##' ## close to 1 when beta and LD are aligned to the same alleles
##' check_alignment(coloc_test_data$D3, do_plot = FALSE)
check_alignment <- function(D,thr=0.2,do_plot=TRUE) {
  check_dataset(D)
  ## z_i * z_j / LD_ij for each pair of snps i < j with |LD_ij| above thr,
  ## positive when the two z-scores agree in sign with their LD. LD is matched
  ## to the snps by name, and read in blocks of 500 columns to bound the memory
  z <- D$beta/sqrt(D$varbeta)
  at <- match(D$snp, colnames(D$LD))
  tmp <- list()
  for (j in split(seq_along(at), ceiling(seq_along(at) / 500))) {
    i <- seq_len(max(j))
    R <- D$LD[at[i], at[j], drop = FALSE]
    pair <- outer(i, j, "<") & abs(R) > thr
    tmp[[length(tmp) + 1L]] <- (outer(z[i], z[j]) / R)[pair]
  }
  tmp <- unlist(tmp)
  if(do_plot) {
  hist(tmp,
       xlab="ratio of product of Z scores to LD",
       main="alignment check plot",
       sub="expect most values to be positive\nsymmetry is a warning sign of potentially poor alignment")
  legend("topright",legend=paste0("% positive = ",100*round(mean( tmp > 0 ),3)))
  abline(v=0,col="red")
  }
  return(mean(tmp > 0))
}
#'@rdname check_alignment
#'@export
#' @param ... arguments passed to check_alignment()
check.alignment=function(...) {
  warning("Deprecated, use check_alignment() in future")
  check_alignment(...)
}
