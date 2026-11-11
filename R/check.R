##' Check a coloc dataset
##'
##' `check_dataset()` checks that a dataset has the elements coloc needs, and
##' that coloc can use their values. It stops at the first problem it finds,
##' with a message that names the element and says what is wrong. You seldom
##' need to call it, as [coloc_dataset()] and [as_coloc_dataset()] check a
##' dataset when they make it, and [coloc()], [wakefield_bf()], [runsusie()]
##' and the older functions, such as [coloc.abf()], check every dataset they
##' are given. Call it yourself to check your datasets before you analyse
##' them, for example every region of a study before a loop over the regions,
##' or to check a dataset you have changed.
##'
##' A coloc dataset is a list that describes one trait in one region. Some of
##' its elements have one value per snp, and these must all have the same
##' length:
##'
##' * `snp`, the snp ids, which coloc matches the two traits on;
##' * `beta`, the regression coefficient of each snp, and `varbeta`, its
##'   variance, the standard error squared;
##' * `pvalues`, the p-value of each snp, and `MAF`, its minor allele
##'   frequency;
##' * `position`, the position of each snp, which [plot_dataset()] needs;
##' * `z`, the z-score of each snp, which [runsusie()] fine-maps from in place
##'   of `beta / sqrt(varbeta)`.
##'
##' The other elements are single values that describe the study, and the LD
##' matrix:
##'
##' * `type`, `"quant"` for a quantitative trait or `"cc"` for a case-control
##'   trait;
##' * `N`, the number of samples;
##' * `s`, for a case-control trait, the proportion of samples that are cases;
##' * `sdY`, for a quantitative trait, the standard deviation of the trait;
##' * `LD`, the correlation matrix of the snps, with snp ids as its row and
##'   column names, which fine-mapping with SuSiE needs.
##'
##' Which elements are needed depends on the data. Always give `type` and
##' `snp`. Give `beta` and `varbeta` if you have them. Without them, coloc
##' approximates them from `pvalues`, `MAF` and `N`, and a case-control trait
##' also needs `s`. A quantitative trait given as `beta` and `varbeta` needs
##' `sdY`, or `MAF` and `N` so coloc can estimate `sdY`.
##'
##' `check_dataset()` checks in three steps, and stops at the first problem:
##'
##' 1. Each element, wherever it is given. The elements that `req` names must
##'    be present, and no element but `z` may have missing values. `type` must
##'    be given, as one value, `"quant"` or `"cc"`. `s` must be one number
##'    between 0 and 1, `MAF` numbers above 0 and below 1, and `N` one finite
##'    number above 1. Each snp id must appear once, and `snp` must not be a
##'    factor. The per-snp elements must have the same length. `beta` and
##'    `varbeta`, when both are given, must have no infinite values, and
##'    `varbeta` must have no zeros. `LD` must have two dimensions and be
##'    square, with row names that match its column names, in the same order,
##'    and cover every snp.
##' 2. What coloc reads on the dataset's route. Without both `beta` and
##'    `varbeta`, coloc reads `pvalues`, which must be numbers above 0 and at
##'    most 1, `MAF` and `N`, and `s` for a case-control trait. With them,
##'    `beta` and `varbeta` must be numeric and `varbeta` not negative, and a
##'    quantitative trait needs `sdY`, one finite number above 0, or `MAF` and
##'    `N` so coloc can estimate it.
##' 3. What SuSiE reads, when `req` includes `"LD"`, as [runsusie()]'s does.
##'    `LD` must be a numeric correlation matrix for the dataset's snps, with
##'    values from -1 to 1, symmetric, and 1 on the diagonal. A `z`, if given,
##'    must have a finite number for each snp. SuSiE doesn't read `sdY`, so
##'    step 2 doesn't check it on this route.
##'
##' It also warns, without stopping, in two cases. When no p-value is below
##' `warn.minp`, the trait may have no association in the region; a
##' [coloc_dataset()] gets this warning once, when it is made. When a
##' quantitative trait given as `beta` and `varbeta` has an `sdY` more than
##' twice or under half its estimate from `varbeta`, `MAF` and `N`, `varbeta`
##' may hold standard errors in place of variances; this warning comes each
##' time the dataset is checked, except on the SuSiE route.
##'
##' @param d a dataset, as a list or a [coloc_dataset()], but not a data.frame
##' @param suffix a label for the dataset, such as 1 or 2, which begins each
##'     message as "dataset 1: ". With "" messages have no label
##' @param req the elements that must be present. When it includes `"LD"`, as
##'     [runsusie()]'s does, `check_dataset()` also checks what SuSiE reads,
##'     step 3 in Details
##' @param warn.minp warn when no p-value is below this
##' @return `NULL` when it finds no problem. Otherwise it stops with an error
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

    if (!is.list(d) || is.data.frame(d)) # a data.frame repeats type, sdY and N on every row
        .stop_not_a_dataset(d, dataset)
    nd <- intersect(names(d), .dataset_items)          # the elements coloc reads
    p_route <- !("beta" %in% nd && "varbeta" %in% nd)  # p-values unless both are given
    susie <- "LD" %in% req                             # SuSiE reads LD, N and z, not sdY
    ## for messages: "pvalues is missing", "pvalues and MAF are missing"
    is_missing <- function(v)
        paste(if (length(v) < 2L) v else paste(paste(v[-length(v)], collapse = ", "), "and", v[length(v)]),
              ngettext(length(v), "is missing", "are missing"))

    ## The checks stop at the first problem, in three blocks: what each element
    ## must be wherever it is given, then what each route reads, then what
    ## SuSiE reads. So an element's shape is checked before a route reads it.

    ## 1. Each element, wherever it is given
    absent <- setdiff(req, nd)
    if (length(absent))
        fail("missing required element(s) ", paste(absent, collapse = ", "))
    ## no missing values - make people clean their own data rather than make
    ## assumptions here for datasets I don't know. z is checked where SuSiE
    ## reads it, below
    for (v in setdiff(nd, "z"))
        if (any(is.na(d[[v]])))
            fail(v, " contains missing values")
    if (!("type" %in% nd))
        fail('type is missing; set type to "quant" for a quantitative trait or "cc" for a case-control trait')
    if (length(d$type) != 1L)
        fail("type must be a single value, quant or cc")
    if (!(d$type %in% c("quant", "cc")))
        fail("type must be quant or cc")
    if ("s" %in% nd) {
        if (length(d$s) != 1L)
            fail("s must be a single number between 0 and 1")
        if (!is.numeric(d$s) || d$s <= 0 || d$s >= 1)
            fail("s must be between 0 and 1")
    }
    if ("snp" %in% nd && any(duplicated(d$snp)))
        fail("duplicated snps found")
    if ("snp" %in% nd && is.factor(d$snp))
        fail("snp should be a character vector but is a factor")
    if ("MAF" %in% nd && (!is.numeric(d$MAF) || any(d$MAF <= 0 | d$MAF >= 1)))
        fail("MAF should be a numeric, strictly >0 & <1")
    ## each per-snp element but z has as many values as snp, or, without snp,
    ## as the first of them in .per_snp_items' order; z is compared with the
    ## same element below
    per_snp <- intersect(setdiff(.per_snp_items, "z"), nd)
    ref <- if ("snp" %in% per_snp) "snp" else per_snp[1]
    for (v in per_snp)
        if (length(d[[v]]) != length(d[[ref]]))
            fail(ref, " has ", length(d[[ref]]), if (length(d[[ref]]) == 1L) " value" else " values",
                 " but ", v, " has ", length(d[[v]]))
    if ("N" %in% nd) {
        if (!is.numeric(d$N) || length(d$N) != 1L)
            fail("N must be a single number; coloc does not take per-snp sample sizes")
        if (!is.finite(d$N) || d$N <= 1)
            fail("N must be finite and > 1")
    }
    if (any(is.infinite(d$beta) | is.infinite(d$varbeta))) {
        v <- if (any(is.infinite(d$beta))) "beta" else "varbeta"
        fail(v, " has infinite values")
    }
    ## a varbeta that isn't a number is refused where it is read, below
    if (is.numeric(d$varbeta) && any(d$varbeta == 0))
        fail("varbeta must be > 0")
    if ("LD" %in% nd) {
        if (length(dim(d$LD)) < 2L)
            fail("LD must be a matrix")
        if (nrow(d$LD) != ncol(d$LD))
            fail("LD must be square, but it has ", nrow(d$LD), ngettext(nrow(d$LD), " row", " rows"),
                 " and ", ncol(d$LD), ngettext(ncol(d$LD), " column", " columns"))
        if (!identical(colnames(d$LD), rownames(d$LD)))
            fail("LD rownames != colnames")
        if (length(setdiff(d$snp, colnames(d$LD))))
            fail("LD must cover every snp")
        if (length(dim(d$LD)) > 2L)
            fail("LD must be a matrix")
    }

    ## 2. What each route reads
    if (p_route) {
        ## without beta and varbeta, coloc works from p-values, MAF and N
        if (!("pvalues" %in% nd && "MAF" %in% nd))
            fail(is_missing(setdiff(c("pvalues", "MAF"), nd)),
                 "; without beta and varbeta, coloc needs pvalues, MAF and N")
        if (!is.numeric(d$pvalues))
            fail("pvalues must be numeric")
        if (any(d$pvalues <= 0 | d$pvalues > 1))
            fail("pvalues must be above 0 and at most 1")
        if (d$type == "cc" && !("s" %in% nd))
            fail("s is missing; a case-control dataset without beta and varbeta needs s, ",
                 "the proportion of samples who are cases")
        if (!("N" %in% nd))
            fail("N is missing; without beta and varbeta, coloc needs N, the sample size")
    } else {
        for (v in c("beta", "varbeta"))
            if (!is.numeric(d[[v]]))
                fail(v, " must be numeric")
        if (any(d$varbeta < 0))
            fail("varbeta must be > 0")
        if (d$type == "quant" && !("sdY" %in% nd) && !("MAF" %in% nd && "N" %in% nd))
            fail("must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
        ## sdY is read for a quantitative trait, but not by SuSiE
        if ("sdY" %in% nd && d$type == "quant" && !susie &&
            (!is.numeric(d$sdY) || length(d$sdY) != 1L || !is.finite(d$sdY) || d$sdY <= 0))
            fail("sdY must be a single finite number above 0")
    }
    ## a warning when no p-value is below warn.minp; a coloc_dataset had it
    ## when it was made
    if (!inherits(d, "coloc_dataset")) {
        p <- min(if (p_route) d$pvalues else pnorm(-abs(d$beta / sqrt(d$varbeta))) * 2)
        if (p > warn.minp)
            warning(warningCondition(
                paste0(prefix, "the smallest p-value is ", format.pval(p),
                       ", so this dataset may have no association in the region.",
                       "\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument. If that's not the explanation, please check the 02_data vignette."),
                class = "coloc_min_p", call = NULL))
    }
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
                paste0(prefix, sprintf("sdY is %.3g, but varbeta, MAF and N suggest %.3g; ", d$sdY, sqrt(fit)),
                       "check that varbeta holds variances, not standard errors, ",
                       "and that sdY is on the scale of beta"),
                class = "coloc_sdY", call = NULL))
    }

    ## 3. What SuSiE reads; LD is given, as req names it
    if (susie) {
        if (!(is.matrix(d$LD) && is.numeric(d$LD)))
            fail("LD must be a numeric matrix")
        ## only the dataset's snps, which are all SuSiE uses, copied out only
        ## when LD covers more; with as many, they are the same snps
        LD <- if (ncol(d$LD) > length(d$snp)) d$LD[d$snp, d$snp] else d$LD
        if (!.is_correlation(LD))
            fail("LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal")
        if ("z" %in% nd) {
            if (any(is.na(d$z)))
                fail("z contains missing values")
            if (length(d$z) != length(d[[ref]]))
                fail(ref, " has ", length(d[[ref]]), if (length(d[[ref]]) == 1L) " value" else " values",
                     " but z has ", length(d$z))
            if (!is.numeric(d$z))
                fail("z must be numeric")
            if (any(is.infinite(d$z)))
                fail("z has infinite values")
        }
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


## The refusal of anything but a list of dataset elements, for check_dataset()
## and new_coloc_dataset(). With a label, such as 1, the message begins
## "dataset 1: ".
.stop_not_a_dataset <- function(d, dataset = NULL)
    stop(if (is.null(dataset)) "the dataset" else paste0("dataset ", dataset, ": this"),
         " is ", .what_is(d), ", not a list of coloc_dataset elements (snp, beta, varbeta, ...)",
         call. = FALSE)

## TRUE when LD is a correlation matrix, to within tol: no missing values,
## values in [-1, 1], 1 on the diagonal, and symmetric. The default passes an
## LD worked out in 64-bit arithmetic or stored as 32-bit floats, but one
## worked out in 32-bit arithmetic can be out by about 2e-5, and is refused.
## The transpose is compared in blocks of 500 columns, so a large LD isn't
## copied whole.
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
      stop("LD is missing")
    if(!is.matrix(LD))
      stop("LD must be a matrix")
    ## if(any(LD[upper.tri(LD)]==1))
    ##   stop("LD includes SNPs in perfect LD")
    if(nrow(LD)!=ncol(LD))
        stop("LD must be square, but it has ", nrow(LD), ngettext(nrow(LD), " row", " rows"),
             " and ", ncol(LD), ngettext(ncol(LD), " column", " columns"))
    if(!identical(colnames(LD),rownames(LD)))
      stop("LD rownames != colnames")
    if(length(setdiff(D$snp,colnames(LD))))
        stop("LD must cover every snp")
    ## finemap.signals() and coloc.signals() compute with LD, so it must be a
    ## numeric correlation matrix, as for SuSiE. Only D's snps are checked, or
    ## all of LD when D names none
    if(!is.numeric(LD))
        stop("LD must be a numeric matrix")
    if(length(D$snp) && ncol(LD) > length(D$snp))
        LD <- LD[D$snp,D$snp]
    if(anyNA(LD))
        stop("LD contains missing values")
    if(!.is_correlation(LD))
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
  check_dataset(D, req = c("type", "snp", "beta", "varbeta", "LD"))
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
