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
    ## the rules are .dataset_rules, below; the call is passed on so that an
    ## error still shows this call, not .gate()'s
    .gate(d, suffix, req, warn.minp, call = sys.call())
    ## if we reach here, no badness detected
    NULL
}

#'@rdname check_dataset
#' @param ... arguments passed to check_dataset()
#'@export
check.dataset=function(...) {
  warning("Deprecated, use check_dataset() in future")
  check_dataset(...)
}


## ---- check_dataset()'s rules ------------------------------------------------
## One entry of .dataset_rules per rule, in the order check_dataset() checks
## them: the order decides which message a dataset with two mistakes gets, and
## test-check.R pins it against 6.0.3's. 6.0.3's conditions are kept, some
## with a test in front for input they failed on with R's error, so a dataset
## 6.0.3 refused gets its message. A rule added since goes where 6.0.3 failed
## with R's error on what it refuses, or after 6.0.3's rules if 6.0.3 took it,
## so it doesn't take over one of 6.0.3's messages; test-check.R adds each new
## mistake to its cases to check that.

## One rule. check is function(d, ctx): NULL when the dataset passes, else the
## message without its prefix, which .gate() adds. A check may also return
## .passed(detail), a pass with what a report shows for it. about names the
## elements the rule reads, so check_ld() can pick out the LD rules and a
## report can say a rule's elements aren't given. level "warning" is the
## minimum p-value step alone: .gate() warns and goes on. level "report" is
## for the rules only the reports run
.rule <- function(id, check, about = character(0),
                  level = c("error", "warning", "report"))
    list(id = id, check = check, about = about, level = match.arg(level))

## the p-value route: no beta, or no varbeta
.p_route <- function(ctx) !("beta" %in% ctx$nd) || !("varbeta" %in% ctx$nd)

## the SuSiE route: a req naming LD, as .dataset_uses$susie does. SuSiE reads
## LD, N and z, and not sdY
.susie_route <- function(ctx) "LD" %in% ctx$req

## The message for per-snp elements v whose lengths differ, NULL if none: each
## against snp, or without snp against the first, in .per_snp_items' order, so
## the message names the same elements however the list was built
.lengths_differ <- function(d, v) {
    if (length(v) < 2L)
        return(NULL)
    n <- lengths(d[v])
    ref <- if ("snp" %in% v) "snp" else v[1]
    bad <- v[n != n[[ref]]]
    if (length(bad))
        paste(ref, "has", n[[ref]], if (n[[ref]] == 1L) "value" else "values",
              "but", bad[1], "has", n[[bad[1]]])
}

## LD as a correlation matrix: values in [-1, 1], 1 on the diagonal, and equal
## to its transpose, each to within tol. LDs worked out by different software,
## or in 32-bit arithmetic, differ from exact by about 1e-7; 1e-6 is ten times
## that, while a covariance matrix or a value of 1.5 is far outside it. The
## transpose is compared in blocks of 500 columns, so a large LD isn't copied
## whole.
.is_correlation <- function(LD, tol = 1e-6) {
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

## z_i * z_j / LD_ij for each pair of snps i < j with |LD_ij| above thr, and
## with z_thr above 0, both |z| above z_thr. The ratio is positive when the
## two z-scores agree in sign with their LD, as they should when beta and LD
## are for the same alleles. LD is matched to snp by name, and read in blocks
## of 500 columns, as .is_correlation() reads it
.z_ld_ratios <- function(z, LD, snp, thr = 0.2, z_thr = 0) {
    if (z_thr > 0) {
        keep <- abs(z) > z_thr
        z <- z[keep]
        snp <- snp[keep]
    }
    at <- match(snp, colnames(LD))
    out <- list()
    for (j in split(seq_along(at), ceiling(seq_along(at) / 500))) {
        i <- seq_len(max(j))
        R <- LD[at[i], at[j], drop = FALSE]
        pair <- outer(i, j, "<") & abs(R) > thr
        out[[length(out) + 1L]] <- (outer(z[i], z[j]) / R)[pair]
    }
    unlist(out)
}

## Where each element is read, which decides where a rule added since 6.0.3
## refuses a value coloc can't use: refusing it where nothing reads it would
## stop a dataset whose answer is right. Wakefield's beta route (beta and
## varbeta both given) reads beta, varbeta, and for a quantitative trait sdY,
## or MAF and N to estimate it. Its p-value route reads pvalues, MAF, N, and s
## for a case-control trait. SuSiE, the use whose req names LD, reads beta,
## varbeta, LD and N, and z when given, and not sdY. N is refused wherever it
## is given, as no sample size is 0 or below.

## ctx is list(nd, req, warn.minp), nd being the elements of d coloc knows
.dataset_rules <- list(
    .rule("required", function(d, ctx)
        if (length(setdiff(ctx$req, ctx$nd)))
            paste0("missing required element(s) ",
                   paste(setdiff(ctx$req, ctx$nd), collapse = ", "))),
    ## no missing values - make people clean their own data rather than make
    ## assumptions here for datasets I don't know. z, which 6.0.3 didn't know,
    ## is checked by z_valid below, after 6.0.3's rules
    .rule("no_na", function(d, ctx) {
        for (v in setdiff(ctx$nd, "z"))
            if (any(is.na(d[[v]])))
                return(paste0(v, " contains missing values"))
        NULL
    }),
    ## guarded, so a dataset with no type reaches "variable type not set"
    ## below, as runsusie() leaves type out of req. A type of any length but
    ## one failed here with R's error, and gets its own sentence
    .rule("type_valid", about = "type", function(d, ctx)
        if ("type" %in% ctx$nd) {
            if (length(d$type) != 1L)
                "type must be a single value, quant or cc"
            else if (!(d$type %in% c("quant", "cc")))
                "type must be quant or cc"
        }),
    .rule("snp_unique", about = "snp", function(d, ctx)
        if ("snp" %in% ctx$nd && any(duplicated(d$snp)))
            "duplicated snps found"),
    .rule("snp_character", about = "snp", function(d, ctx)
        if ("snp" %in% ctx$nd && is.factor(d$snp))
            "snp should be a character vector but is a factor"),
    .rule("maf_range", about = "MAF", function(d, ctx)
        if ("MAF" %in% ctx$nd && (!is.numeric(d$MAF) || any(is.na(d$MAF)) ||
                                  any(d$MAF <= 0) || any(d$MAF >= 1)))
            "MAF should be a numeric, strictly >0 & <1"),
    ## the per-snp elements (.per_snp_items, in R/dataset.R) but z, which
    ## z_valid checks
    .rule("lengths", function(d, ctx)
        .lengths_differ(d, intersect(setdiff(.per_snp_items, "z"), ctx$nd))),
    .rule("N_single", about = "N", function(d, ctx)
        if ("N" %in% ctx$nd && (!is.numeric(d[["N"]]) || length(d[["N"]]) != 1))
            "N must be a single number; coloc does not take per-snp sample sizes"),
    .rule("type_set", function(d, ctx)
        if (!("type" %in% ctx$nd))
            "variable type not set"),
    ## an s of any length but one failed here with R's error, as type did
    .rule("s_range", about = "s", function(d, ctx)
        if ("s" %in% ctx$nd) {
            if (length(d$s) != 1L)
                "s must be a single number between 0 and 1"
            else if (!is.numeric(d$s) || d$s <= 0 || d$s >= 1)
                "s must be between 0 and 1"
        }),
    ## the p-value route: beta and varbeta are needed to estimate var(Y) without them
    .rule("route", function(d, ctx)
        if (.p_route(ctx) && (!("pvalues" %in% ctx$nd) || !("MAF" %in% ctx$nd)))
            "require p values and MAF if beta, varbeta are unavailable"),
    ## character p-values are left to the numeric rule below: here they would
    ## be compared as strings
    .rule("p_range", about = "pvalues", function(d, ctx)
        if (.p_route(ctx) && is.numeric(d$pvalues) && any(d$pvalues <= 0 | d$pvalues > 1))
            "pvalues must be above 0 and at most 1"),
    .rule("s_needed", function(d, ctx)
        if (.p_route(ctx) && d$type == "cc" && !("s" %in% ctx$nd))
            "require, s, proportion of samples who are cases, if beta, varbeta are unavailable"),
    .rule("N_p_route", function(d, ctx)
        if (.p_route(ctx) && (!("N" %in% ctx$nd) || is.null(d$N) || any(d$N <= 0)))
            "sample size N <=0 or not set"),
    ## the values the route computes with, here where 6.0.3 worked out the
    ## p-values from beta and varbeta and failed on a character one with R's
    ## error. On the p-value route, route has made sure the p-values are there
    .rule("numeric", about = c("beta", "varbeta", "pvalues"), function(d, ctx) {
        for (v in if (.p_route(ctx)) "pvalues" else c("beta", "varbeta"))
            if (!is.numeric(d[[v]]))
                return(paste(v, "must be numeric"))
        NULL
    }),
    .rule("finite", function(d, ctx)
        if (any(is.infinite(d$beta) | is.infinite(d$varbeta)))
            "Infinite values in beta and/or varbeta"),
    .rule("zero_varbeta", function(d, ctx)
        if (any(d$varbeta == 0))
            "Zero values in varbeta"),
    ## after zero_varbeta, so a zero keeps its message. 6.0.3 failed on a
    ## negative one at the minimum p-value step, comparing a NaN, and a
    ## coloc_dataset, which skips that comparison, got NA posteriors
    .rule("varbeta_positive", about = "varbeta", function(d, ctx)
        if (!.p_route(ctx) && any(d$varbeta < 0))
            "varbeta must be positive"),
    ## A coloc_dataset gave this warning when it was made, so it is not
    ## repeated each time the dataset is used, and its p-values aren't worked
    ## out. A report unclasses the dataset, so it shows the minimum p-value
    .rule("min_p", level = "warning", function(d, ctx) {
        if (inherits(d, "coloc_dataset"))
            return(NULL)
        p <- min(if (.p_route(ctx)) d$pvalues else pnorm(-abs(d$beta / sqrt(d$varbeta))) * 2)
        if (p > ctx$warn.minp)
            paste0("minimum p value is: ", format.pval(p),
                   "\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument. If that's not the explanation, please check the 02_data vignette.")
        else .passed(paste("minimum p value is", format.pval(p)))
    }),
    .rule("sdY", function(d, ctx)
        if (d$type == "quant" && !("sdY" %in% ctx$nd))
            if (!("MAF" %in% ctx$nd && "N" %in% ctx$nd))
                "must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    ## the rules below need two dimensions. A vector failed in LD_square with
    ## R's error; a data frame has two, and passes where LD isn't required,
    ## as in 6.0.3
    .rule("LD_dims", about = "LD", function(d, ctx)
        if ("LD" %in% ctx$nd && length(dim(d$LD)) != 2L)
            "LD must be a matrix"),
    .rule("LD_square", about = "LD", function(d, ctx)
        if ("LD" %in% ctx$nd && nrow(d$LD) != ncol(d$LD))
            "LD not square"),
    .rule("LD_dimnames", about = "LD", function(d, ctx)
        if ("LD" %in% ctx$nd && !identical(colnames(d$LD), rownames(d$LD)))
            "LD rownames != colnames"),
    .rule("LD_covers", about = "LD", function(d, ctx)
        if ("LD" %in% ctx$nd && length(setdiff(d$snp, colnames(d$LD))))
            "colnames in LD do not contain all SNPs"),
    ## on the SuSiE route, which reads LD. After 6.0.3's rules, as its
    ## check_ld() refused an LD that wasn't a matrix after check_dataset();
    ## check_ld() keeps that test itself, as it runs these rules with no req
    .rule("LD_numeric", about = "LD", function(d, ctx)
        if ("LD" %in% ctx$nd && .susie_route(ctx) &&
            !(is.matrix(d$LD) && is.numeric(d$LD)))
            "LD must be a numeric matrix"),
    ## wherever N is given. After 6.0.3's rules, which took any N beside beta,
    ## so a dataset 6.0.3 refused keeps its message, and a dataset given as
    ## p-values "sample size N <=0 or not set"
    .rule("N_positive", about = "N", function(d, ctx)
        if ("N" %in% ctx$nd && (!is.finite(d$N) || d$N <= 0))
            "N must be a finite number above 0"),
    ## where sdY is read: a quantitative trait on the beta route, and not for
    ## SuSiE. After 6.0.3's rules, which took any sdY
    .rule("sdY_valid", about = "sdY", function(d, ctx)
        if ("sdY" %in% ctx$nd && !.p_route(ctx) && d$type == "quant" &&
            !.susie_route(ctx) &&
            (!is.numeric(d$sdY) || length(d$sdY) != 1L || !is.finite(d$sdY) || d$sdY <= 0))
            "sdY must be a single finite number above 0"),
    ## where z is read: SuSiE fine-maps from it in place of beta /
    ## sqrt(varbeta). After 6.0.3's rules, which didn't know z, so no_na and
    ## lengths leave it here
    .rule("z_valid", about = "z", function(d, ctx) {
        if (!("z" %in% ctx$nd && .susie_route(ctx)))
            return(NULL)
        if (any(is.na(d$z)))
            return("z contains missing values")
        ## z against snp, or the first other per-snp element without snp
        others <- setdiff(intersect(.per_snp_items, ctx$nd), "z")
        problem <- .lengths_differ(d, intersect(c(if ("snp" %in% others) "snp" else others[1], "z"),
                                                ctx$nd))
        if (!is.null(problem))
            return(problem)
        if (!is.numeric(d$z))
            return("z must be numeric")
        if (any(is.infinite(d$z)))
            "Infinite values in z"
    }),
    ## on the SuSiE route, which reads LD: coloc.abf() never reads it, and the
    ## test takes about 2 s at 10,000 snps. Last, after LD_numeric, which it
    ## needs. Only the dataset's snps, the part runsusie() uses, copied out
    ## only when LD covers more; with as many, LD_covers and LD_dimnames make
    ## them the same snps, and their order doesn't matter
    .rule("LD_correlation", about = "LD", function(d, ctx) {
        if (!("LD" %in% ctx$nd && .susie_route(ctx)))
            return(NULL)
        LD <- if (ncol(d$LD) > length(d$snp)) d$LD[d$snp, d$snp] else d$LD
        if (!.is_correlation(LD))
            "LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal"
    })
)
names(.dataset_rules) <- vapply(.dataset_rules, `[[`, "", "id")

## the rules that read one element, in list order
.rules_about <- function(element)
    Filter(function(rule) element %in% rule$about, .dataset_rules)

## "dataset {suffix}: " for a dataset checked with a label, such as 1 or 2 from
## coloc.abf() or a user's "eQTL"; nothing for one checked on its own (suffix
## "", NULL, NA or more than one value), as when a dataset is made. coloc()
## passes the word "trait", so its messages begin "trait {suffix}: "
.prefix <- function(suffix, word = "dataset") {
    if (length(suffix) != 1L || is.na(suffix) || identical(as.character(suffix), ""))
        return("")
    paste0(word, " ", suffix, ": ")
}

## Run the rules in order and stop at the first that fails. Every message, the
## warning's too, starts with .prefix(suffix, word), and carries `call`, the
## call of the function the user called. Returns invisible NULL.
.gate <- function(d, suffix = "", req = character(0), warn.minp = 1e-6,
                  call = NULL, rules = .dataset_rules, word = "dataset") {
    prefix <- .prefix(suffix, word)
    ## without a prefix, the one message with no subject gains one
    if (!is.list(d))
        stop(simpleError(if (nzchar(prefix)) paste0(prefix, "is not a list")
                         else "the dataset is not a list", call))
    ctx <- list(nd = intersect(names(d), .dataset_items), # defined in R/dataset.R
                req = req, warn.minp = warn.minp)
    for (rule in rules) {
        problem <- rule$check(d, ctx)
        if (is.null(problem) || inherits(problem, "passed"))
            next
        if (rule$level == "warning") {
            warning(simpleWarning(paste0(prefix, problem), call))
            next
        }
        stop(simpleError(paste0(prefix, problem), call))
    }
    invisible(NULL)
}

check_ld <- function(D,LD) {
    if(is.null(LD))
      stop("LD required")
    if(!is.matrix(LD))
      stop("LD must be of class matrix")
    ## if(any(LD[upper.tri(LD)]==1))
    ##   stop("LD includes SNPs in perfect LD")
    ## the rest are check_dataset()'s LD rules
    .gate(list(snp = D$snp, LD = LD), call = sys.call(), rules = .rules_about("LD"))
}

## ---- qc_dataset(): every rule, without stopping ------------------------------

## The rules only qc_dataset() runs: data coloc takes, but which looks wrong.
## Their ctx adds classed, whether the dataset as given was a coloc_dataset,
## and qc_dataset()'s thresholds
.report_rules <- list(
    .rule("coloc_dataset", level = "report", function(d, ctx)
        if (!ctx$classed)
            "a plain list: make it with coloc_dataset() or as_coloc_dataset()"),
    ## where sdY is read: a quantitative trait on the beta route. A varbeta
    ## holding standard errors gives an estimate far from the sdY given
    .rule("sdY_plausible", about = c("sdY", "varbeta", "MAF", "N"), level = "report",
          function(d, ctx) {
              if (d$type != "quant" || .p_route(ctx))
                  return(NULL)
              if (!is.numeric(d$sdY) || length(d$sdY) != 1L)
                  return(.not_checked("sdY isn't a single number"))
              fit <- .sdY_fit(d$varbeta, d$MAF, d$N)
              if (!is.finite(fit) || fit <= 0)
                  return("sdY can't be estimated from varbeta, MAF and N")
              estimate <- sqrt(fit)
              if (max(d$sdY / estimate, estimate / d$sdY) > ctx$sdY_ratio)
                  sprintf("sdY is %.3g, but varbeta, MAF and N suggest %.3g: does varbeta hold variances?",
                          d$sdY, estimate)
              else .passed(sprintf("sdY %.3g, estimated %.3g", d$sdY, estimate))
          }),
    ## rounding beta and varbeta moves beta / sqrt(varbeta) by well under 1%
    ## of |z| at 3 significant figures; a z for the other allele, or from
    ## another sample or model, moves it by far more
    .rule("z_matches", about = c("z", "beta", "varbeta"), level = "report",
          function(d, ctx) {
              gap <- abs(d$z - d$beta / sqrt(d$varbeta))
              n <- sum(gap > ctx$z_tol * pmax(1, abs(d$z)))
              if (n)
                  sprintf("z differs from beta / sqrt(varbeta) by more than %g%% at %d snp%s, by up to %.3g",
                          100 * ctx$z_tol, n, if (n == 1L) "" else "s", max(gap))
          }),
    ## only pairs where both snps have a signal: a z-score near 0 has a sign
    ## that says little, and where most snps have none, aligned data scores
    ## not far above half
    .rule("z_ld_signs", about = c("LD", "beta", "varbeta"), level = "report",
          function(d, ctx) {
              r <- .z_ld_ratios(.report_z(d, ctx), d$LD, d$snp, ctx$ld_thr, ctx$z_thr)
              if (length(r) < ctx$min_pairs)
                  return(.not_checked(sprintf(
                      "too few pairs to judge: %d with |LD| > %g and both |z| > %g",
                      length(r), ctx$ld_thr, ctx$z_thr)))
              agree <- mean(r > 0)
              if (agree < ctx$ld_agree)
                  sprintf("%.0f%% of %d pairs agree in sign (|LD| > %g, both |z| > %g): are beta and LD for the same alleles?",
                          100 * agree, length(r), ctx$ld_thr, ctx$z_thr)
              else .passed(sprintf("%.0f%% of %d pairs agree", 100 * agree, length(r)))
          })
)
names(.report_rules) <- vapply(.report_rules, `[[`, "", "id")

## the z-scores SuSiE reads: z when given, else beta / sqrt(varbeta)
.report_z <- function(d, ctx)
    if ("z" %in% ctx$nd) d$z else d$beta / sqrt(d$varbeta)

## a reason a report row couldn't be worked out, which a rule's check can
## return in place of a message
.not_checked <- function(why)
    structure(why, class = "not_checked")

## a pass, with what a report shows for it, which a rule's check can return
## in place of NULL; the gate takes it as a pass
.passed <- function(detail)
    structure(detail, class = "passed")

## One rule's row, as c(status, detail), never an error or a warning. A rule
## is "not given" when the dataset has none of its elements, or for a rule
## only the reports run, not all of them. A rule R fails on, as one after a
## rule that failed can, is "not checked", and so is one whose check gives
## anything but one sentence, which no row could show
.report_row <- function(rule, d, ctx) {
    if (length(rule$about)) {
        given <- rule$about %in% ctx$nd
        if (if (rule$level == "report") !all(given) else !any(given))
            return(c("not given", ""))
    }
    got <- tryCatch(
        withCallingHandlers(rule$check(d, ctx), warning = function(w) invokeRestart("muffleWarning")),
        error = function(e) .not_checked(paste("could not be checked:", conditionMessage(e))))
    if (is.null(got))
        return(c("ok", ""))
    if (!is.character(got) || length(got) != 1L || is.na(got))
        return(c("not checked", "could not be checked: no single message"))
    if (inherits(got, "not_checked"))
        return(c("not checked", unclass(got)))
    if (inherits(got, "passed"))
        return(c("ok", unclass(got)))
    ## the first line only, of the minimum p-value's three
    c(if (rule$level == "error") "FAIL" else "warn", sub("\n.*", "", got))
}

## every rule's row, in list order
.report <- function(d, rules, ctx) {
    rows <- vapply(rules, .report_row, character(2), d = d, ctx = ctx)
    .coloc_qc(names(rules), unname(rows[1, ]), unname(rows[2, ]))
}

.coloc_qc <- function(check, status, detail)
    structure(data.frame(check = check, status = status, detail = detail),
              class = c("coloc_qc", "data.frame"))

##' Every check of a dataset, without stopping
##'
##' `qc_dataset()` runs the rules [check_dataset()] runs, in the same order,
##' and lists every one, where `check_dataset()` stops at the first that
##' fails. It then adds four checks of data coloc takes but which looks wrong.
##' It takes any list, a plain one or a `coloc_dataset`, and never stops or
##' warns.
##'
##' Each row is one rule, named by its id, with one of five statuses:
##' * "ok": the rule found nothing.
##' * "FAIL": checking the dataset for `engine`, as `coloc()` does, would stop
##'   here, and the detail is its message.
##' * "warn": the minimum p-value, or one of the four checks below.
##' * "not given": the dataset hasn't the elements the rule is about.
##' * "not checked": the rule couldn't be worked out, usually because of a
##'   mistake a row above reports, or there were too few pairs for
##'   `z_ld_signs`.
##'
##' The checks only `qc_dataset()` makes:
##' * `coloc_dataset`: the dataset is a plain list, not made with
##'   [coloc_dataset()] or [as_coloc_dataset()].
##' * `sdY_plausible`: for a quantitative trait given as `beta` and `varbeta`,
##'   `sdY` and its estimate from `varbeta`, `MAF` and `N` are more than
##'   `sdY_ratio` apart, as when `varbeta` holds standard errors.
##' * `z_matches`: `z` differs from `beta / sqrt(varbeta)` by more than
##'   `z_tol` times `|z|`, or `z_tol` where `|z|` is below 1. Rounding `beta`
##'   and `varbeta` to 3 significant figures stays within the default.
##' * `z_ld_signs`: among the pairs of snps in absolute LD above `ld_thr`
##'   whose z-scores are both above `z_thr` in size, fewer than `ld_agree`
##'   agree in sign with their LD, as when `beta` and `LD` are for different
##'   alleles. It reads `z` when given, as SuSiE does. With fewer than
##'   `min_pairs` such pairs it is not checked. [check_alignment()] counts
##'   every pair.
##'
##' The checks on `LD` and `z` that SuSiE needs run with `engine = "susie"`,
##' as only SuSiE reads them. With "wakefield" they read "ok".
##' @param d a dataset: a list, or a `coloc_dataset`
##' @param engine the fine-mapper the dataset is for, "wakefield" or "susie",
##'   which decides the elements required, as `coloc()`'s `method` "single"
##'   or "susie" does
##' @param warn.minp the minimum p-value above which its row warns
##' @param sdY_ratio how far apart `sdY` and its estimate may be, as a ratio
##' @param ld_thr only pairs of snps in absolute LD above this count for
##'   `z_ld_signs`
##' @param z_thr and only pairs whose z-scores are both above this in size
##' @param ld_agree the share of those pairs that must agree in sign
##' @param min_pairs the fewest pairs `z_ld_signs` judges
##' @param z_tol how far `z` may be from `beta / sqrt(varbeta)`, as a share of
##'   `|z|`
##' @return a data frame of class `"coloc_qc"`, with columns `check`, `status`
##'   and `detail` and one row per rule. Printed, it shows a count and only
##'   the rows that are "FAIL", "warn" or "not checked", one line each;
##'   `print(x, all = TRUE)` shows every row. A part of it, such as
##'   `x[x$status == "FAIL", ]`, is a plain data frame.
##' @export
##' @examples
##' data(coloc_test_data)
##' qc_dataset(coloc_test_data$D1)
##' D <- as_coloc_dataset(coloc_test_data$D1)
##' D$varbeta <- sqrt(D$varbeta)     # standard errors where variances belong
##' qc_dataset(D)
qc_dataset <- function(d, engine = c("wakefield", "susie"), warn.minp = 1e-6,
                       sdY_ratio = 2, ld_thr = 0.2, z_thr = 3, ld_agree = 0.9,
                       min_pairs = 10, z_tol = 0.01) {
    engine <- match.arg(engine)
    if (!is.list(d))
        return(.coloc_qc("list", "FAIL", "the dataset is not a list"))
    ## "wakefield" is coloc()'s method "single"
    ctx <- list(nd = intersect(names(d), .dataset_items),
                req = .dataset_uses[[if (engine == "wakefield") "single" else engine]],
                warn.minp = warn.minp,
                classed = inherits(d, "coloc_dataset"), sdY_ratio = sdY_ratio,
                ld_thr = ld_thr, z_thr = z_thr, ld_agree = ld_agree,
                min_pairs = min_pairs, z_tol = z_tol)
    ## unclassed, so a coloc_dataset gets its minimum p-value row, which the
    ## gate skips for one
    .report(unclass(d), c(.dataset_rules, .report_rules), ctx)
}

##' @rdname qc_dataset
##' @param x a `coloc_qc`
##' @param all if TRUE, every row, not only those to look at
##' @param ... ignored
##' @export
print.coloc_qc <- function(x, all = FALSE, ...) {
    look <- x$status %in% c("FAIL", "warn", "not checked")
    ok <- sum(x$status == "ok")
    cat(ok, if (ok == 1L) " check ok, " else " checks ok, ",
        sum(x$status == "not given"), " not given; ",
        if (any(look)) paste0(sum(look), " to look at:") else "nothing to look at",
        "\n", sep = "")
    show <- all | look
    ## one line per row, as print.data.frame() would put a long detail in a
    ## block of its own under the check and status
    if (any(show))
        cat(paste0(" ", format(c("check", x$check[show])), " ",
                   format(c("status", x$status[show])), " ", c("detail", x$detail[show])),
            sep = "\n")
    invisible(x)
}

##' @rdname qc_dataset
##' @param i,j rows and columns, as for a data frame
##' @export
`[.coloc_qc` <- function(x, i, j, ...) {
    ## a part of a report is a plain data frame, as the count print.coloc_qc()
    ## gives would describe the part as if it were the whole
    out <- NextMethod()
    if (is.data.frame(out))
        class(out) <- "data.frame"
    out
}

## ---- qc_pair(): two datasets against each other ------------------------------

## element v of dataset d for the snps in shared, in that order
.at <- function(d, v, shared)
    d[[v]][match(shared, d[["snp"]])]

## The rules qc_pair() runs, on x = list(d1, d2, shared), shared being the snps
## in both. Their elements are those both datasets have
.pair_rules <- list(
    .rule("snps_shared", about = "snp", function(x, ctx)
        if (!length(x$shared))
            "no snp is in both datasets"
        else .passed(sprintf("%d of %d and %d", length(x$shared), length(x$d1[["snp"]]),
                             length(x$d2[["snp"]])))),
    .rule("positions_agree", about = c("snp", "position"), level = "report",
          function(x, ctx) {
              if (!length(x$shared))
                  return(.not_checked("no snp is in both datasets"))
              n <- sum(.at(x$d1, "position", x$shared) != .at(x$d2, "position", x$shared))
              if (n)
                  sprintf("%d shared snp%s in position: another genome build?",
                          n, if (n == 1L) " differs" else "s differ")
          }),
    ## MAF is the minor allele's frequency, so alleles flipped between the
    ## studies leave it as it is
    .rule("maf_agree", about = c("snp", "MAF"), level = "report",
          function(x, ctx) {
              if (!length(x$shared))
                  return(.not_checked("no snp is in both datasets"))
              n <- sum(abs(.at(x$d1, "MAF", x$shared) - .at(x$d2, "MAF", x$shared)) > ctx$maf_diff)
              if (n)
                  sprintf("%d shared snp%s in MAF by more than %g: different populations, or frequencies of different alleles?",
                          n, if (n == 1L) " differs" else "s differ", ctx$maf_diff)
          })
)
names(.pair_rules) <- vapply(.pair_rules, `[[`, "", "id")

##' Two datasets checked against each other
##'
##' `qc_pair()` checks what only shows when two datasets meet, before
##' [coloc()] compares them. [qc_dataset()] checks each dataset on its own.
##' Values are matched by snp id, so the datasets' orders don't matter.
##'
##' * `snps_shared`: "FAIL" when no snp is in both datasets, as `coloc()`
##'   works on the snps they share. Its detail gives the counts.
##' * `positions_agree`: how many shared snps have different positions, which
##'   usually means the two use different genome builds.
##' * `maf_agree`: how many shared snps have MAFs more than `maf_diff` apart.
##'   MAF is the minor allele's frequency, so flipped alleles don't change
##'   it; a difference means different populations, or that one dataset gives
##'   the frequency of a named allele as `MAF`.
##'
##' A check is "not given" unless both datasets have its elements. LD isn't
##' compared, as each study has its own.
##' @param d1,d2 datasets: lists, or `coloc_dataset`s
##' @param maf_diff how far apart two MAFs for the same snp may be
##' @return a `coloc_qc`, as [qc_dataset()] returns, printed the same way
##' @export
##' @examples
##' data(coloc_test_data)
##' qc_pair(coloc_test_data$D1, coloc_test_data$D2)
qc_pair <- function(d1, d2, maf_diff = 0.1) {
    if (!is.list(d1) || !is.list(d2))
        return(.coloc_qc("lists", "FAIL", "each dataset must be a list"))
    both <- intersect(names(d1), names(d2))
    ctx <- list(nd = intersect(both, .dataset_items), maf_diff = maf_diff)
    .report(list(d1 = d1, d2 = d2, shared = intersect(d1[["snp"]], d2[["snp"]])),
            .pair_rules, ctx)
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
  ## each pair of snps once, with LD matched to the snps by name
  tmp=.z_ld_ratios(D$beta/sqrt(D$varbeta), D$LD, D$snp, thr)
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
