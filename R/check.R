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
##' Optional vectors are
##'
##' \describe{
##'   \item{position}{a vector of snp positions, required for \code{plot_dataset}}
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
##' @param req names of elements that must be present
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
## test-check.R pins it against 6.0.3's. Each condition is the one
## check_dataset() had, word for word, so R's own errors on input no rule
## expects stay as they were.

## One rule. check is function(d, ctx): NULL when the dataset passes, else the
## message without its prefix, which .gate() adds. about names the elements
## the rule reads, so check_ld() can pick out the LD rules. level "warning" is
## the minimum p-value step alone: .gate() warns and goes on.
.rule <- function(id, check, about = character(0), level = c("error", "warning"))
    list(id = id, check = check, about = about, level = match.arg(level))

## the p-value route: no beta, or no varbeta
.p_route <- function(ctx) !("beta" %in% ctx$nd) || !("varbeta" %in% ctx$nd)

## ctx is list(nd, req, warn.minp), nd being the elements of d coloc knows
.dataset_rules <- list(
    .rule("required", function(d, ctx)
        if (length(setdiff(ctx$req, ctx$nd)))
            paste0("missing required element(s) ",
                   paste(setdiff(ctx$req, ctx$nd), collapse = ", "))),
    ## no missing values - make people clean their own data rather than make
    ## assumptions here for datasets I don't know
    .rule("no_na", function(d, ctx) {
        for (v in ctx$nd)
            if (any(is.na(d[[v]])))
                return(paste0(v, " contains missing values"))
        NULL
    }),
    ## guarded, so a dataset with no type reaches "variable type not set"
    ## below rather than failing here on a zero-length condition: runsusie()
    ## leaves type out of req
    .rule("type_valid", about = "type", function(d, ctx)
        if ("type" %in% ctx$nd && !(d$type %in% c("quant", "cc")))
            "type must be quant or cc"),
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
    ## the per-snp elements, compared in the dataset's own order
    .rule("lengths", function(d, ctx) {
        l <- -1 # impossible length
        for (v in intersect(ctx$nd, c("pvalues", "MAF", "beta", "varbeta", "snp", "position")))
            if (l < 0) {
                l <- length(d[[v]])
            } else if (length(d[[v]]) != l) {
                return("lengths of inputs don't match: ")
            }
        NULL
    }),
    .rule("N_single", about = "N", function(d, ctx)
        if ("N" %in% ctx$nd && (!is.numeric(d[["N"]]) || length(d[["N"]]) != 1))
            "N should be a single numeric. coloc does not accomodate per-snp sample sizes."),
    .rule("type_set", function(d, ctx)
        if (!("type" %in% ctx$nd))
            "variable type not set"),
    .rule("s_range", about = "s", function(d, ctx)
        if (("s" %in% ctx$nd) && (!is.numeric(d$s) || d$s <= 0 || d$s >= 1))
            "s must be between 0 and 1"),
    ## the p-value route: beta and varbeta are needed to estimate var(Y) without them
    .rule("route", function(d, ctx)
        if (.p_route(ctx) && (!("pvalues" %in% ctx$nd) || !("MAF" %in% ctx$nd)))
            "require p values and MAF if beta, varbeta are unavailable"),
    .rule("p_range", about = "pvalues", function(d, ctx)
        if (.p_route(ctx) && any(d$pvalues <= 0 | d$pvalues > 1))
            "pvalues must be above 0 and at most 1"),
    .rule("s_needed", function(d, ctx)
        if (.p_route(ctx) && d$type == "cc" && !("s" %in% ctx$nd))
            "require, s, proportion of samples who are cases, if beta, varbeta are unavailable"),
    .rule("N_p_route", function(d, ctx)
        if (.p_route(ctx) && (!("N" %in% ctx$nd) || is.null(d$N) || any(d$N <= 0)))
            "sample size N <=0 or not set"),
    .rule("finite", function(d, ctx)
        if (any(is.infinite(d$beta) | is.infinite(d$varbeta)))
            "Infinite values in beta and/or varbeta"),
    .rule("zero_varbeta", function(d, ctx)
        if (any(d$varbeta == 0))
            "Zero values in varbeta"),
    ## A coloc_dataset gave this warning when it was made, so it is not
    ## repeated each time the dataset is used. The p-values are worked
    ## out all the same, as before the rules were a list
    .rule("min_p", level = "warning", function(d, ctx) {
        p <- if (.p_route(ctx)) d$pvalues else pnorm(-abs(d$beta / sqrt(d$varbeta))) * 2
        if (!inherits(d, "coloc_dataset") && min(p) > ctx$warn.minp)
            paste0("minimum p value is: ", format.pval(min(p)),
                   "\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument. If that's not the explanation, please check the 02_data vignette.")
    }),
    .rule("sdY", function(d, ctx)
        if (d$type == "quant" && !("sdY" %in% ctx$nd))
            if (!("MAF" %in% ctx$nd && "N" %in% ctx$nd))
                "must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    .rule("LD_square", about = "LD", function(d, ctx)
        if ("LD" %in% ctx$nd && nrow(d$LD) != ncol(d$LD))
            "LD not square"),
    .rule("LD_dimnames", about = "LD", function(d, ctx)
        if ("LD" %in% ctx$nd && !identical(colnames(d$LD), rownames(d$LD)))
            "LD rownames != colnames"),
    .rule("LD_covers", about = "LD", function(d, ctx)
        if ("LD" %in% ctx$nd && length(setdiff(d$snp, colnames(d$LD))))
            "colnames in LD do not contain all SNPs")
)
names(.dataset_rules) <- vapply(.dataset_rules, `[[`, "", "id")

## the rules that read one element, in list order
.rules_about <- function(element)
    Filter(function(rule) element %in% rule$about, .dataset_rules)

## "dataset {suffix}: " for a dataset checked with a label, such as 1 or 2 from
## coloc() or a user's "eQTL"; nothing for one checked on its own (suffix "",
## NULL, NA or more than one value), as when a dataset is made
.prefix <- function(suffix) {
    if (length(suffix) != 1L || is.na(suffix) || identical(as.character(suffix), ""))
        return("")
    paste0("dataset ", suffix, ": ")
}

## Run the rules in order and stop at the first that fails. Every message, the
## warning's too, starts with .prefix(suffix), and carries `call`, the call of
## the function the user called. Returns invisible NULL.
.gate <- function(d, suffix = "", req = character(0), warn.minp = 1e-6,
                  call = NULL, rules = .dataset_rules) {
    prefix <- .prefix(suffix)
    if (!is.list(d))
        stop(simpleError(paste0(prefix, "is not a list"), call))
    ctx <- list(nd = intersect(names(d), .dataset_items), # defined in R/dataset.R
                req = req, warn.minp = warn.minp)
    for (rule in rules) {
        problem <- rule$check(d, ctx)
        if (is.null(problem))
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

##' check alignment between beta and LD
##'
##' @title check alignment
##' @param D a coloc dataset
##' @param thr plot SNP pairs in absolute LD > thr
##' @param do_plot if TRUE (default) plot the diagnostic
##' @export
##' @return proportion of pairs that are positive
##' @author Chris Wallace
##' @examples
##' data(coloc_test_data)
##' ## close to 1 when beta and LD are aligned to the same alleles
##' check_alignment(coloc_test_data$D3, do_plot = FALSE)
check_alignment <- function(D,thr=0.2,do_plot=TRUE) {
  check_dataset(D)
  bprod=outer(D$beta/sqrt(D$varbeta),D$beta/sqrt(D$varbeta),"*")
  ## plot(bprod,D$LD[D$snp,D$snp],xlab="product of z scores",ylab="LD")
  tmp=(bprod/D$LD)[abs(D$LD) > 0.2]
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
