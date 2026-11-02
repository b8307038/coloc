## ---------------------------------------------------------------------------
## coloc(), the entry point. It checks each input, turns it into log Bayes
## factors with the adapter for its kind, and has the engine compare the two.
## Below it are the functions that read its result, and the helpers the old
## functions use to call it.
## ---------------------------------------------------------------------------

##' Colocalisation of two genetic traits
##'
##' Asks whether two traits share a causal variant in a region. `coloc()` is
##' the one entry point: `x1` and `x2` may each be a coloc dataset, a SuSiE
##' fit, or log Bayes factors, in any combination.
##'
##' @section Inputs:
##' Each of `x1` and `x2` is one of three things.
##'
##' * A coloc dataset, which holds one trait's summary statistics for the
##'   region. It is a list with vectors that have one value per SNP (`snp`,
##'   and `beta` with `varbeta`, or `pvalues` with `MAF`) and single values
##'   that describe the study (`type`, which is `"quant"` or `"cc"`; `N`; and
##'   `sdY` or `s`). Fine-mapping with SuSiE also needs `LD`, the correlation
##'   matrix of the SNPs, and `beta` with `varbeta`, as p-values are not
##'   enough. [check_dataset()] lists every element and which combinations
##'   are needed. Make a dataset with [coloc_dataset()], or with
##'   [as_coloc_dataset()] from a list you already have; both check it as
##'   they make it. [coloc_test_data] holds four examples. A plain list works
##'   too, as `coloc()` converts it and says so once per session.
##' * A SuSiE fit, the result of [runsusie()] on a coloc dataset, or of
##'   `susieR::susie_rss()` passed through [annotate_susie()], which names its
##'   SNPs. `coloc()` uses the log Bayes factors of the credible sets the fit
##'   found, one signal per credible set.
##' * Log Bayes factors, as a numeric vector with one value per SNP, named by
##'   SNP id, for a single signal; or as a matrix with one row per signal and
##'   one column per SNP, with SNP ids as its column names. A column named
##'   `"null"`, if there is one, holds the log Bayes factor for no
##'   association. They must have no missing or infinite values, and no SNP
##'   id twice.
##'
##' `coloc()` checks each dataset with [check_dataset()], for what its
##' `method` needs, before it fine-maps either trait. It stops at the first
##' problem, with a message that begins "dataset 1:" or "dataset 2:", so a
##' mistake in trait 2 stops `coloc()` before SuSiE runs on trait 1.
##'
##' The two traits are matched by SNP id, and only the SNPs both have are
##' compared.
##'
##' @section Choosing `method`:
##' `method` says how a dataset is turned into log Bayes factors, and so how
##' many signals a trait can have.
##'
##' * `"single"` assumes the trait has at most one causal variant in the
##'   region, and computes one Bayes factor per SNP from `beta` and `varbeta`
##'   (or `pvalues` and `MAF`), as [coloc.abf()] did. This is the default.
##' * `"susie"` fine-maps the dataset with SuSiE through [runsusie()], which
##'   needs `beta` and `varbeta` (p-values are not enough), `LD` and `N`, and
##'   gives one signal per credible set, so a trait can have several.
##'
##' One value applies to both traits. Two values apply to trait 1 and trait 2
##' in turn, so the traits can be treated differently: `method = c("single",
##' "susie")` fine-maps only trait 2 with SuSiE, for example when only trait 2
##' has an LD matrix. `method` applies only to datasets, because a SuSiE fit
##' is already fine-mapped and log Bayes factors are already computed. Asking
##' for `"single"` for a SuSiE fit is an error, and naming a method for a trait
##' given as log Bayes factors gives a warning.
##'
##' @param x1 trait 1: a coloc dataset, a SuSiE fit or log Bayes factors (see
##'   Inputs)
##' @param x2 trait 2, with the same choices
##' @param method `"single"` or `"susie"`: one value for both traits, or one
##'   for each trait (see Choosing `method`). The default, `c("single",
##'   "susie")`, lists the choices: if `method` is not given, both traits use
##'   `"single"`
##' @param prior the prior probabilities, made with [coloc_prior()]
##' @param control when to leave a pair of signals out, made with
##'   [coloc_control()]
##' @param susie.args other arguments for [runsusie()], for datasets
##'   fine-mapped here
##' @return A `coloc_result`. Its `$summary` has one row per pair of signals
##'   compared, with the posterior probabilities of the five hypotheses,
##'   `PP.H0.abf` to `PP.H4.abf`; `$results` gives, for each SNP, the
##'   posterior probability that it is the shared causal variant; and
##'   `$status` says what happened. [coloc_result] describes every field.
##' @seealso [coloc_prior()], [coloc_control()], [pp()], [is_ok()],
##'   [n_signals()], [log_bf()], [coloc_result]
##' @export
##' @examples
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##' D2 <- as_coloc_dataset(coloc_test_data$D2)
##'
##' ## two datasets, with at most one causal variant assumed in each trait
##' res <- coloc(D1, D2)
##' res
##'
##' ## reading the result
##' is_ok(res)              # the posteriors were computed
##' n_signals(res)          # pairs of signals compared
##' pp(res, row = 1)        # number of snps, and PP.H0 to PP.H4 for the first pair
##' bf1 <- log_bf(res, 1)   # trait 1's log Bayes factors, to fine-map
##' finemap.bf(bf1)$summary
##'
##' if (requireNamespace("susieR", quietly = TRUE)) {
##'     D3 <- as_coloc_dataset(coloc_test_data$D3)
##'     D4 <- as_coloc_dataset(coloc_test_data$D4)
##'
##'     ## fine-map both traits with SuSiE: one row per pair of credible sets
##'     res2 <- coloc(D3, D4, method = "susie")
##'     res2$summary
##'     pp(res2, row = 2)   # the second pair: row picks a row of $summary
##'
##'     ## fine-map only trait 2 with SuSiE. Trait 1 is one signal and trait 2
##'     ## has two credible sets, so trait 1 is compared with each: two rows
##'     res3 <- coloc(D1, D3, method = c("single", "susie"))
##'     res3$summary
##'
##'     ## two matrices of log Bayes factors, one row per signal: here, the
##'     ## credible sets of two SuSiE fits
##'     fit3 <- runsusie(D3)
##'     fit4 <- runsusie(D4)
##'     bf3 <- fit3$lbf_variable[fit3$sets$cs_index, , drop = FALSE]
##'     bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
##'     res4 <- coloc(bf3, bf4)
##'     res4$summary
##' }
coloc <- function(x1, x2,
                  method = c("single", "susie"),
                  prior = coloc_prior(),
                  control = coloc_control(),
                  susie.args = list()) {
    ## method's default lists the choices, c("single", "susie"), which would
    ## read as "single" for trait 1 and "susie" for trait 2, so a method not
    ## given is "single"
    supplied <- !missing(method)
    if (!supplied) method <- "single"

    ## a coloc_prior too, as one edited with $<- keeps its class
    validate_coloc_prior(prior)

    kinds  <- c(.classify_input(x1), .classify_input(x2))
    method <- .resolve_method(method, kinds, supplied = supplied)

    ## Each dataset is checked, with the trait's number and what its method
    ## needs, before either is fine-mapped, and comes back a coloc_dataset. A
    ## mistake in trait 2 then stops coloc() before SuSiE runs on trait 1, and
    ## the checks inside wakefield_bf() and runsusie() find nothing more.
    xs <- list(x1, x2)
    for (i in 1:2)
        if (kinds[i] == "dataset")
            xs[[i]] <- .new_coloc_dataset(xs[[i]], trait = i, use = method[i])

    ## the adapters check what only they can: a fit's or log BFs' values
    adapt <- function(x, kind, m, trait) switch(kind,
        susie_fit = .adapt_susie(x, trait, susie.args = susie.args),
        bf        = .adapt_bf(x, trait),
        dataset   = if (identical(m, "susie")) .adapt_susie(x, trait, susie.args = susie.args)
                    else .adapt_single(x))
    a1 <- adapt(xs[[1]], kinds[1], method[1], 1L)
    a2 <- adapt(xs[[2]], kinds[2], method[2], 2L)

    ## The prior for one snp times the number of snps must stay below 1, or
    ## nothing is left for "no causal variant". adjust_prior() warns and caps
    ## it at 1/(nsnps + 1), as coloc.abf() does; p12 counts the shared snps.
    snps1 <- setdiff(colnames(a1$bf), "null")
    snps2 <- setdiff(colnames(a2$bf), "null")
    prior$p1  <- adjust_prior(prior$p1,  length(snps1), "1")
    prior$p2  <- adjust_prior(prior$p2,  length(snps2), "2")
    prior$p12 <- adjust_prior(prior$p12, length(intersect(snps1, snps2)), "12")

    res <- .coloc_engine(a1$bf, a2$bf, prior = prior, control = control)

    ## the engine numbers each trait's signals by their rows in the adapter's
    ## matrix; meta$idx maps them back to the input's own numbers, such as a
    ## SuSiE fit's components
    if (nrow(res$summary)) {
        res$summary$idx1 <- a1$meta$idx[res$summary$idx1]
        res$summary$idx2 <- a2$meta$idx[res$summary$idx2]
    }
    ## the engine doesn't know which adapter ran, so coloc() adds the methods
    ## and each trait's signal ids to $meta
    res$meta$method <- c(a1$meta$method, a2$meta$method)
    res$meta$idx1 <- a1$meta$idx
    res$meta$idx2 <- a2$meta$idx
    ## each trait's log BFs as the engine got them, before they were cut to the
    ## shared snps, so they can be fine-mapped without computing them again
    res$bf <- list(bf1 = a1$bf, bf2 = a2$bf)
    validate_coloc_result(res)

    ## a plain list that got this far passed its checks and was converted
    plain <- which(vapply(list(x1, x2), .is_plain_list, NA))
    if (length(plain))
        .note_once("plain_list", .plain_list_note("coloc", plain))
    res
}

## a list given as a dataset that isn't a coloc_dataset; a SuSiE fit is a list too
.is_plain_list <- function(x)
    is.list(x) && !inherits(x, c("coloc_dataset", "susie"))

## ---- functions that read a coloc_result -----------------------------------

##' Posterior probabilities from a coloc_result
##'
##' Returns the named vector that `coloc.abf()`'s `$summary` used to be, for
##' one row of the result.
##'
##' @param x a coloc_result
##' @param row which signal pair, from 1 to `n_signals(x)`
##' @return named numeric, length 6
##' @export
##' @examples
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##' D2 <- as_coloc_dataset(coloc_test_data$D2)
##' res <- coloc(D1, D2)
##' pp(res)                   # row = 1, the first pair, by default
##'
##' if (requireNamespace("susieR", quietly = TRUE)) {
##'     ## two credible sets in D3 against one in D4: two pairs
##'     D3 <- as_coloc_dataset(coloc_test_data$D3)
##'     D4 <- as_coloc_dataset(coloc_test_data$D4)
##'     res2 <- coloc(D3, D4, method = "susie")
##'     n_signals(res2)
##'     pp(res2, row = 2)
##' }
pp <- function(x, row = 1) {
    .need_coloc_result(x)
    ## with no rows, nsnps is 0, as n_signals() is, and the posteriors are NA,
    ## as none were computed
    if (!nrow(x$summary))
        return(stats::setNames(c(0, rep(NA_real_, 5L)),
                               c("nsnps", paste0("PP.H", 0:4, ".abf"))))
    ## data.table would return a row of NA for a row past the end, and nothing
    ## for row 0
    if (length(row) != 1 || !row %in% seq_len(nrow(x$summary)))
        stop("row must be between 1 and ", nrow(x$summary), ", the number of signal pairs")
    unlist(x$summary[row, c("nsnps", paste0("PP.H", 0:4, ".abf")), with = FALSE])
}

##' Did a colocalisation analysis produce posteriors for every pair?
##'
##' `is_ok()` is `TRUE` when the status is `"ok"`: at least one pair of
##' signals was compared, and every pair in `$summary` has its posterior
##' probabilities. Pairs left out under `overlap.min`, while others were kept,
##' have no row and do not count against it. It answers for the whole result,
##' so with status `"degenerate"` it is `FALSE` even if some pairs have
##' posteriors. It says nothing about whether the traits colocalise: read
##' the posteriors for that.
##'
##' `if (res$summary$PP.H4.abf > 0.9)` fails when the posteriors were not
##' computed. With no rows, when the traits share no SNPs or a trait has no
##' signal, the test is `logical(0)`. When every pair was trimmed, or a pair's
##' posteriors could not be computed, `PP.H0.abf` to `PP.H4.abf` are `NA`, and
##' so is the test. Ask `is_ok()` first.
##'
##' @param x a coloc_result
##' @return `TRUE` or `FALSE`
##' @export
##' @examples
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##' D2 <- as_coloc_dataset(coloc_test_data$D2)
##' res <- coloc(D1, D2)
##' if (is_ok(res) && pp(res)["PP.H4.abf"] > 0.8)
##'     message("the traits share a causal variant")
is_ok <- function(x) {
    .need_coloc_result(x)
    identical(x$status, "ok")
}

##' How many signal pairs a coloc_result carries
##'
##' @param x a coloc_result
##' @return integer(1)
##' @export
##' @examples
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##' D2 <- as_coloc_dataset(coloc_test_data$D2)
##' res <- coloc(D1, D2)
##' n_signals(res)
n_signals <- function(x) {
    .need_coloc_result(x)
    nrow(x$summary)
}

##' Log Bayes factors a coloc_result was computed from
##'
##' Each trait's log Bayes factors, as they went into the colocalisation:
##' signals in rows, snps in columns, covering the trait's own snps before
##' the two traits were intersected or any signal pair trimmed. Rows are named
##' as in `x$meta$idx1` and `x$meta$idx2`, unless log Bayes factors given to
##' [coloc()] had row names of their own. Use them to fine-map a trait from
##' the same numbers [coloc()] used.
##'
##' @param x a coloc_result
##' @param trait 1 or 2
##' @return numeric matrix, signals x snps
##' @export
##' @examples
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##' D2 <- as_coloc_dataset(coloc_test_data$D2)
##' res <- coloc(D1, D2)
##' bf1 <- log_bf(res, 1)
##' dim(bf1)                         # one signal, 500 SNPs
##' ## fine-map trait 1 from the same numbers coloc() used
##' finemap.bf(bf1)$summary
log_bf <- function(x, trait = 1) {
    .need_coloc_result(x)
    if (length(trait) != 1 || !trait %in% 1:2)
        stop("trait must be 1 or 2")
    if (is.null(x$bf)) {
        if (isTRUE(x$meta$schema < 2))
            stop("this coloc_result was made before coloc() kept the log Bayes ",
                 "factors; rerun coloc() to get them")
        stop("this coloc_result carries no log Bayes factors")
    }
    x$bf[[paste0("bf", trait)]]
}

##' @export
print.coloc_result <- function(x, ...) {
    cat("coloc_result: status", x$status, "with", nrow(x$summary), "signal pair(s)\n")
    if (!is.null(x$meta$method))
        cat("  method:", paste(ifelse(is.na(x$meta$method), "bf", x$meta$method),
                               collapse = ", "), "\n")
    cat("  priors:", paste(names(x$priors), signif(x$priors, 3),
                           sep = " = ", collapse = ", "), "\n")
    if (nrow(x$summary)) print(x$summary)
    for (m in x$messages) cat("  note:", m, "\n")
    invisible(x)
}

## ---- helpers for the old functions ----------------------------------------

## Builds the prior from the old functions' loose arguments in `...`, starting
## from the defaults. It is a plain list, which coloc() checks.
.prior_from_dots <- function(dots) {
    out <- .default_prior()
    for (nm in c("p1", "p2", "p12"))
        if (!is.null(dots[[nm]])) out[[nm]] <- dots[[nm]]
    ## x["nm"] <- list(NULL) keeps the element, as NULL, where
    ## x[["nm"]] <- NULL would delete it, so the prior has .default_prior()'s
    ## elements
    out["weights1"] <- list(dots$prior_weights1)
    out["weights2"] <- list(dots$prior_weights2)
    out
}

.control_from_dots <- function(dots) {
    out <- .default_control()
    for (nm in c("overlap.min", "trim_by_posterior"))
        if (!is.null(dots[[nm]])) out[[nm]] <- dots[[nm]]
    out
}

## Runs coloc() for an old function, without coloc()'s note on a plain list or
## its warning that method doesn't apply to log BFs, as neither is about a
## choice the old function's caller had. The note's key is put back, so a
## user's own coloc() on a plain list later in the session still gets it.
.as_old_function <- function(expr) {
    noted <- .coloc_state$noted
    on.exit(.coloc_state$noted <- noted)
    withCallingHandlers(expr,
                        coloc_plain_list = function(w) invokeRestart("muffleWarning"),
                        coloc_method_ignored = function(w) invokeRestart("muffleWarning"))
}
