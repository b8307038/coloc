## ---------------------------------------------------------------------------
## Adapters: each input becomes a matrix of log Bayes factors.
##
## Each adapter returns
##
##   list(bf   = <signals x snps, dimnames(signal_id, snp)>,
##        meta = list(method, idx, has_null))
##
## so a new fine-mapping method needs one adapter and one `method` value, and
## the engine, which sees only the matrices, doesn't change. .classify_input()
## tells what kind of input each trait is, and .resolve_method() picks the
## method for each dataset.
## ---------------------------------------------------------------------------

## ---- what kind of input each trait is -------------------------------------

## Tells the kind of input from its class alone. Whether its contents can be
## used is checked later, by new_coloc_dataset() for a dataset and by the
## adapter for log BFs, with messages that name the problem, so anything that
## is neither a SuSiE fit nor log BFs is taken as a dataset.
.classify_input <- function(x) {
    if (is.null(x))                    return("dataset")  # is.atomic(NULL) is TRUE before R 4.4
    if (inherits(x, "susie"))          return("susie_fit")
    if (is.matrix(x))                  return("bf")
    if (is.atomic(x))                  return("bf")   ## type checked by .adapt_bf()
    "dataset"
}

## ---- the method for each trait --------------------------------------------

.resolve_method <- function(method, kinds, supplied = TRUE) {
    method <- unname(rep_len(method, 2L))        ## length 1 recycles to both traits
    ## without USE.NAMES = FALSE, vapply() names each value after itself, and
    ## identical(m, "susie") is then FALSE
    method <- vapply(method, match.arg, character(1L),
                     choices = c("single", "susie"), USE.NAMES = FALSE)
    for (i in seq_len(2L)) {
        ## a fit implies "susie", so "single" is an error only when the caller
        ## named it, not when it is the default
        if (kinds[i] == "susie_fit" && method[i] == "single") {
            if (supplied)
                stop("trait ", i, ': method = "single" on a susie fit; a fit carries no ',
                     "beta/varbeta to run the Wakefield calculation on")
            method[i] <- "susie"
        }
        if (kinds[i] == "bf") {
            if (supplied)
                warning(warningCondition(
                    paste0("trait ", i, ": method is inapplicable to precomputed log BFs"),
                    class = "coloc_method_ignored", call = NULL))
            method[i] <- NA_character_
        }
    }
    method
}

## ---- the three adapters ----------------------------------------------------

## Makes Wakefield's log ABF for each snp with wakefield_bf(), as coloc.abf()
## computes them. The matrix has one row, for the one causal variant the
## method assumes, and no null column, as coloc.abf() has none.
.adapt_single <- function(x) {
    w <- wakefield_bf(x)
    bf <- matrix(w$lABF, nrow = 1L, dimnames = list("sig1", w$snp))
    list(bf = bf, meta = list(method = "single", idx = 1L, has_null = FALSE))
}

## For a SuSiE fit, or a dataset fine-mapped with runsusie(), takes the rows of
## lbf_variable that form credible sets, named by $sets$cs_index
.adapt_susie <- function(x, trait = 1L, susie.args = list(), ...) {
    if (inherits(x, "susie")) {
        fit <- x
    } else {
        ## coloc() sets runsusie()'s suffix to the trait's number
        if ("suffix" %in% names(susie.args))
            stop("susie.args cannot set suffix: coloc() sets it to each trait's number")
        ## the dataset goes in by its name, x, so an error from runsusie()
        ## doesn't carry a call holding the whole dataset
        fit <- do.call("runsusie",
                       c(list(d = quote(x), suffix = trait), susie.args))
    }
    if (is.null(fit$lbf_variable))
        stop("trait ", trait, ": the susie fit has no lbf_variable")
    if (is.null(colnames(fit$lbf_variable)))
        stop("trait ", trait, ": the susie fit's lbf_variable needs snp ids as ",
             "colnames, which runsusie() adds")
    idx <- fit$sets$cs_index
    if (is.null(idx) || !length(idx)) {          ## -> engine status "no_signal"
        snps <- colnames(fit$lbf_variable)
        return(list(bf = matrix(numeric(0), nrow = 0L, ncol = length(snps),
                                dimnames = list(NULL, snps)),
                    meta = list(method = "susie", idx = integer(0),
                                ## the fit's colnames are kept above, so this
                                ## branch can carry a "null" column too
                                has_null = "null" %in% snps)))
    }
    bf <- fit$lbf_variable[idx, , drop = FALSE]
    ## a fit given directly was never checked, so the rows the engine gets are
    ## checked here, as log BFs given directly are
    .check_bf_values(bf, trait)
    rownames(bf) <- paste0("sig", idx)           ## lbf_variable carries no rownames
    list(bf = bf, meta = list(method = "susie", idx = idx,
                              has_null = "null" %in% colnames(bf)))
}

## Checks the values of log BFs, for .adapt_bf() and .adapt_susie() alike: no
## missing or infinite values, as a missing one makes every pair using its row
## NA, and no snp id twice, as the engine matches snps by name and would use
## only the first
.check_bf_values <- function(x, trait) {
    if (any(!is.finite(x)))
        stop("trait ", trait, ": log Bayes factors contain missing or infinite values")
    if (anyDuplicated(colnames(x)))
        stop("trait ", trait, ": duplicated snps found in the log Bayes factors")
    invisible(x)
}

## Checks log BFs given directly; every row is a signal
.adapt_bf <- function(x, trait = 1L, ...) {
    ## first, as a character matrix would pass the checks below and fail in
    ## logbf_to_pp() with R's "non-numeric argument to binary operator"
    if (!is.numeric(x))
        stop("trait ", trait, ": log Bayes factors must be numeric, not ", typeof(x))
    if (is.null(dim(x))) {
        if (is.null(names(x)))
            stop("trait ", trait, ": a log-BF vector must be named by snp")
        x <- matrix(x, nrow = 1L, dimnames = list("sig1", names(x)))
    }
    if (is.null(colnames(x)))
        stop("trait ", trait, ": a log-BF matrix needs snp ids as colnames")
    .check_bf_values(x, trait)
    ## a matrix with no rows gets no names: assigning character(0) is an error
    if (is.null(rownames(x)) && nrow(x)) rownames(x) <- paste0("sig", seq_len(nrow(x)))
    list(bf = x, meta = list(method = NA_character_, idx = seq_len(nrow(x)),
                             has_null = "null" %in% colnames(x)))
}
