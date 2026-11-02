## ---------------------------------------------------------------------------
## The colocalisation engine. It takes two matrices of log Bayes factors, one
## per trait, with signals in rows and snps in columns, compares each signal
## of trait 1 with each of trait 2, and returns a coloc_result, whose $status
## says what happened. It sees only the matrices, not how they were made, so
## each trait's log Bayes factors can come from a different method.
## ---------------------------------------------------------------------------

## The defaults as plain lists, from coloc_prior() and coloc_control(), so
## each default is written in one place
.default_prior   <- function() unclass(coloc_prior())
.default_control <- function() unclass(coloc_control())

.coloc_engine <- function(bf1, bf2,
                          prior = .default_prior(),
                          control = .default_control()) {

    msgs <- character(0)
    ## each note is a warning, and is kept in $messages too, so a loop over
    ## many regions can read what happened without catching warnings
    note <- function(...) {
        m <- paste0(...)
        warning(m, call. = FALSE)
        msgs <<- c(msgs, m)
    }

    p1 <- prior$p1
    p2 <- prior$p2
    p12 <- prior$p12
    w1 <- prior$weights1
    w2 <- prior$weights2

    ## ---- 0. two matrices with snp ids, and weights to match them ----------
    stopifnot(is.matrix(bf1), is.matrix(bf2),
              !is.null(colnames(bf1)), !is.null(colnames(bf2)))
    if (!is.null(w1) && ncol(bf1) != length(w1))
        stop("Length of prior_weights1 must match size of dataset 1")
    if (!is.null(w2) && ncol(bf2) != length(w2))
        stop("Length of prior_weights2 must match size of dataset 2")

    ## ---- 1. a trait with no signal, as when SuSiE finds no credible set ---
    if (nrow(bf1) == 0L || nrow(bf2) == 0L) {
        note(if (nrow(bf1) == 0L && nrow(bf2) == 0L) "traits 1 and 2 have"
             else paste("trait", if (nrow(bf1) == 0L) 1 else 2, "has"),
             " no signal, so there is nothing to colocalise; with SuSiE, this means no credible set")
        return(new_coloc_result(status = "no_signal", messages = msgs, prior = prior))
    }

    ## ---- 2. the snps both traits have, leaving out the null column --------
    isnps <- setdiff(intersect(colnames(bf1), colnames(bf2)), "null")
    if (!length(isnps)) {
        note("no snps in common between the two traits")
        return(new_coloc_result(status = "no_overlap", messages = msgs, prior = prior))
    }

    ## ---- 3. the null column, where a matrix has one -----------------------
    ## moved last, and subtracted from every column, so the null's log BF is 0
    nulllast <- function(m) m[, c(setdiff(colnames(m), "null"), "null"), drop = FALSE]
    if ("null" %in% colnames(bf1))
        bf1 <- nulllast(bf1) - matrix(bf1[, "null"], nrow(bf1), ncol(bf1))
    if ("null" %in% colnames(bf2))
        bf2 <- nulllast(bf2) - matrix(bf2[, "null"], nrow(bf2), ncol(bf2))

    ## ---- 4. how much of each signal's posterior is on the shared snps -----
    pp1 <- logbf_to_pp(bf1, p1, last_is_null = "null" %in% colnames(bf1))
    pp2 <- logbf_to_pp(bf2, p2, last_is_null = "null" %in% colnames(bf2))
    prop <- function(pp) rowSums(pp[, isnps, drop = FALSE]) /
                         rowSums(pp[, setdiff(colnames(pp), "null"), drop = FALSE])
    prop1 <- prop(pp1)
    prop2 <- prop(pp2)

    ## ---- 5. every pair of signals, and those trimmed under overlap.min ----
    todo <- as.data.table(expand.grid(i = seq_len(nrow(bf1)), j = seq_len(nrow(bf2))))
    n_in <- nrow(todo)
    drop <- if (isTRUE(control$trim_by_posterior))
                prop1[todo$i] < control$overlap.min | prop2[todo$j] < control$overlap.min
            else rep(FALSE, n_in)
    ## a share is NaN when every snp posterior of a signal underflows to 0,
    ## from complete input; keep such a pair rather than trim on it
    drop[is.na(drop)] <- FALSE

    if (all(drop)) {
        note("every signal pair was left out, because less than overlap.min (",
             control$overlap.min, ") of a signal's posterior falls on snps both traits have")
        ## each trait's own top snp by log BF, shared or not, which shows where
        ## its signal is when no pair overlaps enough
        top <- function(bf, ix) {
            snps <- setdiff(colnames(bf), "null")
            snps[apply(bf[, snps, drop = FALSE], 1, which.max)][ix]
        }
        skel <- data.table(nsnps = length(isnps),
                           hit1 = top(bf1, todo$i), hit2 = top(bf2, todo$j),
                           PP.H0.abf = NA_real_,
                           PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
                           PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
                           idx1 = todo$i, idx2 = todo$j)
        return(new_coloc_result(summary = skel, status = "all_trimmed",
                                messages = msgs, prior = prior,
                                meta = list(n_pairs_in = n_in, n_pairs_kept = 0L,
                                            dropped = "overlap.min")))
    }
    if (any(drop))
        todo <- todo[!drop]

    ## ---- 6. each matrix and its weights, cut to the shared snps -----------
    restrict <- function(bf, w) {
        if (length(isnps) == ncol(bf)) return(list(bf = bf, w = w))
        keep <- match(isnps, colnames(bf))
        list(bf = bf[, keep, drop = FALSE], w = w[keep])
    }
    r1 <- restrict(bf1, w1)
    bf1 <- r1$bf
    w1 <- r1$w
    r2 <- restrict(bf2, w2)
    bf2 <- r2$bf
    w2 <- r2$w

    ## combine_abf_weighted() is used only with weights, as it works out H3 in
    ## another order, which would change unweighted results in the last digits.
    ## Both need single numbers for p1, p2 and p12. coloc() and coloc.bf_bf()
    ## check that, and coloc.abf() stops on a vector in adjust_prior().
    user_weighted <- !is.null(w1) || !is.null(w2)

    ## ---- 7. one signal pair at a time -------------------------------------
    hit_of <- function(bf, i) {
        w <- which.max(bf[i, ])
        if (!length(w)) NA_character_ else colnames(bf)[w]
    }
    pps  <- matrix(NA_real_, nrow(todo), 5L)
    snpp <- vector("list", nrow(todo))
    h1 <- character(nrow(todo))
    h2 <- character(nrow(todo))
    degenerate <- FALSE

    for (k in seq_len(nrow(todo))) {
        l1 <- bf1[todo$i[k], isnps]
        l2 <- bf2[todo$j[k], isnps]
        lsum <- l1 + l2
        sp <- exp(lsum - logsum(lsum))
        pp <- if (user_weighted)
                  combine_abf_weighted(l1, l2, p1, p2, p12, w1, w2, quiet = TRUE)
              else
                  combine.abf(l1, l2, p1, p2, p12, quiet = TRUE)
        h1[k] <- hit_of(bf1, todo$i[k])
        h2[k] <- hit_of(bf2, todo$j[k])
        ## when a pair's posteriors can't be computed, as from an NA log BF,
        ## they are NA, and the result's status is "degenerate"
        if (all(is.na(sp)) || is.na(h1[k]) || is.na(h2[k])) {
            degenerate <- TRUE
            sp <- rep(NA_real_, length(isnps))
            pp <- rep(NA_real_, 5L)
        }
        snpp[[k]] <- sp
        pps[k, ] <- unname(pp)
    }

    ## ---- 8. the result ----------------------------------------------------
    skel <- data.table(nsnps = length(isnps), hit1 = h1, hit2 = h2,
                       PP.H0.abf = pps[, 1], PP.H1.abf = pps[, 2],
                       PP.H2.abf = pps[, 3], PP.H3.abf = pps[, 4],
                       PP.H4.abf = pps[, 5],
                       idx1 = todo$i, idx2 = todo$j)
    PP <- as.data.table(do.call(cbind, snpp))
    colnames(PP) <- if (ncol(PP) > 1) paste0("SNP.PP.H4.row", seq_len(ncol(PP)))
                    else "SNP.PP.H4.abf"

    new_coloc_result(summary = skel,
                     results = cbind(data.table(snp = isnps), PP),
                     prior = prior,
                     weights = if (user_weighted) list(w1 = w1, w2 = w2),
                     status = if (degenerate) "degenerate" else "ok",
                     messages = msgs,
                     meta = list(n_pairs_in = n_in, n_pairs_kept = nrow(todo),
                                 dropped = if (any(drop)) "overlap.min" else character(0)))
}
