## coloc() against 6.0.3. get_match_6_0_3(), in helper-v603.R, runs a call with
## 6.0.3's own functions: get_match_6_0_3(coloc.abf(D1, D2)) is what 6.0.3
## returned for coloc.abf(D1, D2).

D1 <- plain_test_data$D1
D2 <- plain_test_data$D2
D3 <- plain_test_data$D3
D4 <- plain_test_data$D4

fit3 <- suppressMessages(runsusie(D3))
fit4 <- suppressMessages(runsusie(D4))
bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]

## small log BF matrices over 8 snps: two signals in trait 1, one in trait 2
snps <- paste0("s", 1:8)
bf_a <- matrix(c(0, 0, 8, 2, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 7, 1, 0), nrow = 2, byrow = TRUE,
               dimnames = list(c("sig1", "sig2"), snps))
bf_b <- matrix(c(0, 0, 7, 3, 0, 0, 0, 0), nrow = 1, dimnames = list("sig1", snps))
with_null <- function(m) cbind(m, null = 0)

## the same snps listed in the opposite order, as two studies sorted differently would be
reverse <- function(d) {
    for (v in c("beta", "varbeta", "MAF", "snp", "position"))
        d[[v]] <- rev(d[[v]])
    d
}

## coloc.abf()'s SNP.PP.H4, in the snp order of a coloc() result
pp4_of <- function(abf, snps) abf$results$SNP.PP.H4[match(snps, abf$results$snp)]

## Seven pairs of datasets for the two coloc.abf() comparisons below: coloc()
## against 6.0.3's coloc.abf(), and today's coloc.abf() against 6.0.3's.
## The first five cover the kinds of dataset coloc.abf() takes: quant or
## case-control traits, given as beta and varbeta or as p-values, and one
## trait of each type. In the last two, the traits share only some snps (300
## of 500), or list the same snps in a different order, so coloc.abf() must
## match them up by snp name.
D1_sub <- suppressMessages(subset_dataset(D1, 1:400))
D2_sub <- suppressMessages(subset_dataset(D2, 101:500))
abf_pairs <- list(quant = list(D1, D2),                            # beta and varbeta, sdY given
                  cc = list(as_cc(D1), as_cc(D2)),                 # both case-control
                  pvalues = list(to_pvalues(D1), to_pvalues(D2)),  # p-values and MAF only
                  cc_pvalues = list(to_pvalues(as_cc(D1)),         # case-control and p-values
                                    to_pvalues(as_cc(D2))),
                  mixed_type = list(D1, as_cc(D2)),                # one of each type
                  partial_overlap = list(D1_sub, D2_sub),          # 300 of 500 SNPs shared
                  reversed = list(D1, reverse(D2)))                # D2's SNPs in reverse order

## P(no causal variant) for each row of a log BF matrix, with the null at log
## BF 0 and prior p for each snp
p0 <- function(bf, p = 1e-4) {
    bf <- bf[, setdiff(colnames(bf), "null"), drop = FALSE]
    n <- ncol(bf)
    unname((1 - n * p) / ((1 - n * p) + p * rowSums(exp(bf))))
}

## ---- coloc() gives 6.0.3's answers ---------------------------------------------
## Both sides run with their default priors unless a test says otherwise:
## p1 = p2 = 1e-4 and p12 = 5e-6, the same in coloc() and in 6.0.3's
## coloc.bf_bf(), coloc.susie() and coloc.susie_bf(), so each test also checks
## that the defaults agree. coloc.abf()'s p12 default is 1e-5 (AL01), so its
## tests give p12 = 1e-5.

test_that("coloc() matches 6.0.3's coloc.bf_bf(): two signals against one", {
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b))
    got <- coloc(bf_a, bf_b)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
    expect_equal(got$priors, expected$priors)
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): named vectors", {
    v <- bf_a[1, ]
    expected <- get_match_6_0_3(coloc.bf_bf(v, v))
    got <- coloc(v, v)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): a null column", {
    for (pair in list(both = list(with_null(bf_a), with_null(bf_b)),
                      one = list(with_null(bf_a), bf_b),
                      ## nulllast() has to move it
                      not_last = list(with_null(bf_a)[, c("null", snps)], bf_b),
                      ## one row with a null column: a subset of one row
                      ## becomes a vector unless the code says drop = FALSE
                      one_row = list(with_null(bf_b), with_null(bf_b)))) {
        expected <- get_match_6_0_3(coloc.bf_bf(pair[[1]], pair[[2]]))
        got <- coloc(pair[[1]], pair[[2]])
        expect_equal(got$summary, expected$summary)
        expect_equal(got$results, expected$results)
    }
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): partial snp overlap", {
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a[, 1:6], bf_b))
    got <- coloc(bf_a[, 1:6], bf_b)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): trim_by_posterior = FALSE", {
    far <- bf_a[, 5:8]
    expected <- get_match_6_0_3(coloc.bf_bf(far, bf_b, trim_by_posterior = FALSE))
    got <- coloc(far, bf_b, control = coloc_control(trim_by_posterior = FALSE))
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): p12 and prior weights", {
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, p12 = 1e-5))
    got <- coloc(bf_a, bf_b, prior = coloc_prior(p12 = 1e-5))
    expect_equal(got$summary, expected$summary)
    expect_equal(got$priors, expected$priors)

    w <- c(1, 1, 1, 1, 1, 1, 1, 100)
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, prior_weights1 = w, prior_weights2 = w))
    ## coloc_prior() warns when both traits have weights (AL15, below)
    prior <- suppressWarnings(coloc_prior(weights1 = w, weights2 = w))
    got <- coloc(bf_a, bf_b, prior = prior)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.abf()", {
    for (nm in names(abf_pairs)) {
        d1 <- abf_pairs[[nm]][[1]]
        d2 <- abf_pairs[[nm]][[2]]
        ## p12 = 1e-5 is coloc.abf()'s default (AL01)
        expected <- get_match_6_0_3(coloc.abf(d1, d2))
        got <- coloc(d1, d2, prior = coloc_prior(p12 = 1e-5))
        ## as lists, so each posterior meets the tolerance on its own scale
        expect_equal(as.list(pp(got)), as.list(expected$summary), label = nm)
        expect_equal(got$results$SNP.PP.H4.abf, pp4_of(expected, got$results$snp), label = nm)
        expect_equal(got$priors, expected$priors, label = nm)
    }
})

test_that("coloc() matches 6.0.3's coloc.abf(): sdY estimated", {
    d1 <- D1
    d2 <- D2
    d1$sdY <- d2$sdY <- NULL
    ## both estimate sdY, and say so
    expected <- suppressWarnings(get_match_6_0_3(coloc.abf(d1, d2)))
    got <- suppressWarnings(coloc(d1, d2, prior = coloc_prior(p12 = 1e-5)))
    expect_equal(as.list(pp(got)), as.list(expected$summary))
    expect_equal(got$results$SNP.PP.H4.abf, pp4_of(expected, got$results$snp))
})

test_that("coloc() stops on a bad dataset with 6.0.3's coloc.abf() message, naming the trait", {
    bad <- D2
    bad$beta[1] <- NA
    expected <- tryCatch(get_match_6_0_3(coloc.abf(D1, bad)), error = conditionMessage)
    expect_equal(expected, "dataset 2: beta contains missing values")
    expect_error(coloc(D1, bad), "trait 2: beta contains missing values", fixed = TRUE)
})

test_that("SuSiE fits each test dataset the same way every time", {
    ## The SuSiE comparisons fit each dataset twice, once for each side, so they
    ## agree only if SuSiE gives the same fit both times. It has one random step:
    ## for a credible set of more than n_purity (100) snps, it checks a random 100
    ## of them. So every test dataset's credible sets must stay smaller than that.
    threshold <- eval(formals(susieR::susie_get_cs)$n_purity)
    for (nm in c("D1", "D2", "D3", "D4")) {
        fit <- suppressMessages(runsusie(plain_test_data[[nm]]))
        expect_lt(max(lengths(fit$sets$cs)), threshold, label = paste(nm, "largest credible set"))
    }
})

test_that("coloc() matches 6.0.3's coloc.susie(): two coloc datasets", {
    expected <- suppressMessages(get_match_6_0_3(coloc.susie(D3, D4)))
    got <- suppressMessages(coloc(D3, D4, method = "susie"))
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
    expect_equal(got$priors, expected$priors)
})

test_that("coloc() matches 6.0.3's coloc.susie(): two SuSiE fits", {
    ## D4 has one credible set, so as trait 1 it is a one-row matrix: a subset
    ## of one row becomes a vector unless the code says drop = FALSE
    for (pair in list(list(fit4, fit3), list(fit3, fit3))) {
        expected <- get_match_6_0_3(coloc.susie(pair[[1]], pair[[2]]))
        got <- coloc(pair[[1]], pair[[2]])
        expect_equal(got$summary, expected$summary)
        expect_equal(got$results, expected$results)
    }
})

test_that("coloc() matches 6.0.3's coloc.susie_bf()", {
    expected <- get_match_6_0_3(coloc.susie_bf(fit3, bf4))
    got <- coloc(fit3, bf4)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
    expect_equal(got$priors, expected$priors)
})

## ---- today's versions of the changed functions, before and after ---------------
## The old entry points, and process.dataset() and runsusie(), which the data
## stage changes, against their 6.0.3 copies. The tests of coloc.abf(),
## coloc.bf_bf(), coloc.susie() and coloc.susie_bf() go with those functions.

test_that("coloc.abf() returns what 6.0.3's did", {
    for (nm in names(abf_pairs)) {
        d1 <- abf_pairs[[nm]][[1]]
        d2 <- abf_pairs[[nm]][[2]]
        expected <- get_match_6_0_3(coloc.abf(d1, d2))
        invisible(utils::capture.output(got <- coloc.abf(d1, d2)))
        expect_equal(got, expected, label = nm)
        expect_equal(as.list(got$summary), as.list(expected$summary), label = nm)
    }
})

test_that("process.dataset() returns what 6.0.3's did, for every kind of dataset", {
    ## every dataset in abf_pairs; one whose sdY has to be estimated; both
    ## routes without position; one with no small p-value, which 6.0.3's
    ## process.dataset() didn't warn about, as it never checked its dataset;
    ## and named vectors, which gave 6.0.3's result its row names. Then what
    ## those miss: a factor type, as a data frame column gives (with levels cc
    ## and quant, or its own level alone); one snp, given as p-values (with
    ## beta and varbeta it is AL21); no snps; and, unchecked, a type that is
    ## neither quant nor cc, and a dataset with neither route
    D1_no_sdY <- D1
    D1_no_sdY$sdY <- NULL
    weak <- D1
    weak$beta <- weak$beta / 100
    beta_named <- D1
    names(beta_named$beta) <- D1$snp
    both_named <- beta_named
    names(both_named$varbeta) <- paste0("v", D1$snp)
    pvalues_named <- to_pvalues(D1)
    names(pvalues_named$pvalues) <- D1$snp
    factor_type <- D1
    factor_type$type <- factor("quant", levels = c("cc", "quant"))
    factor_cc_pvalues <- to_pvalues(as_cc(D1))
    factor_cc_pvalues$type <- factor("cc")
    one_snp_pvalues <- to_pvalues(D1)
    for (v in c("pvalues", "MAF", "snp", "position"))
        one_snp_pvalues[[v]] <- one_snp_pvalues[[v]][105]
    no_snps <- no_position(D1)
    for (v in c("beta", "varbeta", "MAF", "snp"))
        no_snps[[v]] <- no_snps[[v]][0]
    binary <- D1
    binary$type <- "binary"
    no_route <- to_pvalues(D1)
    no_route$pvalues <- NULL
    datasets <- c(unlist(abf_pairs, recursive = FALSE),
                  list(no_sdY = D1_no_sdY,
                       no_position = no_position(D1),
                       pvalues_no_position = no_position(to_pvalues(D1)),
                       weak = weak,
                       beta_named = beta_named,
                       both_named = both_named,
                       pvalues_named = pvalues_named,
                       factor_type = factor_type,
                       factor_cc_pvalues = factor_cc_pvalues,
                       one_snp_pvalues = one_snp_pvalues,
                       no_snps = no_snps,
                       binary = binary,
                       no_route = no_route))
    ## the value; the warnings, which suppressWarnings() would hide; and an
    ## error's message and call, as 6.0.3's process.dataset() raised its own
    run <- function(expr) {
        warned <- character(0)
        value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
            warned <<- c(warned, conditionMessage(w))
            invokeRestart("muffleWarning")
        }), error = function(e) list(error = conditionMessage(e),
                                     call = deparse(conditionCall(e))))
        list(value = value, warned = warned)
    }
    ## coloc.abf() names the columns df1 and df2, finemap.abf() "", which
    ## leaves a dot; NULL leaves V, z, r and lABF bare, but on the p-value
    ## route gives pvalues., MAF. and N.
    suffixes <- list("df1", "", NULL)
    for (nm in names(datasets)) {
        d <- datasets[[nm]]
        for (suffix in suffixes)
            expect_identical(run(process.dataset(d, suffix = suffix)),
                             run(get_match_6_0_3(process.dataset(d, suffix = suffix))),
                             label = paste0(nm, ", suffix ", deparse(suffix)))
    }
    expect_error(process.dataset(no_route, suffix = "df1"),
                 "Must give, as a minimum, one of", fixed = TRUE)
    ## ... is ignored, as it was (AL06)
    wide <- c(quant = 0.3, cc = 0.4)
    expect_identical(process.dataset(D1, suffix = "df1", effect_priors = wide),
                     get_match_6_0_3(process.dataset(D1, suffix = "df1", effect_priors = wide)))
    expect_identical(process.dataset(D1, suffix = "df1", effect_priors = wide),
                     process.dataset(D1, suffix = "df1"))
})

test_that("runsusie() returns what 6.0.3's did", {
    ## the whole fit, so what annotate_susie() and .susie_setld() add is compared too
    expect_equal(fit3, suppressMessages(get_match_6_0_3(runsusie(D3))))
    expect_equal(fit4, suppressMessages(get_match_6_0_3(runsusie(D4))))
})

test_that("coloc.bf_bf() returns what 6.0.3's did", {
    expect_equal(coloc.bf_bf(bf_a, bf_b), get_match_6_0_3(coloc.bf_bf(bf_a, bf_b)))
})

test_that("coloc.bf_bf() with every pair trimmed returns 6.0.3's result, apart from PP.H0", {
    ## PP.H0 is the one intended change here (AL08 and AL13, below)
    far <- bf_a[, 5:8]
    expect_warning(got <- coloc.bf_bf(far, bf_b), "overlap too small")
    expect_warning(expected <- get_match_6_0_3(coloc.bf_bf(far, bf_b)), "overlap too small")
    expect_equal(names(got), names(expected))
    expect_equal(sapply(got$summary, class), sapply(expected$summary, class))
    cols <- setdiff(names(expected$summary), "PP.H0.abf")
    expect_equal(got$summary[, cols, with = FALSE], expected$summary[, cols, with = FALSE])
})

test_that("coloc.susie() returns what 6.0.3's did", {
    expect_equal(coloc.susie(fit3, fit4), get_match_6_0_3(coloc.susie(fit3, fit4)))
})

test_that("coloc.susie_bf() returns what 6.0.3's did", {
    expect_equal(coloc.susie_bf(fit3, bf4), get_match_6_0_3(coloc.susie_bf(fit3, bf4)))
})

## ---- the deliberate differences, one test each -------------------------------
## Each test shows 6.0.3's answer and the new one. The earlier tests carry the
## id of the difference in the design notes' allowlist; the later ones say what
## changed instead.

test_that("AL01: coloc() defaults p12 to 5e-6, where coloc.abf() used 1e-5", {
    expected <- get_match_6_0_3(coloc.abf(D1, D2))
    expect_equal(expected$priors[["p12"]], 1e-5)
    expect_equal(coloc(D1, D2)$priors[["p12"]], 5e-6)
    expect_equal(as.list(pp(coloc(D1, D2, prior = coloc_prior(p12 = 1e-5)))),
                 as.list(expected$summary))
})

test_that("AL02: the priors given are used, which 6.0.3's coloc.susie_bf() ignored", {
    expected <- get_match_6_0_3(coloc.susie_bf(fit3, bf4, p12 = 1e-4))
    expect_equal(expected$priors[["p12"]], 5e-6)
    got <- coloc(fit3, bf4, prior = coloc_prior(p12 = 1e-4))
    expect_equal(got$priors[["p12"]], 1e-4)
})

test_that("AL03: coloc.susie() now applies the region-size guard, which 6.0.3's never did", {
    ## Today's coloc.susie() runs through coloc(), which applies the region-size
    ## guard; 6.0.3's ran through coloc.bf_bf(), which does not. With 500 SNPs,
    ## p1 = p2 = 0.01 implies five causal variants in each trait, so coloc() caps
    ## both at 1/501. coloc.susie_bf() also runs through coloc() now, but cannot
    ## show this: 6.0.3's ignored the priors it was given (AL02).
    w <- capture_warnings(expected <- get_match_6_0_3(
        coloc.susie(fit3, fit4, p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    expect_length(w, 0)
    expect_equal(expected$priors, c(p1 = 0.01, p2 = 0.01, p12 = 1e-4))
    w <- capture_warnings(got <- coloc.susie(fit3, fit4, p1 = 0.01, p2 = 0.01, p12 = 1e-4))
    expect_equal(w, c("p1 * nsnps >= 1, setting p1=1/(nsnps + 1)",
                      "p2 * nsnps >= 1, setting p2=1/(nsnps + 1)"))
    expect_equal(got$priors, c(p1 = 1 / 501, p2 = 1 / 501, p12 = 1e-4))
    expect_gt(got$summary$PP.H4.abf[1], expected$summary$PP.H4.abf[1])
    ## and nothing else moves: the same as 6.0.3 given the capped priors
    capped <- get_match_6_0_3(coloc.susie(fit3, fit4, p1 = 1 / 501, p2 = 1 / 501, p12 = 1e-4))
    expect_equal(got$summary, capped$summary)
})

test_that("AL08: PP.H0 in a trimmed result is a probability, where 6.0.3 left rounding noise", {
    far <- bf_a[, 5:8]
    expect_warning(expected <- get_match_6_0_3(coloc.bf_bf(far, bf_b)), "overlap too small")
    expect_warning(got <- coloc(far, bf_b), "overlap too small")
    expect_true(all(abs(expected$summary$PP.H0.abf) < 1e-15))
    expect_equal(got$summary$PP.H0.abf, p0(far) * p0(bf_b))
    ## and nothing else moves
    cols <- c("nsnps", "hit1", "hit2", "idx1", "idx2")
    expect_equal(got$summary[, cols, with = FALSE], expected$summary[, cols, with = FALSE])
})

test_that("AL10: a vector p1 or p2 is rejected, where 6.0.3 allowed it, with extra hypotheses", {
    ## p1 and p2 of length 8, one per snp: 6.0.3 accepted them and reported
    ## 4 * 8 + 1 = 33 hypotheses in place of five; the count grows with the
    ## length of the vectors
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, p1 = rep(1e-4, 8), p2 = rep(1e-4, 8)))
    expect_equal(sum(startsWith(names(expected$summary), "PP.H")), 33)
    expect_error(coloc(bf_a, bf_b, prior = coloc_prior(p1 = rep(1e-4, 8))),
                 "must each be a scalar")
    expect_error(coloc.bf_bf(bf_a, bf_b, p1 = rep(1e-4, 8)), "must each be a scalar")
})

test_that("AL11: traits with no snps in common warn, where 6.0.3 returned nsnps = NA silently", {
    ## 6.0.3's coloc.bf_bf() returned early with data.table(nsnps = NA) and no warning;
    ## coloc() warns "no snps in common" and returns status "no_overlap"
    bf_x <- bf_b
    colnames(bf_x) <- paste0("x", 1:8)
    expect_silent(expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_x)))
    expect_true(is.na(expected$nsnps))
    expect_warning(got <- coloc(bf_a, bf_x), "no snps in common")
    expect_equal(got$status, "no_overlap")
})

test_that("AL12: an all-NA signal is rejected, where 6.0.3 reported PP.H2 = 1", {
    bf_na <- bf_a
    bf_na[1, ] <- NA
    expected <- get_match_6_0_3(coloc.bf_bf(bf_na, bf_b, trim_by_posterior = FALSE))
    expect_equal(expected$summary$PP.H2.abf[1], 1)
    expect_equal(expected$summary$hit1[1], "-")
    expect_error(coloc(bf_na, bf_b), "missing or infinite")
})

test_that("AL13: trimmed PP.H0 is the joint of the traits, where 6.0.3 took the smaller", {
    far <- with_null(bf_a[, 5:8])
    b <- with_null(bf_b)
    expect_warning(expected <- get_match_6_0_3(coloc.bf_bf(far, b)), "overlap too small")
    expect_warning(got <- coloc(far, b), "overlap too small")
    expect_equal(expected$summary$PP.H0.abf, pmin(p0(far), p0(b)))
    expect_equal(got$summary$PP.H0.abf, p0(far) * p0(b))
    cols <- c("nsnps", "hit1", "hit2", "idx1", "idx2")
    expect_equal(got$summary[, cols, with = FALSE], expected$summary[, cols, with = FALSE])
})

test_that("AL14: a dataset with a BF matrix warns that method cannot apply to the matrix", {
    expected <- suppressMessages(get_match_6_0_3(coloc.susie_bf(D3, bf4)))
    expect_warning(got <- suppressMessages(coloc(D3, bf4, method = "susie")),
                   "method is inapplicable")
    expect_equal(got$summary, expected$summary)
})

test_that("AL15: coloc_prior() warns when both traits have weights, as coloc.abf() did", {
    w <- c(1, 1, 1, 1, 1, 1, 1, 100)
    expect_silent(expected <- get_match_6_0_3(
        coloc.bf_bf(bf_a, bf_b, prior_weights1 = w, prior_weights2 = w)))
    expect_warning(prior <- coloc_prior(weights1 = w, weights2 = w), "both traits")
    expect_equal(coloc(bf_a, bf_b, prior = prior)$summary, expected$summary)
})

test_that("AL16: a p-value above 1 is refused, which 6.0.3 accepted", {
    msg <- "dataset 1: pvalues must be above 0 and at most 1"
    above <- to_pvalues(D1)
    above$pvalues[5] <- 1.5
    expect_null(get_match_6_0_3(check_dataset(above, 1)))
    expect_error(check_dataset(above, 1), msg, fixed = TRUE)
    expect_error(coloc(above, D2), "trait 1: pvalues must be above 0 and at most 1", fixed = TRUE)
    expect_error(finemap.abf(above), "pvalues must be above 0 and at most 1", fixed = TRUE)
    ## 0 and below were refused already; only the message moved
    zero <- to_pvalues(D1)
    zero$pvalues[5] <- 0
    expect_error(get_match_6_0_3(check_dataset(zero, 1)),
                 "pvalues should not be negative or exactly 0")
    expect_error(check_dataset(zero, 1), msg, fixed = TRUE)
    ## exactly 1 is a p-value
    one <- to_pvalues(D1)
    one$pvalues[5] <- 1
    expect_null(check_dataset(one, 1))
})

test_that("AL17: check_dataset() reports a dataset with no type, where 6.0.3 gave R's error", {
    ## runsusie() leaves type out of req, so check_dataset() has to notice a
    ## dataset without one by itself
    d <- D3
    d$type <- NULL
    req <- c("beta", "varbeta", "LD", "snp", "N")   # runsusie()'s, without type
    expect_error(get_match_6_0_3(check_dataset(d, 1, req = req)), "argument is of length zero")
    expect_error(check_dataset(d, 1, req = req), "dataset 1: variable type not set")
    ## with type required, as coloc.abf() asks, the two agree
    expect_error(get_match_6_0_3(check_dataset(d, 1)), "missing required element")
    expect_error(check_dataset(d, 1), "missing required element")
})

test_that("AL18: a SuSiE fit with a duplicated snp id is refused, where 6.0.3 used the first", {
    ## today's coloc.susie() runs through coloc(), which checks a fit's log
    ## BFs as it checks log BFs given directly; 6.0.3's went to coloc.bf_bf()
    dup <- fit4
    colnames(dup$lbf_variable)[2] <- colnames(dup$lbf_variable)[1]
    expected <- get_match_6_0_3(coloc.susie(fit3, dup))
    expect_equal(nrow(expected$summary), 2L)
    expect_error(coloc.susie(fit3, dup),
                 "trait 2: duplicated snps found in the log Bayes factors")
    ## an NA in a credible set's log BFs stopped 6.0.3 too, with R's own message
    na <- fit4
    na$lbf_variable[na$sets$cs_index[1], 10] <- NA
    expect_error(get_match_6_0_3(coloc.susie(fit3, na)), "missing value where TRUE/FALSE needed")
    expect_error(coloc.susie(fit3, na),
                 "trait 2: log Bayes factors contain missing or infinite values")
})

test_that("AL21: process.dataset() no longer names a one-snp result's row after its type", {
    ## 6.0.3's approx.bf.estimates() took the prior as effect_priors["quant"],
    ## which keeps the name, so with one snp, and no names on beta or varbeta,
    ## data.frame() named the row "quant" or "cc". Its p-value route used an
    ## unnamed prior and never did
    one <- D1
    for (v in c("beta", "varbeta", "MAF", "snp", "position"))
        one[[v]] <- one[[v]][105]
    for (d in list(one, as_cc(one))) {
        got <- process.dataset(d, suffix = "df1")
        expected <- get_match_6_0_3(process.dataset(d, suffix = "df1"))
        expect_identical(rownames(expected), d$type)
        expect_identical(rownames(got), "1")
        ## nothing else differs
        rownames(expected) <- NULL
        expect_identical(got, expected)
    }
})

test_that("check_dataset() has no prefix without a label, and \"dataset N: \" on every message with one, unlike 6.0.3", {
    ## 6.0.3 began its messages "dataset N: ", or "dataset : " without a
    ## label; four of them named no dataset, and nor did the minimum p-value
    ## warning. The expected text is written out, so a mistake shared by
    ## reprefix() and .prefix() can't pass
    msg <- function(expr) tryCatch(expr, error = conditionMessage)
    bad <- D2
    bad$beta[1] <- NA
    ## with a label, as 6.0.3 had it
    expect_identical(msg(get_match_6_0_3(coloc.abf(D1, bad))),
                     "dataset 2: beta contains missing values")
    expect_identical(msg(coloc.abf(D1, bad)), "dataset 2: beta contains missing values")
    expect_identical(msg(check_dataset(bad, 2)), "dataset 2: beta contains missing values")
    bad3 <- D3
    bad3$beta[1] <- NA
    expect_identical(msg(get_match_6_0_3(runsusie(bad3))),
                     "dataset 1: beta contains missing values")
    expect_identical(msg(runsusie(bad3)), "dataset 1: beta contains missing values")
    ## without a label (helper-v603.R has no copy of finemap.abf(), so 6.0.3's
    ## check_dataset() stands for it)
    expect_identical(msg(get_match_6_0_3(check_dataset(bad))),
                     "dataset : beta contains missing values")
    expect_identical(msg(finemap.abf(bad)), "beta contains missing values")
    expect_identical(msg(check_dataset(bad)), "beta contains missing values")
    ## without a label, the one message with no subject gains one
    expect_identical(msg(get_match_6_0_3(check_dataset(1:3))), "dataset : is not a list")
    expect_identical(msg(check_dataset(1:3)), "the dataset is not a list")
    expect_identical(msg(finemap.abf(1:3)), "the dataset is not a list")
    expect_identical(msg(check_dataset(1:3, 2)), "dataset 2: is not a list")
    ## the messages that named no dataset (the N message is reworded too, below)
    per_snp_N <- D2
    per_snp_N$N <- rep(D2$N, length(D2$snp))
    expect_identical(msg(get_match_6_0_3(coloc.abf(D1, per_snp_N))),
                     "N should be a single numeric. coloc does not accomodate per-snp sample sizes.")
    expect_identical(msg(coloc.abf(D1, per_snp_N)),
                     "dataset 2: N must be a single number; coloc does not take per-snp sample sizes")
    not_square <- D2
    not_square$LD <- not_square$LD[, -1]
    expect_identical(msg(get_match_6_0_3(coloc.abf(D1, not_square))), "LD not square")
    expect_identical(msg(coloc.abf(D1, not_square)), "dataset 2: LD not square")
    ## the minimum p-value warning, its first line
    weak <- D1
    weak$beta <- weak$beta / 100
    weak2 <- D2
    weak2$beta <- weak2$beta / 100
    warned <- function(expr) sub("\n.*", "", outcome(expr)$warned)
    expect_identical(warned(get_match_6_0_3(coloc.abf(weak, weak2))),
                     c("minimum p value is: 0.93383", "minimum p value is: 0.93863"))
    expect_identical(warned(coloc.abf(weak, weak2)),
                     c("dataset 1: minimum p value is: 0.93383",
                       "dataset 2: minimum p value is: 0.93863"))
    expect_identical(warned(check_dataset(weak, 2)), "dataset 2: minimum p value is: 0.93383")
    expect_identical(warned(check_dataset(weak)), "minimum p value is: 0.93383")
})

test_that("check_dataset() doesn't repeat the minimum p-value warning for a coloc_dataset, where 6.0.3 did", {
    ## PR #2: a coloc_dataset gave the warning when it was made; a plain list
    ## still gives it each time, as test-dataset.R checks
    weak <- D1
    weak$beta <- weak$beta / 100
    made <- suppressWarnings(as_coloc_dataset(weak))
    now <- outcome(check_dataset(made, 1))
    then <- outcome(get_match_6_0_3(check_dataset(made, 1)))
    expect_length(now$warned, 0)
    expect_match(then$warned, "minimum p value is: ", fixed = TRUE)
    ## the warning is the only difference
    expect_identical(now[c("value", "error")], then[c("value", "error")])
    ## a coloc_dataset changed to hold a negative varbeta is refused before the
    ## p-values are worked out, where 6.0.3 went on to compare a NaN and stopped
    made$varbeta[1] <- -made$varbeta[1]
    now <- outcome(check_dataset(made, 1))
    then <- outcome(get_match_6_0_3(check_dataset(made, 1)))
    expect_identical(now$error, "dataset 1: varbeta must be positive")
    expect_length(now$warned, 0)
    expect_identical(then$error, "missing value where TRUE/FALSE needed")
    expect_identical(then$warned, "NaNs produced")
})

test_that("coloc.susie_bf() refuses a bf2 that isn't log Bayes factors, where 6.0.3 failed with R's errors", {
    ## coloc() sorts its inputs by shape, so it would fine-map a dataset or a
    ## fit given as bf2. 6.0.3 failed on both inside coloc.bf_bf()
    msg <- "bf2 must be log Bayes factors: a named numeric vector or a matrix"
    expect_error(suppressMessages(get_match_6_0_3(coloc.susie_bf(fit3, D4))))
    expect_error(suppressMessages(get_match_6_0_3(coloc.susie_bf(fit3, fit4))))
    expect_error(coloc.susie_bf(fit3, D4), msg, fixed = TRUE)
    expect_error(coloc.susie_bf(fit3, fit4), msg, fixed = TRUE)
    chr <- bf4
    storage.mode(chr) <- "character"
    expect_error(coloc.susie_bf(fit3, chr), msg, fixed = TRUE)
})

## check_alignment()'s share counted directly: the pairs of snps above the
## diagonal of the snps' own LD with |LD| above thr, and of those, the share
## whose product of z-scores has the sign of their LD
alignment_share <- function(d, thr = 0.2) {
    z <- d$beta / sqrt(d$varbeta)
    ld <- d$LD[d$snp, d$snp]
    pair <- upper.tri(ld) & abs(ld) > thr
    mean((outer(z, z) * ld)[pair] > 0)
}

test_that("check_alignment() counts each pair of snps once and uses thr, where 6.0.3 counted every snp with itself and took 0.2", {
    for (d in list(D1, D2, D3, D4)) {
        now <- check_alignment(d, do_plot = FALSE)
        expect_equal(now, alignment_share(d))
        ## 6.0.3 also counted each snp paired with itself, which always agrees
        expect_equal(now, get_match_6_0_3(check_alignment(d, do_plot = FALSE)), tolerance = 1e-3)
    }
    expect_equal(check_alignment(D1, thr = 0.5, do_plot = FALSE), alignment_share(D1, 0.5))
    expect_identical(get_match_6_0_3(check_alignment(D1, thr = 0.5, do_plot = FALSE)),
                     get_match_6_0_3(check_alignment(D1, do_plot = FALSE)))
    ## the plot is drawn from the same pairs
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    expect_equal(check_alignment(D1), alignment_share(D1))
})

test_that("check_alignment() matches LD to the snps by name, where 6.0.3 divided by position", {
    reversed <- D1
    reversed$LD <- D1$LD[500:1, 500:1]
    wider <- D1
    wider$LD <- rbind(cbind(D1$LD, 0), 0)
    wider$LD[501, 501] <- 1
    dimnames(wider$LD) <- rep(list(c(D1$snp, "extra")), 2)
    expect_equal(check_alignment(reversed, do_plot = FALSE), alignment_share(D1))
    expect_equal(check_alignment(wider, do_plot = FALSE), alignment_share(D1))
    expect_false(isTRUE(all.equal(get_match_6_0_3(check_alignment(reversed, do_plot = FALSE)),
                                  alignment_share(D1))))
    expect_error(get_match_6_0_3(check_alignment(wider, do_plot = FALSE)), "non-conformable arrays")
})

test_that("check_alignment() counts the pairs across blocks of a large LD", {
    ## D1's snps twice over, each copy interleaved with the other, so that the
    ## pairs within a copy fall across the 500-column blocks
    twice <- function(v) rep(v, each = 2)
    big <- list(beta = twice(D1$beta), varbeta = twice(D1$varbeta), type = D1$type,
                sdY = D1$sdY, snp = as.vector(rbind(paste0(D1$snp, "a"), paste0(D1$snp, "b"))),
                LD = kronecker(D1$LD, diag(2)))
    dimnames(big$LD) <- list(big$snp, big$snp)
    expect_equal(check_alignment(big, do_plot = FALSE), alignment_share(big))
    expect_equal(check_alignment(big, do_plot = FALSE), alignment_share(D1))
})

## ---- values coloc can't use, at each way a dataset is checked ---------------
## check_dataset() checks a dataset for coloc.abf(), finemap.abf(), the
## constructors, wakefield_bf(), runsusie() and coloc(). at_each_door() gives
## the error each of them stops with, NA where it goes on, and everywhere() the
## errors they should give: one sentence, with each one's prefix. runsusie() is
## asked only about a dataset it could fine-map. On 6.0.3's side an R error is
## tested with no text, as R's wording changes between its versions.

at_each_door <- function(d, other = D1, susie = FALSE) {
    msg <- function(expr) tryCatch({
        utils::capture.output(suppressWarnings(suppressMessages(expr)))
        NA_character_
    }, error = conditionMessage)
    doors <- c(check_dataset = msg(check_dataset(d, 2)),
               coloc.abf = msg(coloc.abf(other, d)),
               finemap.abf = msg(finemap.abf(d)),
               as_coloc_dataset = msg(as_coloc_dataset(d)),
               coloc = msg(coloc(other, d)),
               wakefield_bf = msg(wakefield_bf(d)))
    if (susie)
        doors <- c(doors, runsusie = msg(runsusie(d)))
    doors
}

everywhere <- function(sentence, susie = FALSE) {
    prefix <- c(check_dataset = "dataset 2: ", coloc.abf = "dataset 2: ", finemap.abf = "",
                as_coloc_dataset = "", coloc = "trait 2: ", wakefield_bf = "")
    if (susie)
        prefix <- c(prefix, runsusie = "dataset 1: ")
    prefix[] <- if (is.na(sentence)) NA_character_ else paste0(prefix, sentence)
    prefix
}

## coloc.abf()'s summary, without its printing
abf_summary <- function(d1, d2) {
    utils::capture.output(res <- coloc.abf(d1, d2))
    res$summary
}

## the error a call stops with, NA if none, or "SuSiE ran first" if SuSiE
## said anything before it, as it does when it starts fine-mapping
error_before_susie <- function(expr) {
    said <- FALSE
    err <- withCallingHandlers(
        tryCatch({
            expr
            NA_character_
        }, error = conditionMessage),
        message = function(m) {
            said <<- TRUE
            invokeRestart("muffleMessage")
        })
    if (said) "SuSiE ran first" else err
}

test_that("beta, varbeta and p-values must be numeric where they are read, where 6.0.3 failed with R's errors", {
    chr_beta <- D1
    chr_beta$beta <- as.character(chr_beta$beta)
    chr_varbeta <- D1
    chr_varbeta$varbeta <- as.character(chr_varbeta$varbeta)
    chr_p <- to_pvalues(D1)
    chr_p$pvalues <- as.character(chr_p$pvalues)
    expect_error(get_match_6_0_3(check_dataset(chr_beta, 2)))
    expect_error(get_match_6_0_3(check_dataset(chr_varbeta, 2)))
    expect_identical(at_each_door(chr_beta, susie = TRUE),
                     everywhere("beta must be numeric", susie = TRUE))
    expect_identical(at_each_door(chr_varbeta, susie = TRUE),
                     everywhere("varbeta must be numeric", susie = TRUE))
    ## 6.0.3's check let character p-values through, and coloc.abf() failed after
    expect_null(get_match_6_0_3(check_dataset(chr_p, 2)))
    expect_error(get_match_6_0_3(coloc.abf(D1, chr_p)))
    expect_identical(at_each_door(chr_p), everywhere("pvalues must be numeric"))
    ## 6.0.3 computed with a logical beta
    lgl_beta <- D1
    lgl_beta$beta <- lgl_beta$beta > 0
    expect_null(get_match_6_0_3(check_dataset(lgl_beta, 2)))
    expect_identical(at_each_door(lgl_beta), everywhere("beta must be numeric"))
    ## where nothing reads them they pass, with 6.0.3's answer
    p_beside_beta <- D1
    p_beside_beta$pvalues <- chr_p$pvalues
    beta_beside_p <- to_pvalues(D1)
    beta_beside_p$beta <- chr_beta$beta
    varbeta_beside_p <- to_pvalues(D1)
    varbeta_beside_p$varbeta <- chr_varbeta$varbeta
    for (d in list(p_beside_beta, beta_beside_p, varbeta_beside_p)) {
        expect_identical(at_each_door(d), everywhere(NA))
        expect_equal(abf_summary(D1, d), get_match_6_0_3(coloc.abf(D1, d))$summary)
    }
    ## with a zero varbeta as well: before, check_dataset() reached the zero first
    zero_too <- chr_beta
    zero_too$varbeta[1] <- 0
    expect_error(get_match_6_0_3(check_dataset(zero_too, 2)))
    expect_identical(at_each_door(zero_too), everywhere("beta must be numeric"))
})

test_that("a negative varbeta is refused, where 6.0.3 failed with R's error", {
    neg <- D1
    neg$varbeta[5] <- -neg$varbeta[5]
    expect_error(suppressWarnings(get_match_6_0_3(check_dataset(neg, 2))))
    expect_identical(at_each_door(neg, susie = TRUE),
                     everywhere("varbeta must be positive", susie = TRUE))
    ## with a zero as well, 6.0.3's message for the zero
    zero_too <- neg
    zero_too$varbeta[1] <- 0
    expect_identical(outcome(get_match_6_0_3(check_dataset(zero_too, 2)))$error,
                     "dataset 2: Zero values in varbeta")
    expect_identical(at_each_door(zero_too), everywhere("Zero values in varbeta"))
    ## beside p-values varbeta isn't read, and passes as in 6.0.3
    beside_p <- to_pvalues(D1)
    beside_p$varbeta <- -D1$varbeta
    expect_null(get_match_6_0_3(check_dataset(beside_p, 2)))
    expect_identical(at_each_door(beside_p), everywhere(NA))
    ## a coloc_dataset changed to hold one gave NA posteriors, as its minimum
    ## p-value comparison is skipped
    made <- as_coloc_dataset(D1)
    made$varbeta[5] <- -made$varbeta[5]
    expect_error(coloc(made, D2), "^trait 1: varbeta must be positive$")
    expect_error(coloc.abf(made, D2), "^dataset 1: varbeta must be positive$")
    expect_error(wakefield_bf(made), "^varbeta must be positive$")
})

test_that("type and s must be single values, where 6.0.3 failed with R's errors", {
    for (type in list(c("quant", "cc"), rep("quant", length(D1$snp)), character(0))) {
        d <- D1
        d$type <- type
        expect_error(get_match_6_0_3(check_dataset(d, 2)))
        expect_identical(at_each_door(d), everywhere("type must be a single value, quant or cc"))
    }
    two_s <- as_cc(D1)
    two_s$s <- c(0.5, 0.5)
    expect_error(get_match_6_0_3(check_dataset(two_s, 2)))
    expect_identical(at_each_door(two_s), everywhere("s must be a single number between 0 and 1"))
})

test_that("LD must be a matrix, and for SuSiE a numeric one, where 6.0.3 failed with R's errors", {
    vector_ld <- D1
    vector_ld$LD <- as.vector(D1$LD)
    expect_error(get_match_6_0_3(check_dataset(vector_ld, 2)))
    expect_identical(at_each_door(vector_ld, susie = TRUE),
                     everywhere("LD must be a matrix", susie = TRUE))
    ## 6.0.3 took an LD with a third dimension
    three_d <- D1
    three_d$LD <- array(D1$LD, dim = c(dim(D1$LD), 1), dimnames = c(dimnames(D1$LD), list(NULL)))
    expect_null(get_match_6_0_3(check_dataset(three_d, 2)))
    expect_identical(at_each_door(three_d, susie = TRUE),
                     everywhere("LD must be a matrix", susie = TRUE))
    ## a data frame LD has two dimensions: where LD isn't required it passes,
    ## with 6.0.3's answer
    df_ld <- D1
    df_ld$LD <- as.data.frame(D1$LD)
    expect_identical(at_each_door(df_ld), everywhere(NA))
    expect_equal(abf_summary(D1, df_ld), get_match_6_0_3(coloc.abf(D1, df_ld))$summary)
    ## SuSiE needs a numeric matrix. 6.0.3's check_ld() refused a data frame
    ## after check_dataset(), and a character matrix failed in R; now the check
    ## refuses both, before either trait is fine-mapped
    df_ld3 <- D3
    df_ld3$LD <- as.data.frame(D3$LD)
    chr_ld3 <- D3
    storage.mode(chr_ld3$LD) <- "character"
    expect_error(get_match_6_0_3(runsusie(df_ld3)), "LD must be of class matrix", fixed = TRUE)
    expect_error(suppressMessages(get_match_6_0_3(runsusie(chr_ld3))))
    for (d in list(df_ld3, chr_ld3)) {
        expect_error(runsusie(d), "^dataset 1: LD must be a numeric matrix$")
        expect_identical(error_before_susie(coloc(D4, d, method = "susie")),
                         "trait 2: LD must be a numeric matrix")
    }
})

test_that("N must be a finite number above 0 wherever it is given, where 6.0.3 took any", {
    ## beside beta and sdY nothing reads it, and 6.0.3 gave its answer
    for (n in c(0, -5, Inf)) {
        d <- D1
        d$N <- n
        expect_null(get_match_6_0_3(check_dataset(d, 2)))
        expect_identical(at_each_door(d, susie = TRUE),
                         everywhere("N must be a finite number above 0", susie = TRUE))
    }
    ## without sdY, which is estimated from N, 6.0.3 went on all the same
    no_sdY <- D1
    no_sdY$sdY <- NULL
    no_sdY$N <- 0
    expect_no_error(suppressWarnings(get_match_6_0_3(coloc.abf(D1, no_sdY))))
    expect_identical(at_each_door(no_sdY), everywhere("N must be a finite number above 0"))
    ## given as p-values, 0 keeps 6.0.3's message, and Inf, which 6.0.3 took,
    ## is refused
    p_zero <- to_pvalues(D1)
    p_zero$N <- 0
    expect_identical(outcome(get_match_6_0_3(check_dataset(p_zero, 2)))$error,
                     "dataset 2: sample size N <=0 or not set")
    expect_identical(at_each_door(p_zero), everywhere("sample size N <=0 or not set"))
    p_inf <- to_pvalues(D1)
    p_inf$N <- Inf
    expect_null(get_match_6_0_3(check_dataset(p_inf, 2)))
    expect_identical(at_each_door(p_inf), everywhere("N must be a finite number above 0"))
    ## for SuSiE, refused before either trait is fine-mapped; 6.0.3 left it to
    ## susieR, after the other trait was
    zero4 <- D4
    zero4$N <- 0
    expect_error(suppressMessages(get_match_6_0_3(runsusie(zero4))))
    expect_identical(error_before_susie(coloc(D3, zero4, method = "susie")),
                     "trait 2: N must be a finite number above 0")
})

test_that("sdY must be a single finite number above 0 where it is read, where 6.0.3 took any", {
    ## a quantitative trait given as beta and varbeta
    for (sdY in list(0, -1, c(1, 1), Inf, "1")) {
        d <- D1
        d$sdY <- sdY
        expect_null(get_match_6_0_3(check_dataset(d, 2)))
        expect_identical(at_each_door(d), everywhere("sdY must be a single finite number above 0"))
    }
    ## where nothing reads it, it passes with 6.0.3's answer: on a case-control
    ## trait, and with p-values
    cc <- as_cc(D1)
    cc$sdY <- 0
    p <- to_pvalues(D1)
    p$sdY <- 0
    for (d in list(cc, p)) {
        expect_identical(at_each_door(d), everywhere(NA))
        expect_equal(abf_summary(D1, d), get_match_6_0_3(coloc.abf(D1, d))$summary)
    }
    ## and for SuSiE, which never reads it
    susie3 <- D3
    susie3$sdY <- 0
    expect_null(check_dataset(susie3, 1, req = .dataset_uses$susie))
})

test_that("for SuSiE LD must be a correlation matrix, where 6.0.3 took any numeric matrix", {
    msg <- "LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal"
    susie_req <- c("beta", "varbeta", "LD", "snp", "N")   # runsusie()'s, in 6.0.3 too
    above_one <- D3
    above_one$LD[1, 2] <- above_one$LD[2, 1] <- 1.5
    covariance <- D3
    covariance$LD <- D3$LD * 0.5
    asymmetric <- D3
    asymmetric$LD[1, 2] <- asymmetric$LD[1, 2] + 0.1
    for (d in list(above_one, covariance, asymmetric)) {
        ## 6.0.3's runsusie() checked the dataset and its LD, and fine-mapped it
        expect_null(get_match_6_0_3(check_dataset(d, 1, req = susie_req)))
        expect_null(get_match_6_0_3(check_ld(d, d$LD)))
        expect_identical(error_before_susie(runsusie(d)), paste0("dataset 1: ", msg))
        expect_identical(error_before_susie(coloc(D4, d, method = "susie")), paste0("trait 2: ", msg))
        ## where LD isn't required it isn't tested, and the answer is 6.0.3's
        expect_null(check_dataset(d, 2))
        expect_equal(abf_summary(D4, d), get_match_6_0_3(coloc.abf(D4, d))$summary)
    }
    ## to within 1e-6, on the diagonal and off it
    susie_check <- function(d) outcome(check_dataset(d, 1, req = .dataset_uses$susie))$error
    for (off in c(5e-7, 2e-6)) {
        diagonal <- D3
        diag(diagonal$LD) <- 1 + off
        pair <- D3
        pair$LD[1, 2] <- pair$LD[1, 2] + off
        expected <- if (off < 1e-6) NA_character_ else paste0("dataset 1: ", msg)
        expect_identical(susie_check(diagonal), expected)
        expect_identical(susie_check(pair), expected)
    }
    ## only the dataset's snps, the part runsusie() uses: an LD covering more
    ## snps passes with a fault among the others, and not with one among them
    fewer <- D3
    for (v in c("beta", "varbeta", "snp", "position", "MAF"))
        fewer[[v]] <- fewer[[v]][-(1:10)]
    outside <- fewer
    outside$LD[1, 2] <- outside$LD[2, 1] <- 1.5
    inside <- fewer
    inside$LD[11, 12] <- inside$LD[12, 11] <- 1.5
    expect_identical(susie_check(outside), NA_character_)
    expect_identical(susie_check(inside), paste0("dataset 1: ", msg))
    ## a non-square LD keeps 6.0.3's message for SuSiE too
    not_square <- D3
    not_square$LD <- D3$LD[, -1]
    expect_identical(susie_check(not_square), "dataset 1: LD not square")
})

test_that(".is_correlation() compares every block of a large LD", {
    set.seed(1)
    ld <- cor(matrix(rnorm(300 * 1200), 300))
    expect_true(.is_correlation(ld))
    ld[1100, 900] <- ld[1100, 900] + 0.1
    expect_false(.is_correlation(ld))
})

test_that("for SuSiE z must have a finite number for each snp, where 6.0.3 failed with R's errors or fit it as given", {
    z3 <- D3$beta / sqrt(D3$varbeta)
    ## each z, whether 6.0.3's runsusie() failed on it with R's or susieR's
    ## error or fine-mapped it as given, and the message it gets now
    broken <- list(list(z = replace(z3, 5, NA), failed = FALSE, msg = "z contains missing values"),
                   list(z = z3[-1], failed = TRUE, msg = "snp has 500 values but z has 499"),
                   list(z = c(z3, 1), failed = TRUE, msg = "snp has 500 values but z has 501"),
                   list(z = as.character(z3), failed = TRUE, msg = "z must be numeric"),
                   list(z = z3 > 0, failed = FALSE, msg = "z must be numeric"),
                   list(z = replace(z3, 5, Inf), failed = FALSE, msg = "Infinite values in z"))
    for (b in broken) {
        d <- D3
        d$z <- b$z
        then <- outcome(suppressMessages(get_match_6_0_3(runsusie(d))))$error
        expect_identical(!is.na(then), b$failed, label = b$msg)
        expect_error(runsusie(d), paste0("^dataset 1: ", b$msg, "$"))
        expect_identical(error_before_susie(coloc(D4, d, method = "susie")),
                         paste0("trait 2: ", b$msg))
        ## where nothing reads z it passes, with 6.0.3's answer
        d1 <- D1
        d1$z <- b$z
        expect_identical(at_each_door(d1), everywhere(NA))
    }
    expect_equal(abf_summary(D1, d1), get_match_6_0_3(coloc.abf(D1, d1))$summary)
})

test_that("runsusie() fine-maps from the dataset's z, as 6.0.3 did", {
    half <- D3
    half$z <- D3$beta / sqrt(D3$varbeta) / 2
    fit_half <- suppressMessages(runsusie(half))
    expect_false(isTRUE(all.equal(fit_half$pip, fit3$pip)))
    expect_identical(fit_half$pip, suppressMessages(get_match_6_0_3(runsusie(half)))$pip)
})

## ---- messages reworded since 6.0.3 ------------------------------------------
## The same datasets are refused, at the same rule, with a sentence that says
## more. at_each_door() and everywhere() are above.

test_that("lengths that differ name the elements, where 6.0.3 said only that they don't match", {
    short <- D1
    short$beta <- short$beta[-1]
    expect_identical(outcome(get_match_6_0_3(check_dataset(short, 2)))$error,
                     "dataset 2: lengths of inputs don't match: ")
    expect_identical(at_each_door(short, susie = TRUE),
                     everywhere("snp has 500 values but beta has 499", susie = TRUE))
    msg <- function(d, ...) outcome(check_dataset(d, ...))$error
    ## the elements are compared in a fixed order, whatever the list's
    two_short <- short
    two_short$MAF <- two_short$MAF[-(1:2)]
    expect_identical(msg(two_short), "snp has 500 values but beta has 499")
    expect_identical(msg(rev(two_short)), "snp has 500 values but beta has 499")
    ## each per-snp element is compared with snp, position too
    short_position <- D1
    short_position$position <- short_position$position[-1]
    expect_identical(msg(short_position), "snp has 500 values but position has 499")
    short_snp <- D1
    short_snp$snp <- short_snp$snp[-1]
    expect_identical(msg(short_snp), "snp has 499 values but beta has 500")
    ## without snp, against the first element given
    no_snp <- short
    no_snp$snp <- NULL
    expect_identical(msg(no_snp, req = "type"), "beta has 499 values but varbeta has 500")
    one <- list(snp = "s1", beta = c(0.1, 0.2), varbeta = c(0.01, 0.01), type = "quant", sdY = 1)
    expect_identical(msg(one), "snp has 1 value but beta has 2")
})

test_that("a per-snp N gets a reworded message, where 6.0.3's misspelt a word", {
    for (n in list(rep(D1$N, length(D1$snp)), as.character(D1$N))) {
        d <- D1
        d$N <- n
        expect_identical(outcome(get_match_6_0_3(check_dataset(d, 2)))$error,
                         "N should be a single numeric. coloc does not accomodate per-snp sample sizes.")
        expect_identical(at_each_door(d, susie = TRUE),
                         everywhere("N must be a single number; coloc does not take per-snp sample sizes",
                                    susie = TRUE))
    }
})

## ---- the dataset functions, pinned as numbers ------------------------------
## The tests above compare process.dataset() and runsusie() with their 6.0.3
## copies. These pin their values as plain numbers too, worked out with 6.0.3's
## code. Both sides of every comparison use the installed susieR, so a change
## to susieR shows only here.

test_that("process.dataset() and sdY.est() still give 6.0.3's numbers", {
    at_s105 <- function(d, suffix) {
        pd <- process.dataset(d, suffix = suffix)
        pd[[paste0("lABF.", suffix)]][pd$snp == "s105"]
    }
    ## each value is written to 12 digits, so compared to 10
    expect_equal(at_s105(D1, "df1"), 28.343423232, tolerance = 1e-10)
    expect_equal(at_s105(to_pvalues(D1), "df1"), 29.9818369694, tolerance = 1e-10)
    expect_equal(at_s105(as_cc(D2), "df2"), 25.698511433, tolerance = 1e-10)
    expect_equal(suppressWarnings(sdY.est(D1$varbeta, D1$MAF, D1$N)), 1.28815854798,
                 tolerance = 1e-10)
})

test_that("runsusie() still finds D3's credible sets", {
    rows <- fit3$lbf_variable[fit3$sets$cs_index, ]
    expect_equal(fit3$sets$cs_index, c(1, 2))
    expect_equal(unname(lengths(fit3$sets$cs)), c(3L, 13L))
    expect_equal(colnames(rows)[apply(rows, 1, which.max)], c("s105", "s89"))
    ## a SuSiE fit can differ across machines by about 1e-8, so 6 digits
    expect_equal(apply(rows, 1, max), c(25.82456833, 13.32711173), tolerance = 1e-6)
})
