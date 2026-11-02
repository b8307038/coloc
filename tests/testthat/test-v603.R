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

## ---- coloc() gives 6.0.3's answers ---------------------------------------------
## Both sides run with their default priors unless a test says otherwise:
## p1 = p2 = 1e-4 and p12 = 5e-6, the same in coloc() and in 6.0.3's
## coloc.bf_bf(), coloc.susie() and coloc.susie_bf(), so each test also checks
## that the defaults agree. coloc.abf()'s p12 default is 1e-5, so its tests
## give p12 = 1e-5.

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
    ## coloc_prior() warns when both traits have weights
    prior <- suppressWarnings(coloc_prior(weights1 = w, weights2 = w))
    got <- coloc(bf_a, bf_b, prior = prior)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.abf()", {
    for (nm in names(abf_pairs)) {
        d1 <- abf_pairs[[nm]][[1]]
        d2 <- abf_pairs[[nm]][[2]]
        ## p12 = 1e-5 is coloc.abf()'s default
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
    d1$sdY <- NULL
    d2$sdY <- NULL
    ## both estimate sdY, and say so
    expected <- suppressWarnings(get_match_6_0_3(coloc.abf(d1, d2)))
    got <- suppressWarnings(coloc(d1, d2, prior = coloc_prior(p12 = 1e-5)))
    expect_equal(as.list(pp(got)), as.list(expected$summary))
    expect_equal(got$results$SNP.PP.H4.abf, pp4_of(expected, got$results$snp))
})

test_that("coloc() stops on a bad dataset, as 6.0.3's coloc.abf() did, naming the trait", {
    bad <- D2
    bad$beta[1] <- NA
    expect_error(get_match_6_0_3(coloc.abf(D1, bad)))
    expect_error(coloc(D1, bad), "dataset 2: beta contains missing values", fixed = TRUE)
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

## ---- the old functions against their 6.0.3 copies ------------------------------
## coloc.abf(), coloc.bf_bf(), coloc.susie(), coloc.susie_bf(),
## process.dataset() and runsusie(), each against its 6.0.3 copy. Each test
## goes with the function it tests.

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
    ## and quant, or its own level alone); one snp, given as p-values, as a test
    ## below covers one snp with beta and varbeta; no snps; and, unchecked, a
    ## type that is neither quant nor cc, and a dataset with neither route
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
    ## the value, or the call that stopped it, as 6.0.3's process.dataset()
    ## raised its own error; and how many warnings it gave, which
    ## suppressWarnings() would hide
    result_of <- function(expr) {
        warned <- 0
        value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
            warned <<- warned + 1
            invokeRestart("muffleWarning")
        }), error = function(e) list(stopped_in = deparse(conditionCall(e))))
        list(value = value, warned = warned)
    }
    ## coloc.abf() names the columns df1 and df2, finemap.abf() "", which
    ## leaves a dot; NULL leaves V, z, r and lABF bare, but on the p-value
    ## route gives pvalues., MAF. and N.
    suffixes <- list("df1", "", NULL)
    for (nm in names(datasets)) {
        d <- datasets[[nm]]
        for (suffix in suffixes)
            expect_identical(result_of(process.dataset(d, suffix = suffix)),
                             result_of(get_match_6_0_3(process.dataset(d, suffix = suffix))),
                             label = paste0(nm, ", suffix ", deparse(suffix)))
    }
    expect_error(process.dataset(no_route, suffix = "df1"),
                 "Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)",
                 fixed = TRUE)
    ## ... is ignored, as in 6.0.3, so effect_priors changes nothing
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
    ## PP.H0 is NA, where 6.0.3 gave a number
    far <- bf_a[, 5:8]
    expect_warning(got <- coloc.bf_bf(far, bf_b),
                   "snp overlap too small between datasets: too few snps with high posterior in one trait represented in other",
                   fixed = TRUE)
    expect_warning(expected <- get_match_6_0_3(coloc.bf_bf(far, bf_b)))
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

test_that("coloc.susie() refuses a dataset without snp, as 6.0.3's did, naming it", {
    no_snp <- D3
    no_snp$snp <- NULL
    expect_error(suppressMessages(get_match_6_0_3(coloc.susie(no_snp, D4))))
    expect_error(coloc.susie(no_snp, D4),
                 "dataset 1: missing required element(s) snp", fixed = TRUE)
})

## ---- the deliberate differences, one test each -------------------------------
## Each test shows the new answer and, for most, 6.0.3's value, or whether
## 6.0.3 stopped or warned.

test_that("coloc() defaults p12 to 5e-6, where coloc.abf() used 1e-5", {
    expected <- get_match_6_0_3(coloc.abf(D1, D2))
    expect_equal(expected$priors[["p12"]], 1e-5)
    expect_equal(coloc(D1, D2)$priors[["p12"]], 5e-6)
    expect_equal(as.list(pp(coloc(D1, D2, prior = coloc_prior(p12 = 1e-5)))),
                 as.list(expected$summary))
})

test_that("the priors given are used, which 6.0.3's coloc.susie_bf() ignored", {
    expected <- get_match_6_0_3(coloc.susie_bf(fit3, bf4, p12 = 1e-4))
    expect_equal(expected$priors[["p12"]], 5e-6)
    got <- coloc(fit3, bf4, prior = coloc_prior(p12 = 1e-4))
    expect_equal(got$priors[["p12"]], 1e-4)
})

test_that("coloc.susie() now applies the region-size guard, which 6.0.3's never did", {
    ## With 500 snps, p1 = p2 = 0.01 means five causal variants in each trait,
    ## so coloc(), which coloc.susie() now runs through, caps both at 1/501.
    ## 6.0.3's coloc.susie() ran through coloc.bf_bf(), which has no cap.
    expect_no_warning(expected <- get_match_6_0_3(
        coloc.susie(fit3, fit4, p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    expect_equal(expected$priors, c(p1 = 0.01, p2 = 0.01, p12 = 1e-4))
    ## one warning for each prior capped
    expect_warning(expect_warning(
        got <- coloc.susie(fit3, fit4, p1 = 0.01, p2 = 0.01, p12 = 1e-4),
        "p1 * nsnps >= 1, setting p1=1/(nsnps + 1)", fixed = TRUE),
        "p2 * nsnps >= 1, setting p2=1/(nsnps + 1)", fixed = TRUE)
    expect_equal(got$priors, c(p1 = 1 / 501, p2 = 1 / 501, p12 = 1e-4))
    expect_gt(got$summary$PP.H4.abf[1], expected$summary$PP.H4.abf[1])
    ## and nothing else moves: the same as 6.0.3 given the capped priors
    capped <- get_match_6_0_3(coloc.susie(fit3, fit4, p1 = 1 / 501, p2 = 1 / 501, p12 = 1e-4))
    expect_equal(got$summary, capped$summary)
})

test_that("a vector p1 or p2 is rejected, where 6.0.3 allowed it, with extra hypotheses", {
    ## p1 and p2 of length 8, one per snp: 6.0.3 accepted them and reported
    ## 4 * 8 + 1 = 33 hypotheses in place of five; the count grows with the
    ## length of the vectors
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, p1 = rep(1e-4, 8), p2 = rep(1e-4, 8)))
    expect_equal(sum(startsWith(names(expected$summary), "PP.H")), 33)
    expect_error(coloc(bf_a, bf_b, prior = coloc_prior(p1 = rep(1e-4, 8))),
                 "p1, p2 and p12 must each be a scalar; to vary the prior across SNPs, use weights1/weights2",
                 fixed = TRUE)
    expect_error(coloc.bf_bf(bf_a, bf_b, p1 = rep(1e-4, 8)),
                 "p1, p2 and p12 must each be a scalar; to vary the prior across SNPs, use weights1/weights2",
                 fixed = TRUE)
})

test_that("the old functions refuse a p12 outside [p1 * p2, min(p1, p2)], which 6.0.3 took", {
    ## p12 = 1e-3 is above p1 = 1e-4, outside the range sensitivity() sweeps
    expect_no_error(get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, p12 = 1e-3)))
    expect_no_error(get_match_6_0_3(coloc.susie(fit3, fit4, p12 = 1e-3)))
    expect_error(coloc.bf_bf(bf_a, bf_b, p12 = 1e-3),
                 "p12 must lie in [p1*p2, min(p1, p2)]",
                 fixed = TRUE)
    expect_error(coloc.susie(fit3, fit4, p12 = 1e-3),
                 "p12 must lie in [p1*p2, min(p1, p2)]",
                 fixed = TRUE)
    expect_error(coloc.susie_bf(fit3, bf4, p12 = 1e-3),
                 "p12 must lie in [p1*p2, min(p1, p2)]",
                 fixed = TRUE)
})

test_that("traits with no snps in common warn, where 6.0.3 returned nsnps = NA silently", {
    ## 6.0.3's coloc.bf_bf() returned early with data.table(nsnps = NA) and no warning;
    ## coloc() warns "no snps in common" and returns status "no_overlap"
    bf_x <- bf_b
    colnames(bf_x) <- paste0("x", 1:8)
    expect_silent(expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_x)))
    expect_true(is.na(expected$nsnps))
    expect_warning(got <- coloc(bf_a, bf_x),
                   "no snps in common between the two traits", fixed = TRUE)
    expect_equal(got$status, "no_overlap")
})

test_that("an all-NA signal is rejected, where 6.0.3 reported PP.H2 = 1", {
    bf_na <- bf_a
    bf_na[1, ] <- NA
    expected <- get_match_6_0_3(coloc.bf_bf(bf_na, bf_b, trim_by_posterior = FALSE))
    expect_equal(expected$summary$PP.H2.abf[1], 1)
    expect_equal(expected$summary$hit1[1], "-")
    expect_error(coloc(bf_na, bf_b),
                 "trait 1: log Bayes factors contain missing or infinite values", fixed = TRUE)
})

test_that("coloc(), given a dataset and log BFs, warns that method does not apply to the log BFs", {
    expected <- suppressMessages(get_match_6_0_3(coloc.susie_bf(D3, bf4)))
    expect_warning(got <- suppressMessages(coloc(D3, bf4, method = "susie")),
                   class = "coloc_method_ignored")
    expect_equal(got$summary, expected$summary)
})

test_that("coloc_prior() warns when both traits have weights, as coloc.abf() did", {
    w <- c(1, 1, 1, 1, 1, 1, 1, 100)
    expect_silent(expected <- get_match_6_0_3(
        coloc.bf_bf(bf_a, bf_b, prior_weights1 = w, prior_weights2 = w)))
    expect_warning(prior <- coloc_prior(weights1 = w, weights2 = w),
                   "Prior weights specified for both traits.\n  The two weight vectors should be derived from independent sources of\n  information, i.e. not computed using the same method.",
                   fixed = TRUE)
    expect_equal(coloc(bf_a, bf_b, prior = prior)$summary, expected$summary)
})

test_that("a p-value above 1 is refused, which 6.0.3 accepted", {
    above <- to_pvalues(D1)
    above$pvalues[5] <- 1.5
    expect_error(check_dataset(above, 1),
                 "dataset 1: pvalues must be above 0 and at most 1", fixed = TRUE)
    expect_error(coloc(above, D2), "dataset 1: pvalues must be above 0 and at most 1", fixed = TRUE)
    expect_error(finemap.abf(above), "pvalues must be above 0 and at most 1", fixed = TRUE)
    ## 0 and below were refused already
    zero <- to_pvalues(D1)
    zero$pvalues[5] <- 0
    expect_error(check_dataset(zero, 1),
                 "dataset 1: pvalues must be above 0 and at most 1", fixed = TRUE)
    ## exactly 1 is a p-value
    one <- to_pvalues(D1)
    one$pvalues[5] <- 1
    expect_null(check_dataset(one, 1))
})

test_that("check_dataset() reports a dataset with no type, where 6.0.3 gave R's error", {
    ## runsusie() leaves type out of req, so check_dataset() has to notice a
    ## dataset without one by itself
    d <- D3
    d$type <- NULL
    req <- c("beta", "varbeta", "LD", "snp", "N")   # runsusie()'s, without type
    expect_error(check_dataset(d, 1, req = req),
                 'dataset 1: variable type not set',
                 fixed = TRUE)
    ## 6.0.3 refused it too when type was required, as coloc.abf() asks
    expect_error(check_dataset(d, 1), "dataset 1: missing required element(s) type", fixed = TRUE)
})

test_that("a SuSiE fit with a duplicated snp id is refused, where 6.0.3 used the first", {
    ## today's coloc.susie() runs through coloc(), which checks a fit's log
    ## BFs as it checks log BFs given directly; 6.0.3's went to coloc.bf_bf()
    dup <- fit4
    colnames(dup$lbf_variable)[2] <- colnames(dup$lbf_variable)[1]
    expected <- get_match_6_0_3(coloc.susie(fit3, dup))
    expect_equal(nrow(expected$summary), 2L)
    expect_error(coloc.susie(fit3, dup),
                 "trait 2: duplicated snps found in the log Bayes factors", fixed = TRUE)
    ## an NA in a credible set's log BFs stopped 6.0.3 too, with R's own error
    na <- fit4
    na$lbf_variable[na$sets$cs_index[1], 10] <- NA
    expect_error(get_match_6_0_3(coloc.susie(fit3, na)))
    expect_error(coloc.susie(fit3, na),
                 "trait 2: log Bayes factors contain missing or infinite values", fixed = TRUE)
})

test_that("process.dataset() no longer names a one-snp result's row after its type", {
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

test_that("check_dataset()'s errors and warnings hold the dataset's number, and none without a label", {
    ## 6.0.3 named no dataset in four of its messages, among them N's and LD's,
    ## or in the minimum p-value warning. Each error and warning now begins
    ## with the dataset's number, as the caller gave it
    bad <- D2
    bad$beta[1] <- NA
    expect_error(coloc.abf(D1, bad), "dataset 2: beta contains missing values", fixed = TRUE)
    expect_error(check_dataset(bad, 2), "dataset 2: beta contains missing values", fixed = TRUE)
    bad3 <- D3
    bad3$beta[1] <- NA
    expect_error(runsusie(bad3), "dataset 1: beta contains missing values", fixed = TRUE)
    expect_error(check_dataset(1:3, 2),
                 "dataset 2: is not a list",
                 fixed = TRUE)
    ## without a label
    expect_error(finemap.abf(bad), "beta contains missing values", fixed = TRUE)
    expect_error(check_dataset(bad), "beta contains missing values", fixed = TRUE)
    ## the messages that named no dataset
    per_snp_N <- D2
    per_snp_N$N <- rep(D2$N, length(D2$snp))
    expect_error(coloc.abf(D1, per_snp_N),
                 "dataset 2: N must be a single number; coloc does not take per-snp sample sizes",
                 fixed = TRUE)
    not_square <- D2
    not_square$LD <- not_square$LD[, -1]
    expect_error(coloc.abf(D1, not_square),
                 "dataset 2: LD not square", fixed = TRUE)
    ## the minimum p-value warning, given once for each dataset
    weak <- D1
    weak$beta <- weak$beta / 100
    weak2 <- D2
    weak2$beta <- weak2$beta / 100
    expect_warning(expect_warning(coloc.abf(weak, weak2),
                                  "dataset 1: minimum p value is: 0.93383",
                                  fixed = TRUE, class = "coloc_min_p"),
                   "dataset 2: minimum p value is: 0.93863",
                   fixed = TRUE, class = "coloc_min_p")
    expect_warning(check_dataset(weak, 2),
                   "dataset 2: minimum p value is: 0.93383",
                   fixed = TRUE, class = "coloc_min_p")
    expect_warning(check_dataset(weak),
                   "minimum p value is: 0.93383",
                   fixed = TRUE, class = "coloc_min_p")
})

test_that("check_dataset() doesn't repeat the minimum p-value warning for a coloc_dataset, where 6.0.3 did", {
    ## a coloc_dataset gave the warning when it was made; a plain list still
    ## gives it each time, as test-dataset.R checks
    weak <- D1
    weak$beta <- weak$beta / 100
    made <- suppressWarnings(as_coloc_dataset(weak))
    expect_no_warning(check_dataset(made, 1))
    ## a coloc_dataset changed to hold a negative varbeta is refused before the
    ## p-values are worked out, where 6.0.3 went on to compare a NaN and stopped
    made$varbeta[1] <- -made$varbeta[1]
    expect_no_warning(expect_error(check_dataset(made, 1),
                                   "dataset 1: varbeta must be positive", fixed = TRUE))
})

test_that("coloc.susie_bf() refuses a bf2 that isn't log Bayes factors, where 6.0.3 failed with R's errors", {
    ## coloc() sorts its inputs by shape, so it would fine-map a dataset or a
    ## fit given as bf2. 6.0.3 failed on both inside coloc.bf_bf()
    expect_error(suppressMessages(get_match_6_0_3(coloc.susie_bf(fit3, D4))))
    expect_error(suppressMessages(get_match_6_0_3(coloc.susie_bf(fit3, fit4))))
    expect_error(coloc.susie_bf(fit3, D4),
                 "bf2 must be log Bayes factors: a named numeric vector or a matrix",
                 fixed = TRUE)
    expect_error(coloc.susie_bf(fit3, fit4),
                 "bf2 must be log Bayes factors: a named numeric vector or a matrix",
                 fixed = TRUE)
    chr <- bf4
    storage.mode(chr) <- "character"
    expect_error(coloc.susie_bf(fit3, chr),
                 "bf2 must be log Bayes factors: a named numeric vector or a matrix",
                 fixed = TRUE)
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
    expect_error(get_match_6_0_3(check_alignment(wider, do_plot = FALSE)))
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
