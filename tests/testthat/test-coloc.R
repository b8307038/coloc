## coloc(), the entry point: the method each trait gets; that each kind of
## input reaches its adapter; the checks before the engine runs, on the prior,
## on each dataset before either is fine-mapped, and on log BFs; and the cap on
## a prior in a large region, as coloc.abf() has it. Some tests compare with
## 6.0.3: get_match_6_0_3(), in helper-v603.R, runs a call as 6.0.3 made it.

D1 <- plain_test_data$D1
D2 <- plain_test_data$D2
D3 <- plain_test_data$D3
D4 <- plain_test_data$D4

## the warnings an expression gives, in order: each one's own class, the first
## of its classes, and its message
warnings_of <- function(expr) {
    own <- character(0)
    msg <- character(0)
    withCallingHandlers(expr, warning = function(cond) {
        own <<- c(own, class(cond)[1])
        msg <<- c(msg, conditionMessage(cond))
        invokeRestart("muffleWarning")
    })
    data.frame(class = own, message = msg)
}

test_that("an unsupplied method means single for datasets", {
    ## method's default is a vector, c("single", "susie"), not a request
    res <- coloc(D1, D2)
    expect_equal(res$meta$method, c("single", "single"))
    expect_equal(res$summary, coloc(D1, D2, method = "single")$summary)
})

test_that("method is per trait: a trait without LD can stay single", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D1_noLD <- D1
    D1_noLD$LD <- NULL
    ## SuSiE needs LD, so "susie" for both cannot run
    expect_error(suppressMessages(coloc(D1_noLD, D3, method = "susie")),
                 "dataset 1: LD is missing, and coloc needs it here", fixed = TRUE)
    mixed <- suppressMessages(coloc(D1_noLD, D3, method = c("single", "susie")))
    expect_equal(mixed$meta$method, c("single", "susie"))
    expect_equal(n_signals(mixed), 2L)
    ## "single" for both runs, but takes D3's two signals as one
    expect_equal(n_signals(coloc(D1_noLD, D3)), 1L)
})

test_that("every kind of input reaches its adapter", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    fit4 <- suppressMessages(runsusie(D4))
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
    from_fit <- suppressMessages(coloc(D3, fit4, method = "susie"))
    ## the method meant for D3 reaches the log BFs too, and warns
    expect_warning(from_bf <- suppressMessages(coloc(D3, bf4, method = "susie")),
                   "trait 2: method is inapplicable to precomputed log BFs",
                   fixed = TRUE, class = "coloc_method_ignored")
    expect_equal(from_bf$summary$PP.H4.abf, from_fit$summary$PP.H4.abf)
    expect_equal(from_bf$meta$method, c("susie", NA))
    ## a fit needs no method: with none named, it is taken as "susie". One
    ## method applies to both traits, so "single" here reaches the fit too,
    ## and stops coloc(): a fit has no beta or varbeta for "single" to use
    expect_equal(coloc(fit4, fit4)$meta$method, c("susie", "susie"))
    expect_error(coloc(fit4, D3, method = "single"),
                 'trait 1: method = "single" on a susie fit; a fit carries no beta/varbeta to run the Wakefield calculation on',
                 fixed = TRUE)
    ## a method named for log BFs warns once for each trait
    expect_equal(warnings_of(coloc(bf4, bf4, method = "single")),
                 data.frame(class = "coloc_method_ignored",
                            message = c("trait 1: method is inapplicable to precomputed log BFs",
                                        "trait 2: method is inapplicable to precomputed log BFs")))
})

test_that("a prior given as a plain list is checked for scalars and weights", {
    bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    expect_error(coloc(bf, bf, prior = list(p1 = rep(1e-4, 4), p2 = 1e-4, p12 = 5e-6)),
                 "p1, p2 and p12 must each be a scalar; to vary the prior across SNPs, use weights1/weights2",
                 fixed = TRUE)
    ## a bad weight would otherwise reach the engine, and end in an internal
    ## error that blames coloc
    expect_error(coloc(bf, bf, prior = list(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                                            weights1 = c(1, -1, 5, 1))),
                 "weights1 must be above 0 with no missing values; weight 2 is -1", fixed = TRUE)
})

test_that("coloc() checks p12's bounds in a plain list, and in a coloc_prior edited after it was made", {
    bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    expect_error(coloc(bf, bf, prior = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-3)),
                 "p12 must be at least p1 * p2 and at most the smaller of p1 and p2; it is 0.001",
                 fixed = TRUE)
    edited <- coloc_prior()
    edited$p12 <- 1e-3
    expect_error(coloc(bf, bf, prior = edited),
                 "p12 must be at least p1 * p2 and at most the smaller of p1 and p2; it is 0.001",
                 fixed = TRUE)
})

## one mistake each, made to a copy of a good dataset
mistakes <- list(
    "no type" = function(d) { d$type <- NULL; d },
    "a wrong type" = function(d) { d$type <- "quantitative"; d },
    "an NA beta" = function(d) { d$beta[5] <- NA; d },
    "a duplicated snp" = function(d) { d$snp[2] <- d$snp[1]; d },
    "mismatched lengths" = function(d) { d$beta <- d$beta[-1]; d },
    "a zero varbeta" = function(d) { d$varbeta[5] <- 0; d })

## the error each mistake gets, after its "dataset 1: " label. The SuSiE route
## doesn't require type, so a missing one gets the check on type's message
refusals <- c("no type" = "type is missing, and coloc needs it here",
              "a wrong type" = "type must be quant or cc",
              "an NA beta" = "beta contains missing values",
              "a duplicated snp" = "snp has duplicated ids, such as s1; each snp must appear once",
              "mismatched lengths" = "snp has 500 values but beta has 499",
              "a zero varbeta" = "varbeta has zeros, such as at s5; a variance must be above 0")
susie_refusals <- replace(refusals, "no type",
                          'type is missing; set type to "quant" for a quantitative trait or "cc" for a case-control trait')

## each dataset as a plain list, and made with as_coloc_dataset() before the
## mistake: $<- keeps the class, so coloc() checks the changed dataset again
kinds <- list("as a plain list" = identity,
              "as a coloc_dataset" = as_coloc_dataset)

test_that("a bad dataset stops coloc() as it stopped 6.0.3's coloc.abf(), naming the trait, for either trait", {
    for (kind in names(kinds)) for (what in names(mistakes)) for (trait in 1:2) {
        d1 <- kinds[[kind]](D1)
        d2 <- kinds[[kind]](D2)
        if (trait == 1)
            d1 <- mistakes[[what]](d1)
        else
            d2 <- mistakes[[what]](d2)
        label <- paste(what, "in trait", trait, kind)
        expect_error(get_match_6_0_3(coloc.abf(d1, d2)), label = label)
        expect_error(coloc(d1, d2), paste0("dataset ", trait, ": ", refusals[[what]]),
                     fixed = TRUE, label = label)
    }
})

test_that("on the SuSiE route a bad dataset is reported with its own trait number", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    fit3 <- suppressMessages(runsusie(D3))
    fit4 <- suppressMessages(runsusie(D4))
    for (kind in names(kinds)) for (what in names(mistakes)) {
        d3 <- mistakes[[what]](kinds[[kind]](D3))
        d4 <- mistakes[[what]](kinds[[kind]](D4))
        expect_error(suppressMessages(coloc(d3, fit4, method = "susie")),
                     paste0("dataset 1: ", susie_refusals[[what]]),
                     fixed = TRUE, label = paste(what, "in trait 1", kind))
        expect_error(suppressMessages(coloc(fit3, d4, method = "susie")),
                     paste0("dataset 2: ", susie_refusals[[what]]),
                     fixed = TRUE, label = paste(what, "in trait 2", kind))
    }
})

test_that("bad log BFs stop coloc() with the trait named", {
    good <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    bad <- good
    bad[1, 2] <- NA
    expect_error(coloc(bad, good),
                 "trait 1: log Bayes factors contain missing or infinite values", fixed = TRUE)
    expect_error(coloc(good, bad),
                 "trait 2: log Bayes factors contain missing or infinite values", fixed = TRUE)
    expect_error(coloc(c(1, 2, 3), good),
                 "trait 1: a log-BF vector must be named by snp", fixed = TRUE)
})

## ---- both datasets are checked before either is fine-mapped ------------
## coloc() passes each dataset through .new_coloc_dataset(), as coloc_dataset()
## does, with the trait's number and what the method needs, before any adapter
## runs. Its errors begin with the trait's number, as "dataset 2: ", and show
## no call.

test_that("a mistake in a dataset names the dataset, and shows no call", {
    bad <- D2
    bad$beta[1] <- NA
    err <- expect_error(coloc(D1, bad), "dataset 2: beta contains missing values", fixed = TRUE)
    expect_null(conditionCall(err))
    ## two of the refusals that named no dataset in 6.0.3
    per_snp_N <- D2
    per_snp_N$N <- rep(D2$N, length(D2$snp))
    expect_error(coloc(D1, per_snp_N),
                 "dataset 2: N must be a single number; coloc does not take per-snp sample sizes",
                 fixed = TRUE)
    not_square <- D2
    not_square$LD <- not_square$LD[, -1]
    expect_error(coloc(D1, not_square),
                 "dataset 2: LD must be square, but it has 500 rows and 499 columns", fixed = TRUE)
})

test_that("a coloc_dataset changed after it was made is checked again", {
    made <- as_coloc_dataset(D2)
    made$beta[1] <- NA                  # still a coloc_dataset: $<- keeps the class
    expect_error(coloc(as_coloc_dataset(D1), made),
                 "dataset 2: beta contains missing values", fixed = TRUE)
})

test_that("both datasets are checked before either is fine-mapped", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    no_ld <- D4
    no_ld$LD <- NULL
    ## runsusie() says "running max iterations" as it starts
    said <- character(0)
    err <- withCallingHandlers(
        tryCatch(coloc(D3, no_ld, method = "susie"), error = identity),
        message = function(m) {
            said <<- c(said, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_identical(conditionMessage(err), "dataset 2: LD is missing, and coloc needs it here")
    expect_length(said, 0)
})

test_that("a weak dataset gives one minimum p-value warning, on either route", {
    weak <- function(d) {
        d$beta <- d$beta / 100
        d
    }
    ## the first line of each minimum p-value warning
    minp <- function(w) sub("\n.*", "", w$message[w$class == "coloc_min_p"])
    expect_equal(minp(warnings_of(coloc(weak(D1), weak(D2)))),
                 c("dataset 1: the smallest p-value is 0.93383, so this dataset may have no association in the region.",
                   "dataset 2: the smallest p-value is 0.93863, so this dataset may have no association in the region."))
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    w <- warnings_of(suppressMessages(coloc(weak(D3), weak(D4), method = "susie")))
    expect_equal(minp(w),
                 c("dataset 1: the smallest p-value is 0.91098, so this dataset may have no association in the region.",
                   "dataset 2: the smallest p-value is 0.94878, so this dataset may have no association in the region."))
})

test_that("coloc.susie() and coloc.susie_bf() name the dataset, and show no call", {
    bad3 <- D3
    bad3$beta[1] <- NA
    err <- expect_error(coloc.susie(bad3, D4),
                        "dataset 1: beta contains missing values", fixed = TRUE)
    expect_null(conditionCall(err))
    bf4 <- setNames(rep(0, length(D4$snp)), D4$snp)
    err <- expect_error(coloc.susie_bf(bad3, bf4),
                        "dataset 1: beta contains missing values", fixed = TRUE)
    expect_null(conditionCall(err))
})

test_that("a factor or named type comes back as a string once checked", {
    as_factor <- D1
    as_factor$type <- factor("quant")
    expect_identical(as_coloc_dataset(as_factor)$type, "quant")
    named <- D1
    named$type <- c(trait = "quant")
    expect_identical(as_coloc_dataset(named)$type, "quant")
    expect_identical(coloc_dataset(snp = D1$snp, type = factor("cc"), beta = D1$beta,
                                   varbeta = D1$varbeta, s = 0.5)$type, "cc")
    ## a coloc_dataset is returned as it is
    made <- as_coloc_dataset(D1)
    made$type <- factor("quant")
    expect_identical(as_coloc_dataset(made)$type, factor("quant"))
    ## and coloc() gives the same results either way
    expect_equal(pp(coloc(as_factor, D2)), pp(coloc(D1, D2)))
    expect_equal(pp(coloc(made, D2)), pp(coloc(D1, D2)))
})

test_that("a data frame given as a dataset is refused, naming the trait", {
    df <- as.data.frame(D1[c("snp", "beta", "varbeta")])
    expect_error(coloc(df, D2),
                 "dataset 1: this is a data.frame, not a list of dataset elements (snp, beta, varbeta, ...); make a dataset with coloc_dataset() or as_coloc_dataset()",
                 fixed = TRUE)
})

test_that("what is neither a fit nor log BFs is checked as a dataset, naming the trait", {
    no_snp <- D1
    no_snp$snp <- NULL
    expect_error(coloc(no_snp, D2), "dataset 1: snp is missing, and coloc needs it here", fixed = TRUE)
    expect_error(coloc(D1, no_snp), "dataset 2: snp is missing, and coloc needs it here", fixed = TRUE)
    ## NULL is refused the same way whatever the version of R, where
    ## is.atomic(NULL) changed
    expect_error(coloc(NULL, D2),
                 "dataset 1: this is NULL, not a list of dataset elements (snp, beta, varbeta, ...); make a dataset with coloc_dataset() or as_coloc_dataset()",
                 fixed = TRUE)
    expect_error(coloc(D1, sum),
                 "dataset 2: this is a function, not a list of dataset elements (snp, beta, varbeta, ...); make a dataset with coloc_dataset() or as_coloc_dataset()",
                 fixed = TRUE)
})

test_that("a dataset made with coloc_dataset() is checked for check_dataset()'s default req", {
    expect_identical(.dataset_uses$single, eval(formals(check_dataset)$req))
})

test_that("the region-size guard caps a prior as 6.0.3's coloc.abf() did", {
    ## 500 snps at p1 = p2 = 0.01 imply five causal variants per trait
    old <- suppressWarnings(get_match_6_0_3(
        coloc.abf(D1, D2, p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    prior <- coloc_prior(p1 = 0.01, p2 = 0.01, p12 = 1e-4)
    w <- warnings_of(coloc(D1, D2, prior = prior))
    new <- suppressWarnings(coloc(D1, D2, prior = prior))
    ## one warning for p1 and one for p2
    expect_equal(w$message, c("p1 * nsnps >= 1, setting p1=1/(nsnps + 1)",
                              "p2 * nsnps >= 1, setting p2=1/(nsnps + 1)"))
    expect_equal(new$priors, old$priors)
    expect_equal(new$priors[["p1"]], 1 / 501)
    expect_equal(as.list(pp(new)), as.list(old$summary))
})

test_that("the region-size guard starts where nsnps * p reaches 1", {
    bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    ## 4 snps at 0.25 is exactly 1: capped at 1 / (4 + 1)
    prior <- coloc_prior(p1 = 0.25, p2 = 0.25, p12 = 0.0625)
    w <- warnings_of(coloc(bf, bf, prior = prior))
    res <- suppressWarnings(coloc(bf, bf, prior = prior))
    expect_equal(w$message, c("p1 * nsnps >= 1, setting p1=1/(nsnps + 1)",
                              "p2 * nsnps >= 1, setting p2=1/(nsnps + 1)"))
    expect_equal(res$priors[["p1"]], 1 / 5)
    ## 4 snps at 0.2 is 0.8: nothing is capped
    expect_silent(res <- coloc(bf, bf, prior = coloc_prior(p1 = 0.2, p2 = 0.2, p12 = 0.04)))
    expect_equal(res$priors[["p1"]], 0.2)
})
