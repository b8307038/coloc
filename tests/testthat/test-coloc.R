## coloc(), the entry point, tested through its own steps:
## - which method each trait gets: "single" when none is given, and one
##   method per trait when two are given;
## - that each kind of input (a dataset, a SuSiE fit, log BFs) reaches the
##   adapter for it;
## - the checks made before the engine runs: a prior given as a plain list
##   must hold scalar p1, p2 and p12, and any weights1 and weights2 must be
##   positive and not NA (vectors, one per snp, as in coloc_prior()); every
##   dataset is checked, with the trait's number and coloc()'s call, before
##   either is fine-mapped; bad log BFs stop coloc() with the trait named;
## - the region-size guard, which caps p1, p2 and p12 when a prior times its
##   number of snps reaches 1, as coloc.abf() did.
## Some tests compare with 6.0.3: get_match_6_0_3(), in helper-v603.R, runs a
## call as 6.0.3 made it.

D1 <- plain_test_data$D1
D2 <- plain_test_data$D2
D3 <- plain_test_data$D3
D4 <- plain_test_data$D4

## the messages of every warning an expression gives, in order
warnings_of <- function(expr) {
    w <- character(0)
    withCallingHandlers(expr, warning = function(cond) {
        w <<- c(w, conditionMessage(cond))
        invokeRestart("muffleWarning")
    })
    w
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
                 "missing required element")
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
    ## AL14: the method meant for D3 reaches the log BFs too, and warns
    expect_warning(from_bf <- suppressMessages(coloc(D3, bf4, method = "susie")),
                   "trait 2: method is inapplicable")
    expect_equal(from_bf$summary$PP.H4.abf, from_fit$summary$PP.H4.abf)
    expect_equal(from_bf$meta$method, c("susie", NA))
    ## a fit needs no method: with none named, it is taken as "susie". One
    ## method applies to both traits, so "single" here reaches the fit too,
    ## and stops coloc(): a fit has no beta or varbeta for "single" to use
    expect_equal(coloc(fit4, fit4)$meta$method, c("susie", "susie"))
    expect_error(coloc(fit4, D3, method = "single"), "susie fit")
    ## a method named for log BFs warns once for each trait
    expect_length(warnings_of(coloc(bf4, bf4, method = "single")), 2L)
})

test_that("a prior given as a plain list is checked for scalars and weights", {
    bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    expect_error(coloc(bf, bf, prior = list(p1 = rep(1e-4, 4), p2 = 1e-4, p12 = 5e-6)),
                 "must each be a scalar")
    ## a bad weight would otherwise reach the engine, and end in an internal
    ## error that blames coloc
    expect_error(coloc(bf, bf, prior = list(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                                            weights1 = c(1, -1, 5, 1))),
                 "prior_weights1 must contain non-NA/non-negative weights")
})

## one mistake each, made to a copy of a good dataset
mistakes <- list(
    "no type" = function(d) { d$type <- NULL; d },
    "a wrong type" = function(d) { d$type <- "quantitative"; d },
    "an NA beta" = function(d) { d$beta[5] <- NA; d },
    "a duplicated snp" = function(d) { d$snp[2] <- d$snp[1]; d },
    "mismatched lengths" = function(d) { d$beta <- d$beta[-1]; d },
    "a zero varbeta" = function(d) { d$varbeta[5] <- 0; d })

## each dataset as a plain list, and made with as_coloc_dataset() before the
## mistake: $<- keeps the class, so coloc() checks the changed dataset again
kinds <- list("as a plain list" = identity,
              "as a coloc_dataset" = as_coloc_dataset)

test_that("a bad dataset stops coloc() with 6.0.3's coloc.abf() message, naming the trait, for either trait", {
    ## the lengths message is reworded since 6.0.3
    reworded <- c("mismatched lengths" = "snp has 500 values but beta has 499")
    for (kind in names(kinds)) for (what in names(mistakes)) for (trait in 1:2) {
        d1 <- kinds[[kind]](D1)
        d2 <- kinds[[kind]](D2)
        if (trait == 1) d1 <- mistakes[[what]](d1) else d2 <- mistakes[[what]](d2)
        expected <- tryCatch(get_match_6_0_3(coloc.abf(d1, d2)), error = conditionMessage)
        expect_match(expected, paste("dataset", trait), fixed = TRUE,
                     label = paste(what, "in trait", trait, kind))
        expected <- if (what %in% names(reworded)) paste0("trait ", trait, ": ", reworded[[what]])
                    else reprefix(expected, trait, word = "trait")
        expect_error(coloc(d1, d2), expected, fixed = TRUE,
                     label = paste(what, "in trait", trait, kind))
    }
})

test_that("on the SuSiE route a bad dataset is reported with its own trait number", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    fit3 <- suppressMessages(runsusie(D3))
    fit4 <- suppressMessages(runsusie(D4))
    problem <- c("no type" = "variable type not set",
                 "a wrong type" = "type must be quant or cc",
                 "an NA beta" = "beta contains missing values",
                 "a duplicated snp" = "duplicated snps found",
                 "mismatched lengths" = "snp has 500 values but beta has 499",
                 "a zero varbeta" = "Zero values in varbeta")
    for (kind in names(kinds)) for (what in names(mistakes)) {
        d3 <- mistakes[[what]](kinds[[kind]](D3))
        d4 <- mistakes[[what]](kinds[[kind]](D4))
        expect_error(suppressMessages(coloc(d3, fit4, method = "susie")),
                     paste0("trait 1: ", problem[[what]]), fixed = TRUE,
                     label = paste(what, "in trait 1", kind))
        expect_error(suppressMessages(coloc(fit3, d4, method = "susie")),
                     paste0("trait 2: ", problem[[what]]), fixed = TRUE,
                     label = paste(what, "in trait 2", kind))
    }
})

test_that("bad log BFs stop coloc() with the trait named", {
    good <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    bad <- good
    bad[1, 2] <- NA
    expect_error(coloc(bad, good), "trait 1: log Bayes factors contain missing")
    expect_error(coloc(good, bad), "trait 2: log Bayes factors contain missing")
    expect_error(coloc(c(1, 2, 3), good), "trait 1: a log-BF vector must be named by snp")
})

## ---- the door: both datasets checked before either is fine-mapped -------
## coloc() passes each dataset through .new_coloc_dataset(), as coloc_dataset()
## does, with the trait's number and what the method needs, before any adapter
## runs. Its messages begin "trait N: " and carry coloc()'s call.

test_that("a mistake in a dataset names the trait and coloc()'s call", {
    bad <- D2
    bad$beta[1] <- NA
    err <- tryCatch(coloc(D1, bad), error = identity)
    expect_identical(conditionMessage(err), "trait 2: beta contains missing values")
    expect_identical(conditionCall(err), quote(coloc(D1, bad)))
    ## two of the messages that named no dataset in 6.0.3
    per_snp_N <- D2
    per_snp_N$N <- rep(D2$N, length(D2$snp))
    expect_error(coloc(D1, per_snp_N), "trait 2: N must be a single number", fixed = TRUE)
    not_square <- D2
    not_square$LD <- not_square$LD[, -1]
    expect_error(coloc(D1, not_square), "trait 2: LD not square", fixed = TRUE)
})

test_that("a coloc_dataset changed after it was made is checked again", {
    made <- as_coloc_dataset(D2)
    made$beta[1] <- NA                  # still a coloc_dataset: $<- keeps the class
    expect_error(coloc(as_coloc_dataset(D1), made),
                 "trait 2: beta contains missing values", fixed = TRUE)
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
    expect_identical(conditionMessage(err), "trait 2: missing required element(s) LD")
    expect_length(said, 0)
})

test_that("a weak dataset gives one minimum p-value warning, on either route", {
    weak <- function(d) {
        d$beta <- d$beta / 100
        d
    }
    minp <- function(w) grep("minimum p value", sub("\n.*", "", w), value = TRUE)
    expect_identical(minp(warnings_of(coloc(weak(D1), weak(D2)))),
                     c("trait 1: minimum p value is: 0.93383",
                       "trait 2: minimum p value is: 0.93863"))
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    w <- minp(warnings_of(suppressMessages(coloc(weak(D3), weak(D4), method = "susie"))))
    expect_length(w, 2)
    expect_match(w[1], "^trait 1: minimum p value is: ")
    expect_match(w[2], "^trait 2: minimum p value is: ")
})

test_that("coloc.susie() and coloc.susie_bf() keep \"dataset N: \" and show their own call", {
    bad3 <- D3
    bad3$beta[1] <- NA
    err <- tryCatch(coloc.susie(bad3, D4), error = identity)
    expect_identical(conditionMessage(err), "dataset 1: beta contains missing values")
    expect_identical(conditionCall(err), quote(coloc.susie(bad3, D4)))
    bf4 <- setNames(rep(0, length(D4$snp)), D4$snp)
    err <- tryCatch(coloc.susie_bf(bad3, bf4), error = identity)
    expect_identical(conditionMessage(err), "dataset 1: beta contains missing values")
    expect_identical(conditionCall(err), quote(coloc.susie_bf(bad3, bf4)))
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
                 "trait 1: expected a coloc dataset: a list, or the result of coloc_dataset(), not a data.frame",
                 fixed = TRUE)
})

test_that("what is neither a fit nor log BFs is checked as a dataset, naming the trait", {
    no_snp <- D1
    no_snp$snp <- NULL
    expect_error(coloc(no_snp, D2), "trait 1: missing required element(s) snp", fixed = TRUE)
    expect_error(coloc(D1, no_snp), "trait 2: missing required element(s) snp", fixed = TRUE)
    ## NULL is refused the same way whatever the version of R, where
    ## is.atomic(NULL) changed
    expect_error(coloc(NULL, D2),
                 "trait 1: expected a coloc dataset: a list, or the result of coloc_dataset(), not a NULL",
                 fixed = TRUE)
    expect_error(coloc(D1, sum),
                 "trait 2: expected a coloc dataset: a list, or the result of coloc_dataset(), not a function",
                 fixed = TRUE)
})

test_that("a dataset made with coloc_dataset() is checked for check_dataset()'s default req", {
    expect_identical(.dataset_uses$single, eval(formals(check_dataset)$req))
})

test_that("the region-size guard caps a prior as 6.0.3's coloc.abf() did", {
    ## 500 snps at p1 = p2 = 0.01 imply five causal variants per trait
    old_w <- warnings_of(old <- get_match_6_0_3(
        coloc.abf(D1, D2, p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    new_w <- warnings_of(new <- coloc(D1, D2,
                                      prior = coloc_prior(p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    expect_equal(new_w, old_w)
    expect_equal(new$priors, old$priors)
    expect_equal(new$priors[["p1"]], 1 / 501)
    expect_equal(as.list(pp(new)), as.list(old$summary))
})

test_that("the region-size guard starts where nsnps * p reaches 1", {
    bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    ## 4 snps at 0.25 is exactly 1: capped at 1 / (4 + 1)
    w <- warnings_of(res <- coloc(bf, bf, prior = coloc_prior(p1 = 0.25, p2 = 0.25,
                                                               p12 = 0.0625)))
    expect_equal(w, c("p1 * nsnps >= 1, setting p1=1/(nsnps + 1)",
                      "p2 * nsnps >= 1, setting p2=1/(nsnps + 1)"))
    expect_equal(res$priors[["p1"]], 1 / 5)
    ## 4 snps at 0.2 is 0.8: nothing is capped
    expect_silent(res <- coloc(bf, bf, prior = coloc_prior(p1 = 0.2, p2 = 0.2, p12 = 0.04)))
    expect_equal(res$priors[["p1"]], 0.2)
})
