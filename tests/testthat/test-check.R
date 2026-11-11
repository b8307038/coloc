library(coloc)
D1 <- coloc_test_data$D1
D3 <- coloc_test_data$D3

B1=D1
C1=D1
C1[c("type","s")]=list(type="cc",s=0.5)
B1[c("type","s")]=list(type="cc",s=2)

test_that("missing or bad required elements throws error", {
  expect_null(check_dataset(D1, req = names(D1)))
  e <- expect_error(check_dataset(D1, req = "test"))
  expect_identical(conditionMessage(e), "missing required element(s) test")
  ## type="cc" and s
  expect_null(check_dataset(C1))
  e <- expect_error(check_dataset(B1))
  expect_identical(conditionMessage(e), "s must be between 0 and 1")
  e <- expect_error(check_dataset(list(), ""))
  expect_identical(conditionMessage(e), "missing required element(s) type, snp")
  e <- expect_error(check_dataset(list(beta=1,p=2,type="blah"), ""))
  expect_identical(conditionMessage(e), "missing required element(s) snp")
})

test_that("LD matrix must have dimnames", {
  expect_null(check_ld(D3, D3$LD))
  ld_no_dimnames <- D3$LD
  attr(ld_no_dimnames, "dimnames") <- NULL
  expect_error(check_ld(D3, ld_no_dimnames),
               "LD must cover every snp",
               fixed = TRUE)
})

test_that("issue 79", {
  d1=list(snp=letters[1:5],
          position=1:5,
          N=200000,
          MAF=runif(5)/2,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="cc")
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="quant",
          MAF=rep(0.5,5))
  expect_error(coloc.abf(d1, d2),
               "dataset 2: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated",
               fixed = TRUE)
})

test_that("issue 79, with coloc() in place of coloc.abf()", {
  d1=list(snp=letters[1:5],
          position=1:5,
          N=200000,
          MAF=runif(5)/2,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="cc")
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="quant",
          MAF=rep(0.5,5))
  expect_error(coloc(d1, d2),
               "dataset 2: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated",
               fixed = TRUE)
})

test_that("issue 160", {
  d1=list(snp=letters[1:5],
          position=1:5,
          N=200000,
          MAF=runif(5)/2,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="cc")
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          N=rep(150000,5),
          type="quant",
          MAF=rep(0.5,5))
  expect_error(coloc.abf(d1, d2),
               "dataset 2: N must be a single number; coloc does not take per-snp sample sizes",
               fixed = TRUE)
})

test_that("issue 160, with coloc() in place of coloc.abf()", {
  d1=list(snp=letters[1:5],
          position=1:5,
          N=200000,
          MAF=runif(5)/2,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="cc")
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          N=rep(150000,5),
          type="quant",
          MAF=rep(0.5,5))
  expect_error(coloc(d1, d2),
               "dataset 2: N must be a single number; coloc does not take per-snp sample sizes",
               fixed = TRUE)
})

test_that("Infinite values in beta/varbeta triggers specific error", {
  D1_beta_inf <- D1
  D1_beta_inf$beta[[1]] <- Inf
  e <- expect_error(check_dataset(D1_beta_inf))
  expect_identical(conditionMessage(e), "beta has infinite values")

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- Inf
  e <- expect_error(check_dataset(D1_varbeta_inf))
  expect_identical(conditionMessage(e), "varbeta has infinite values")

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- 0
  e <- expect_error(check_dataset(D1_varbeta_inf))
  expect_identical(conditionMessage(e), "varbeta must be > 0")
})

## ---- check_dataset(), one test per check -----------------------------------
## Most of check_dataset()'s checks have a case here: a dataset that fails it,
## which check_dataset() refuses; test-refusals.R tests the others. Other cases
## are valid datasets, datasets that fail two checks, and the arguments req and
## warn.minp. Each case runs with the label 1, as coloc() checks a dataset, and
## with none, as a dataset is made. A refused case holds the message
## check_dataset() stops with, as coloc() shows it with the label 1, so the
## test shows which check refused the dataset. A case that gives a warning holds
## its first line, which says what is wrong.

## the datasets the cases are made from: D1 as given (beta and varbeta), and
## D1 as p-values, MAF and N
pv1 <- to_pvalues(D1)
## each case: the dataset, any arguments besides the dataset and its suffix,
## the message check_dataset() refuses it with, and the first line of the
## warning it gives
cases <- list(
    ## valid datasets
    "a valid dataset passes" =
        list(d = D1),
    "a dataset given as p-values passes" =
        list(d = pv1),
    "a case-control trait with beta and varbeta needs neither s nor sdY" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$sdY <- NULL
                 d
             })),
    "a dataset without LD passes" =
        list(d = local({
                 d <- D1
                 d$LD <- NULL
                 d
             })),
    "an element check_dataset() doesn't know may hold missing values" =
        list(d = local({
                 d <- D1
                 d$extra <- NA
                 d
             })),
    ## a dataset for each check
    "the dataset is a list" =
        list(d = 1:3,
             error = "dataset 1: this is a numeric vector, not a list of coloc_dataset elements (snp, beta, varbeta, ...)"),
    "the required elements are there" =
        list(d = local({
                 d <- D1
                 d$snp <- NULL
                 d
             }),
             error = "dataset 1: missing required element(s) snp"),
    "two required elements missing" =
        list(d = local({
                 d <- D1
                 d$type <- NULL
                 d$snp <- NULL
                 d
             }),
             error = "dataset 1: missing required element(s) type, snp"),
    "no missing values in beta" =
        list(d = local({
                 d <- D1
                 d$beta[1] <- NA
                 d
             }),
             error = "dataset 1: beta contains missing values"),
    "no missing values in MAF" =
        list(d = local({
                 d <- D1
                 d$MAF[1] <- NA
                 d
             }),
             error = "dataset 1: MAF contains missing values"),
    "type is quant or cc" =
        list(d = local({
                 d <- D1
                 d$type <- "binary"
                 d
             }),
             error = "dataset 1: type must be quant or cc"),
    "the snp ids are unique" =
        list(d = local({
                 d <- D1
                 d$snp[2] <- d$snp[1]
                 d
             }),
             error = "dataset 1: duplicated snps found"),
    "snp is not a factor" =
        list(d = local({
                 d <- D1
                 d$snp <- factor(d$snp)
                 d
             }),
             error = "dataset 1: snp should be a character vector but is a factor"),
    "MAF is above 0" =
        list(d = local({
                 d <- D1
                 d$MAF[1] <- 0
                 d
             }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "MAF is below 1" =
        list(d = local({
                 d <- D1
                 d$MAF[1] <- 1
                 d
             }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "MAF is numeric" =
        list(d = local({
                 d <- D1
                 d$MAF <- as.character(d$MAF)
                 d
             }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "the vectors have the same length" =
        list(d = local({
                 d <- D1
                 d$beta <- d$beta[-1]
                 d
             }),
             error = "dataset 1: snp has 500 values but beta has 499"),
    "N is a single number" =
        list(d = local({
                 d <- D1
                 d$N <- c(1000, 1000)
                 d
             }),
             error = "dataset 1: N must be a single number; coloc does not take per-snp sample sizes"),
    "N is numeric" =
        list(d = local({
                 d <- D1
                 d$N <- as.character(d$N)
                 d
             }),
             error = "dataset 1: N must be a single number; coloc does not take per-snp sample sizes"),
    "s is between 0 and 1" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$s <- 1.5
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is above 0" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$s <- 0
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is below 1" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$s <- 1
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is numeric" =
        list(d = local({
                 d <- D1
                 d$type <- "cc"
                 d$s <- "0.5"
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "beta without varbeta needs p-values and MAF" =
        list(d = local({
                 d <- D1
                 d$varbeta <- NULL
                 d
             }),
             error = "dataset 1: pvalues is missing; without beta and varbeta, coloc needs pvalues, MAF and N"),
    "without beta and varbeta, there are p-values" =
        list(d = local({
                 d <- D1
                 d$beta <- NULL
                 d$varbeta <- NULL
                 d
             }),
             error = "dataset 1: pvalues is missing; without beta and varbeta, coloc needs pvalues, MAF and N"),
    "without beta and varbeta, there is MAF" =
        list(d = local({
                 d <- pv1
                 d$MAF <- NULL
                 d
             }),
             error = "dataset 1: MAF is missing; without beta and varbeta, coloc needs pvalues, MAF and N"),
    "a case-control trait given as p-values has s" =
        list(d = local({
                 d <- pv1
                 d$type <- "cc"
                 d
             }),
             error = "dataset 1: s is missing; a case-control dataset without beta and varbeta needs s, the proportion of samples who are cases"),
    "a trait given as p-values has N" =
        list(d = local({
                 d <- pv1
                 d$N <- NULL
                 d
             }),
             error = "dataset 1: N is missing; without beta and varbeta, coloc needs N, the sample size"),
    "a trait given as p-values has N above 1" =
        list(d = local({
                 d <- pv1
                 d$N <- 1
                 d
             }),
             error = "dataset 1: N must be finite and > 1"),
    "beta and varbeta are finite" =
        list(d = local({
                 d <- D1
                 d$beta[1] <- Inf
                 d
             }),
             error = "dataset 1: beta has infinite values"),
    "varbeta has no zeros" =
        list(d = local({
                 d <- D1
                 d$varbeta[1] <- 0
                 d
             }),
             error = "dataset 1: varbeta must be > 0"),
    "varbeta has no negative values" =
        list(d = local({
                 d <- D1
                 d$varbeta[1] <- -1
                 d
             }),
             error = "dataset 1: varbeta must be > 0"),
    "a warning when no p-value is small" =
        list(d = local({
                 d <- D1
                 d$beta <- d$beta / 100
                 d
             }),
             warning = "dataset 1: the smallest p-value is 0.93383, so this dataset may have no association in the region."),
    "a warning when no p-value is small, given as p-values" =
        list(d = local({
                 d <- pv1
                 d$pvalues[] <- 0.5
                 d
             }),
             warning = "dataset 1: the smallest p-value is 0.5, so this dataset may have no association in the region."),
    "a quantitative trait without sdY has MAF" =
        list(d = local({
                 d <- D1
                 d$sdY <- NULL
                 d$MAF <- NULL
                 d
             }),
             error = "dataset 1: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    "a quantitative trait without sdY has N" =
        list(d = local({
                 d <- D1
                 d$sdY <- NULL
                 d$N <- NULL
                 d
             }),
             error = "dataset 1: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    "LD is square" =
        list(d = local({
                 d <- D1
                 d$LD <- d$LD[, -1]
                 d
             }),
             error = "dataset 1: LD must be square, but it has 500 rows and 499 columns"),
    "LD's row and column names are the same" =
        list(d = local({
                 d <- D1
                 rownames(d$LD)[1] <- "x"
                 d
             }),
             error = "dataset 1: LD rownames != colnames"),
    "LD covers every snp" =
        list(d = local({
                 d <- D1
                 d$LD <- d$LD[-1, -1]
                 d
             }),
             error = "dataset 1: LD must cover every snp"),
    ## the arguments
    "warn.minp sets when the warning is given" =
        list(d = D1, args = list(warn.minp = 1e-30),
             warning = "dataset 1: the smallest p-value is < 2.22e-16, so this dataset may have no association in the region."),
    "req sets the required elements" =
        list(d = D1, args = list(req = c("snp", "pvalues")),
             error = "dataset 1: missing required element(s) pvalues"),
    ## datasets that fail two checks
    "a wrong type and duplicated snps" =
        list(d = local({
                 d <- D1
                 d$type <- "binary"
                 d$snp[2] <- d$snp[1]
                 d
             }),
             error = "dataset 1: type must be quant or cc"),
    "an s out of range, without beta or p-values" =
        list(d = local({
                 d <- D1
                 d$beta <- NULL
                 d$varbeta <- NULL
                 d$type <- "cc"
                 d$s <- 1.5
                 d
             }),
             error = "dataset 1: s must be between 0 and 1"),
    "a bad LD stops a weak dataset before the p-value warning" =
        list(d = local({
                 d <- D1
                 d$beta <- d$beta / 100
                 d$LD <- d$LD[, -1]
                 d
             }),
             error = "dataset 1: LD must be square, but it has 500 rows and 499 columns"))

for (what in names(cases)) {
    test_that(paste0("check_dataset(): ", what), {
        case <- cases[[what]]
        for (suffix in list(1, "")) {
            check <- function() do.call(check_dataset, c(list(case$d, suffix), case$args))
            ## without a label a message has no "dataset 1: ", and the
            ## not-a-list message begins "the dataset is"
            expected <- function(m) {
                if (identical(suffix, 1))
                    return(m)
                sub("^dataset 1: ", "", sub("^dataset 1: this is", "the dataset is", m))
            }
            ## a warning comes before any refusal
            until_refused <- function() tryCatch(check(), error = function(e) NULL)
            if (is.null(case$warning)) {
                expect_no_warning(until_refused())
            } else {
                w <- expect_warning(until_refused(), class = "coloc_min_p")
                expect_identical(sub("\n.*", "", conditionMessage(w)), expected(case$warning))
            }
            if (!is.null(case$error)) {
                e <- expect_error(suppressWarnings(check()))
                expect_identical(conditionMessage(e), expected(case$error))
            } else {
                expect_null(suppressWarnings(check()))
            }
        }
    })
}

test_that("an sdY far from its estimate from varbeta, MAF and N gets a warning, as when varbeta holds standard errors", {
    ## standard errors in varbeta also make the minimum p-value warning
    no_min_p <- function(expr)
        withCallingHandlers(expr, coloc_min_p = function(w) invokeRestart("muffleWarning"))
    se <- D1
    se$varbeta <- sqrt(se$varbeta)
    expect_warning(no_min_p(check_dataset(se, 1)),
                   "dataset 1: sdY is 1.1, but varbeta, MAF and N suggest 4.55; check that varbeta holds variances, not standard errors, and that sdY is on the scale of beta",
                   fixed = TRUE, class = "coloc_sdY")
    expect_no_warning(check_dataset(D1, 1))
    ## a coloc_dataset edited after it was made, and the old functions
    made <- suppressWarnings(as_coloc_dataset(D1))
    made$varbeta <- sqrt(made$varbeta)
    w <- expect_warning(check_dataset(made), class = "coloc_sdY")
    expect_identical(conditionMessage(w),
                     "sdY is 1.1, but varbeta, MAF and N suggest 4.55; check that varbeta holds variances, not standard errors, and that sdY is on the scale of beta")
    expect_warning(no_min_p(utils::capture.output(coloc.abf(coloc_test_data$D2, se))),
                   "dataset 2: sdY is 1.1, but varbeta, MAF and N suggest 4.55; check that varbeta holds variances, not standard errors, and that sdY is on the scale of beta",
                   fixed = TRUE, class = "coloc_sdY")
    ## coloc() gives it once, with the dataset's number, though wakefield_bf()
    ## checks the dataset again
    classes <- character(0)
    withCallingHandlers(coloc(made, as_coloc_dataset(coloc_test_data$D2)),
                        warning = function(w) {
                            classes <<- c(classes, class(w)[1])
                            invokeRestart("muffleWarning")
                        })
    expect_equal(sum(classes == "coloc_sdY"), 1)
    ## not where sdY isn't read
    expect_no_warning(no_min_p(check_dataset(se, 1, req = .dataset_uses$susie)), class = "coloc_sdY")
    expect_no_warning(no_min_p(check_dataset(to_pvalues(se), 1)), class = "coloc_sdY")
})

## ---- what else check_dataset() and check_ld() must keep -----------------------
## check_dataset()'s errors and warnings show no call, as the message names the
## dataset; check_ld()'s show its own (call_of() is in helper-datasets.R)

test_that("check_dataset()'s errors and warnings show no call", {
    dup <- D1
    dup$snp[2] <- dup$snp[1]
    weak <- D1
    weak$beta <- weak$beta / 100
    per_snp_N <- D1
    per_snp_N$N <- rep(D1$N, length(D1$snp))
    e <- expect_error(check_dataset(dup, 1), "dataset 1: duplicated snps found",
                      fixed = TRUE)
    expect_null(conditionCall(e))
    w <- expect_warning(check_dataset(weak, 1),
                        "dataset 1: the smallest p-value is 0.93383, so this dataset may have no association in the region.",
                        fixed = TRUE, class = "coloc_min_p")
    expect_null(conditionCall(w))
    e <- expect_error(check_dataset(1:3, 1),
                      "dataset 1: this is a numeric vector, not a list of coloc_dataset elements (snp, beta, varbeta, ...)",
                      fixed = TRUE)
    expect_null(conditionCall(e))
    e <- expect_error(check_dataset(per_snp_N, 1),
                      "dataset 1: N must be a single number; coloc does not take per-snp sample sizes",
                      fixed = TRUE)
    expect_null(conditionCall(e))
    ## and through a function that calls it
    e <- expect_error(suppressMessages(coloc.abf(dup, D1)), "dataset 1: duplicated snps found",
                      fixed = TRUE)
    expect_null(conditionCall(e))
})

test_that("check_dataset() refuses a data.frame, as coloc() does", {
    df <- data.frame(snp = D1$snp, beta = D1$beta, varbeta = D1$varbeta,
                     type = "quant", sdY = D1$sdY)
    e <- expect_error(check_dataset(df))
    expect_identical(conditionMessage(e),
                     "the dataset is a data.frame, not a list of coloc_dataset elements (snp, beta, varbeta, ...)")
    ## one row, so type is a single value and only the data.frame test refuses it
    e <- expect_error(check_dataset(df[1, ], 1))
    expect_identical(conditionMessage(e),
                     "dataset 1: this is a data.frame, not a list of coloc_dataset elements (snp, beta, varbeta, ...)")
    ## and through a function that calls it
    e <- expect_error(suppressMessages(coloc.abf(df, D1)))
    expect_identical(conditionMessage(e),
                     "dataset 1: this is a data.frame, not a list of coloc_dataset elements (snp, beta, varbeta, ...)")
})

test_that("check_ld() takes a valid LD, and refuses the others with its own call", {
    expect_null(expect_no_warning(check_ld(D3, D3$LD)))
    names_differ <- D3$LD
    rownames(names_differ) <- paste0("r", seq_len(nrow(names_differ)))
    lds <- list(missing = list(LD = NULL, error = "LD is missing"),
                data_frame = list(LD = as.data.frame(D3$LD), error = "LD must be a matrix"),
                not_square = list(LD = D3$LD[, -1],
                                  error = "LD must be square, but it has 500 rows and 499 columns"),
                names_differ = list(LD = names_differ, error = "LD rownames != colnames"),
                snp_missing = list(LD = D3$LD[-1, -1], error = "LD must cover every snp"))
    for (nm in names(lds)) {
        LD <- lds[[nm]]$LD
        expect_error(check_ld(D3, LD), lds[[nm]]$error, fixed = TRUE,
                     label = nm)
        expect_identical(call_of(check_ld(D3, LD)), quote(check_ld(D3, LD)), label = nm)
    }
})

test_that("check_ld(), and so finemap.signals() and coloc.signals(), refuse an LD that isn't a numeric correlation matrix", {
    character_ld <- D3$LD
    storage.mode(character_ld) <- "character"
    ## 1.5 times too large off the diagonal
    too_large <- D3$LD * 1.5
    diag(too_large) <- 1
    with_na <- D3$LD
    with_na[1, 2] <- NA
    with_na[2, 1] <- NA
    lds <- list(character = list(LD = character_ld, error = "LD must be a numeric matrix"),
                too_large = list(LD = too_large,
                                 error = "LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal"),
                with_na = list(LD = with_na, error = "LD contains missing values"))
    for (nm in names(lds)) {
        expect_error(check_ld(D3, lds[[nm]]$LD), lds[[nm]]$error, fixed = TRUE,
                     label = nm)
        expect_error(finemap.signals(D3, lds[[nm]]$LD, method = "cond"), lds[[nm]]$error,
                     fixed = TRUE, label = nm)
    }
    ## coloc.signals() takes LD from the dataset, whose own check refuses the NA
    for (nm in c("character", "too_large")) {
        d <- D3
        d$LD <- lds[[nm]]$LD
        expect_error(coloc.signals(d, coloc_test_data$D4, method = "cond", p12 = 1e-5),
                     lds[[nm]]$error, fixed = TRUE, label = nm)
    }
})
