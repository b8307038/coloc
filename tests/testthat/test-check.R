library(coloc)
D1 <- plain_test_data$D1
D3 <- plain_test_data$D3

B1=C1=D1
C1[c("type","s")]=list(type="cc",s=0.5)
B1[c("type","s")]=list(type="cc",s=2)

test_that("missing or bad required elements throws error", {
  expect_null(check_dataset(D1, req = names(D1)))
  expect_error(check_dataset(D1, req = "test"),
               "missing required element(s) test",
               fixed = TRUE)
  ## type="cc" and s
  expect_null(check_dataset(C1))
  expect_error(check_dataset(B1),
               "s must be between 0 and 1",
               fixed = TRUE)
  expect_error(check_dataset(list(), ""),
               "missing required element(s) type, snp",
               fixed = TRUE)
  expect_error(check_dataset(list(beta=1,p=2,type="blah"), ""),
               "missing required element(s) snp",
               fixed = TRUE)
})

test_that("LD matrix must have dimnames", {
  expect_null(check_ld(D3, D3$LD))
  ld_no_dimnames <- D3$LD
  attr(ld_no_dimnames, "dimnames") <- NULL
  expect_error(check_ld(D3, ld_no_dimnames),
               "colnames in LD do not contain all SNPs",
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
               "trait 2: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated",
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
               "trait 2: N must be a single number; coloc does not take per-snp sample sizes",
               fixed = TRUE)
})

test_that("Infinite values in beta/varbeta triggers specific error", {
  D1_beta_inf <- D1
  D1_beta_inf$beta[[1]] <- Inf
  expect_error(check_dataset(D1_beta_inf),
               "Infinite values in beta and/or varbeta",
               fixed = TRUE)

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- Inf
  expect_error(check_dataset(D1_varbeta_inf),
               "Infinite values in beta and/or varbeta",
               fixed = TRUE)

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- 0
  expect_error(check_dataset(D1_varbeta_inf),
               "Zero values in varbeta",
               fixed = TRUE)
})

## ---- check_dataset(), one test per rule ------------------------------------
## Each rule of check_dataset() has a case for each thing it tests: a dataset
## that breaks it, which check_dataset() refuses. Other cases are valid
## datasets, datasets that break two rules, and the arguments req and
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
        list(d = local({ d <- D1; d$type <- "cc"; d$sdY <- NULL; d })),
    "a dataset without LD passes" =
        list(d = local({ d <- D1; d$LD <- NULL; d })),
    "an element check_dataset() doesn't know may hold missing values" =
        list(d = local({ d <- D1; d$extra <- NA; d })),
    ## the rules, in the order the function checks them
    "the dataset is a list" =
        list(d = 1:3,
             error = "dataset 1: is not a list"),
    "the required elements are there" =
        list(d = local({ d <- D1; d$snp <- NULL; d }),
             error = "dataset 1: missing required element(s) snp"),
    "two required elements missing" =
        list(d = local({ d <- D1; d$type <- NULL; d$snp <- NULL; d }),
             error = "dataset 1: missing required element(s) type, snp"),
    "no missing values in beta" =
        list(d = local({ d <- D1; d$beta[1] <- NA; d }),
             error = "dataset 1: beta contains missing values"),
    "no missing values in MAF" =
        list(d = local({ d <- D1; d$MAF[1] <- NA; d }),
             error = "dataset 1: MAF contains missing values"),
    "type is quant or cc" =
        list(d = local({ d <- D1; d$type <- "binary"; d }),
             error = "dataset 1: type must be quant or cc"),
    "the snp ids are unique" =
        list(d = local({ d <- D1; d$snp[2] <- d$snp[1]; d }),
             error = "dataset 1: duplicated snps found"),
    "snp is not a factor" =
        list(d = local({ d <- D1; d$snp <- factor(d$snp); d }),
             error = "dataset 1: snp should be a character vector but is a factor"),
    "MAF is above 0" =
        list(d = local({ d <- D1; d$MAF[1] <- 0; d }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "MAF is below 1" =
        list(d = local({ d <- D1; d$MAF[1] <- 1; d }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "MAF is numeric" =
        list(d = local({ d <- D1; d$MAF <- as.character(d$MAF); d }),
             error = "dataset 1: MAF should be a numeric, strictly >0 & <1"),
    "the vectors have the same length" =
        list(d = local({ d <- D1; d$beta <- d$beta[-1]; d }),
             error = "dataset 1: snp has 500 values but beta has 499"),
    "N is a single number" =
        list(d = local({ d <- D1; d$N <- c(1000, 1000); d }),
             error = "dataset 1: N must be a single number; coloc does not take per-snp sample sizes"),
    "N is numeric" =
        list(d = local({ d <- D1; d$N <- as.character(d$N); d }),
             error = "dataset 1: N must be a single number; coloc does not take per-snp sample sizes"),
    "s is between 0 and 1" =
        list(d = local({ d <- D1; d$type <- "cc"; d$s <- 1.5; d }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is above 0" =
        list(d = local({ d <- D1; d$type <- "cc"; d$s <- 0; d }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is below 1" =
        list(d = local({ d <- D1; d$type <- "cc"; d$s <- 1; d }),
             error = "dataset 1: s must be between 0 and 1"),
    "s is numeric" =
        list(d = local({ d <- D1; d$type <- "cc"; d$s <- "0.5"; d }),
             error = "dataset 1: s must be between 0 and 1"),
    "beta without varbeta needs p-values and MAF" =
        list(d = local({ d <- D1; d$varbeta <- NULL; d }),
             error = "dataset 1: require p values and MAF if beta, varbeta are unavailable"),
    "without beta and varbeta, there are p-values" =
        list(d = local({ d <- D1; d$beta <- NULL; d$varbeta <- NULL; d }),
             error = "dataset 1: require p values and MAF if beta, varbeta are unavailable"),
    "without beta and varbeta, there is MAF" =
        list(d = local({ d <- pv1; d$MAF <- NULL; d }),
             error = "dataset 1: require p values and MAF if beta, varbeta are unavailable"),
    "a case-control trait given as p-values has s" =
        list(d = local({ d <- pv1; d$type <- "cc"; d }),
             error = "dataset 1: require, s, proportion of samples who are cases, if beta, varbeta are unavailable"),
    "a trait given as p-values has N" =
        list(d = local({ d <- pv1; d$N <- NULL; d }),
             error = "dataset 1: sample size N <=0 or not set"),
    "a trait given as p-values has N above 0" =
        list(d = local({ d <- pv1; d$N <- 0; d }),
             error = "dataset 1: sample size N <=0 or not set"),
    "beta and varbeta are finite" =
        list(d = local({ d <- D1; d$beta[1] <- Inf; d }),
             error = "dataset 1: Infinite values in beta and/or varbeta"),
    "varbeta has no zeros" =
        list(d = local({ d <- D1; d$varbeta[1] <- 0; d }),
             error = "dataset 1: Zero values in varbeta"),
    "a warning when no p-value is small" =
        list(d = local({ d <- D1; d$beta <- d$beta / 100; d }),
             warning = "dataset 1: minimum p value is: 0.93383"),
    "a warning when no p-value is small, given as p-values" =
        list(d = local({ d <- pv1; d$pvalues[] <- 0.5; d }),
             warning = "dataset 1: minimum p value is: 0.5"),
    "a quantitative trait without sdY has MAF" =
        list(d = local({ d <- D1; d$sdY <- NULL; d$MAF <- NULL; d }),
             error = "dataset 1: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    "a quantitative trait without sdY has N" =
        list(d = local({ d <- D1; d$sdY <- NULL; d$N <- NULL; d }),
             error = "dataset 1: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    "LD is square" =
        list(d = local({ d <- D1; d$LD <- d$LD[, -1]; d }),
             error = "dataset 1: LD not square"),
    "LD's row and column names are the same" =
        list(d = local({ d <- D1; rownames(d$LD)[1] <- "x"; d }),
             error = "dataset 1: LD rownames != colnames"),
    "LD covers every snp" =
        list(d = local({ d <- D1; d$LD <- d$LD[-1, -1]; d }),
             error = "dataset 1: colnames in LD do not contain all SNPs"),
    ## the arguments
    "warn.minp sets when the warning is given" =
        list(d = D1, args = list(warn.minp = 1e-30),
             warning = "dataset 1: minimum p value is: < 2.22e-16"),
    "req sets the required elements" =
        list(d = D1, args = list(req = c("snp", "pvalues")),
             error = "dataset 1: missing required element(s) pvalues"),
    ## datasets that break two rules
    "a wrong type and duplicated snps" =
        list(d = local({ d <- D1; d$type <- "binary"
                         d$snp[2] <- d$snp[1]; d }),
             error = "dataset 1: type must be quant or cc"),
    "an s out of range, without beta or p-values" =
        list(d = local({ d <- D1; d$beta <- NULL; d$varbeta <- NULL
                         d$type <- "cc"; d$s <- 1.5; d }),
             error = "dataset 1: s must be between 0 and 1"),
    "the p-value warning comes before the LD checks" =
        list(d = local({ d <- D1; d$beta <- d$beta / 100
                         d$LD <- d$LD[, -1]; d }),
             error = "dataset 1: LD not square",
             warning = "dataset 1: minimum p value is: 0.93383"))

for (what in names(cases)) {
    test_that(paste0("check_dataset(): ", what), {
        case <- cases[[what]]
        for (suffix in list(1, "")) {
            check <- function() do.call(check_dataset, c(list(case$d, suffix), case$args))
            ## a warning comes before any refusal
            until_refused <- function() tryCatch(check(), error = function(e) NULL)
            if (is.null(case$warning)) {
                expect_no_warning(until_refused())
            } else {
                w <- expect_warning(until_refused(), class = "coloc_min_p")
                if (identical(suffix, 1))
                    expect_match(conditionMessage(w), case$warning, fixed = TRUE)
            }
            if (!is.null(case$error)) {
                e <- expect_error(suppressWarnings(check()))
                ## without a label, the not-a-list message changes form
                if (identical(suffix, 1))
                    expect_identical(conditionMessage(e), case$error)
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
                   "dataset 1: sdY is 1.1, but varbeta, MAF and N suggest 4.55: does varbeta hold variances?",
                   fixed = TRUE, class = "coloc_sdY")
    expect_no_warning(check_dataset(D1, 1))
    ## a coloc_dataset edited after it was made, and the old functions
    made <- suppressWarnings(as_coloc_dataset(D1))
    made$varbeta <- sqrt(made$varbeta)
    expect_warning(check_dataset(made),
                   "sdY is 1.1, but varbeta, MAF and N suggest 4.55: does varbeta hold variances?",
                   fixed = TRUE, class = "coloc_sdY")
    expect_warning(no_min_p(utils::capture.output(coloc.abf(plain_test_data$D2, se))),
                   "dataset 2: sdY is 1.1, but varbeta, MAF and N suggest 4.55: does varbeta hold variances?",
                   fixed = TRUE, class = "coloc_sdY")
    ## not where sdY isn't read
    expect_no_warning(no_min_p(check_dataset(se, 1, req = .dataset_uses$susie)), class = "coloc_sdY")
    expect_no_warning(no_min_p(check_dataset(to_pvalues(se), 1)), class = "coloc_sdY")
})

## ---- what else check_dataset() and check_ld() must keep -----------------------
## The rules run in .gate() (R/check.R), but an error or warning shows the call
## the user made, as 6.0.3's did (call_of() is in helper-datasets.R)

test_that("check_dataset() errors and warns with its own call", {
    dup <- D1
    dup$snp[2] <- dup$snp[1]
    weak <- D1
    weak$beta <- weak$beta / 100
    per_snp_N <- D1
    per_snp_N$N <- rep(D1$N, length(D1$snp))
    expect_identical(call_of(check_dataset(dup, 1)), quote(check_dataset(dup, 1)))
    expect_identical(call_of(check_dataset(weak, 1)), quote(check_dataset(weak, 1)))
    expect_identical(call_of(check_dataset(1:3, 1)), quote(check_dataset(1:3, 1)))
    expect_identical(call_of(check_dataset(per_snp_N, 1)), quote(check_dataset(per_snp_N, 1)))
    ## and through a function that calls it
    expect_identical(call_of(suppressMessages(coloc.abf(dup, D1))),
                     quote(check_dataset(d = dataset1, 1)))
})

test_that("check_ld() takes a valid LD, and refuses the others with its own call", {
    expect_null(expect_no_warning(check_ld(D3, D3$LD)))
    names_differ <- D3$LD
    rownames(names_differ) <- paste0("r", seq_len(nrow(names_differ)))
    lds <- list(missing = list(LD = NULL, error = "LD required"),
                data_frame = list(LD = as.data.frame(D3$LD), error = "LD must be of class matrix"),
                not_square = list(LD = D3$LD[, -1],
                                  error = "LD not square"),
                names_differ = list(LD = names_differ, error = "LD rownames != colnames"),
                snp_missing = list(LD = D3$LD[-1, -1], error = "colnames in LD do not contain all SNPs"))
    for (nm in names(lds)) {
        LD <- lds[[nm]]$LD
        expect_error(check_ld(D3, LD), lds[[nm]]$error, fixed = TRUE,
                     label = nm)
        expect_identical(call_of(check_ld(D3, LD)), quote(check_ld(D3, LD)), label = nm)
    }
})
