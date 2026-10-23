library(coloc)
D1 <- plain_test_data$D1
D3 <- plain_test_data$D3

B1=C1=D1
C1[c("type","s")]=list(type="cc",s=0.5)
B1[c("type","s")]=list(type="cc",s=2)

test_that("missing or bad required elements throws error", {
  expect_null(check_dataset(D1, req = names(D1)))
  expect_error(check_dataset(D1, req = "test"))
  ## type="cc" and s
  expect_null(check_dataset(C1))
  expect_error(check_dataset(B1))
  expect_error(check_dataset(list(), ""))
  expect_error(check_dataset(list(beta=1,p=2,type="blah"), ""))
})

test_that("LD matrix must have dimnames", {
  expect_null(check_ld(D3, D3$LD))
  ld_no_dimnames <- D3$LD
  attr(ld_no_dimnames, "dimnames") <- NULL
  expect_error(check_ld(D3, ld_no_dimnames))
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
  expect_error(coloc.abf(d1,d2), "must give sdY")
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
  expect_error(coloc(d1,d2), "must give sdY")
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
  expect_error(coloc.abf(d1,d2), "per-snp sample sizes")
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
  expect_error(coloc(d1,d2), "per-snp sample sizes")
})

test_that("Infinite values in beta/varbeta triggers specific error", {
  D1_beta_inf <- D1
  D1_beta_inf$beta[[1]] <- Inf
  expect_error(check_dataset(D1_beta_inf), "Infinite")

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- Inf
  expect_error(check_dataset(D1_varbeta_inf), "Infinite")

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- 0
  expect_error(check_dataset(D1_varbeta_inf), "Zero")
})

## ---- check_dataset(), today's against 6.0.3's, one test per rule ------------
## The data stage replaces check_dataset() with a list of rules. Before it does,
## each rule of today's function has a case for each thing it tests: a dataset
## that breaks it, checked by today's function and by 6.0.3's copy in
## helper-v603.R, which must give the same value, error and warnings. Other
## cases are valid datasets, datasets that break two rules (these fix the order
## the rules run in), and the arguments req and warn.minp. Each case runs with
## the label 1, as coloc() checks a dataset, and with none, as a dataset is
## made. Its message is 6.0.3's; today's reads as reprefix() words it: nothing
## for "dataset : ", and "dataset 1: " on the four messages that named no
## dataset and on the minimum p-value warning. Two sentences are reworded,
## for lengths that differ and for a per-snp N, and their cases give today's
## sentence as now, without its prefix. Three more things differ from 6.0.3 on
## purpose, each tested in test-v603.R: the range of p-values, above 0 and at
## most 1; a dataset with no type when type isn't required, as on the SuSiE
## route; and the p-value warning, which a coloc_dataset gives when it is made
## and not again (from PR #2).

## the datasets the cases are made from: D1 as given (beta and varbeta), and
## D1 as p-values, MAF and N
pv1 <- to_pvalues(D1)
## each case: the dataset, any arguments besides the dataset and its suffix,
## the error or warning 6.0.3 gives, and for a reworded message, now
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
        list(d = 1:3, error = "dataset 1: is not a list"),
    "the required elements are there" =
        list(d = local({ d <- D1; d$snp <- NULL; d }),
             error = "dataset 1: missing required element(s) snp"),
    "every missing required element is named" =
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
             error = "dataset 1: lengths of inputs don't match: ",
             now = "snp has 500 values but beta has 499"),
    "N is a single number" =
        list(d = local({ d <- D1; d$N <- c(1000, 1000); d }),
             error = "N should be a single numeric. coloc does not accomodate per-snp sample sizes.",
             now = "N must be a single number; coloc does not take per-snp sample sizes"),
    "N is numeric" =
        list(d = local({ d <- D1; d$N <- as.character(d$N); d }),
             error = "N should be a single numeric. coloc does not accomodate per-snp sample sizes.",
             now = "N must be a single number; coloc does not take per-snp sample sizes"),
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
             warning = "minimum p value is: "),
    "a warning when no p-value is small, given as p-values" =
        list(d = local({ d <- pv1; d$pvalues[] <- 0.5; d }),
             warning = "minimum p value is: "),
    "a quantitative trait without sdY has MAF" =
        list(d = local({ d <- D1; d$sdY <- NULL; d$MAF <- NULL; d }),
             error = "dataset 1: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    "a quantitative trait without sdY has N" =
        list(d = local({ d <- D1; d$sdY <- NULL; d$N <- NULL; d }),
             error = "dataset 1: must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated"),
    "LD is square" =
        list(d = local({ d <- D1; d$LD <- d$LD[, -1]; d }),
             error = "LD not square"),
    "LD's row and column names are the same" =
        list(d = local({ d <- D1; rownames(d$LD)[1] <- "x"; d }),
             error = "LD rownames != colnames"),
    "LD covers every snp" =
        list(d = local({ d <- D1; d$LD <- d$LD[-1, -1]; d }),
             error = "colnames in LD do not contain all SNPs"),
    ## the arguments
    "warn.minp sets when the warning is given" =
        list(d = D1, args = list(warn.minp = 1e-30),
             warning = "minimum p value is: "),
    "req sets the required elements" =
        list(d = D1, args = list(req = c("snp", "pvalues")),
             error = "dataset 1: missing required element(s) pvalues"),
    ## datasets that break two rules, where the rule checked first decides
    "type is checked before the snp ids" =
        list(d = local({ d <- D1; d$type <- "binary"
                         d$snp[2] <- d$snp[1]; d }),
             error = "dataset 1: type must be quant or cc"),
    "s is checked before the p-values a dataset without beta needs" =
        list(d = local({ d <- D1; d$beta <- NULL; d$varbeta <- NULL
                         d$type <- "cc"; d$s <- 1.5; d }),
             error = "dataset 1: s must be between 0 and 1"),
    "the p-value warning comes before the LD checks" =
        list(d = local({ d <- D1; d$beta <- d$beta / 100
                         d$LD <- d$LD[, -1]; d }),
             error = "LD not square", warning = "minimum p value is: "))

## today's message for 6.0.3's: through reprefix(), or for a reworded case its
## own sentence with the prefix
today <- function(case, then, suffix) {
    if (is.null(case$now))
        return(reprefix(then, suffix))
    paste0(if (suffix == "") "" else paste0("dataset ", suffix, ": "), case$now)
}

for (what in names(cases)) {
    test_that(paste0("check_dataset(), today's and 6.0.3's: ", what), {
        case <- cases[[what]]
        ## the same value, error and warnings as 6.0.3's, with the prefix as
        ## check_dataset() gives it now, with a label and without
        for (suffix in list(1, "")) {
            args <- c(list(case$d, suffix), case$args)
            now <- outcome(do.call(check_dataset, args))
            then <- outcome(get_match_6_0_3(do.call(check_dataset, args)))
            then$error <- today(case, then$error, suffix)
            then$warned <- reprefix(then$warned, suffix)
            expect_identical(now, then, label = paste0("suffix \"", suffix, "\""))
        }
        ## and 6.0.3 stops at the rule the case is meant to break, or passes,
        ## with the warning it is meant to give
        then <- outcome(get_match_6_0_3(do.call(check_dataset, c(list(case$d, 1), case$args))))
        if (is.null(case$error)) {
            expect_identical(then$error, NA_character_)
        } else {
            expect_identical(then$error, case$error)
        }
        if (is.null(case$warning)) {
            expect_length(then$warned, 0)
        } else {
            expect_match(then$warned, case$warning, fixed = TRUE)
        }
    })
}

## qc_dataset() runs the same rules without stopping, so its first FAIL row is
## the message check_dataset() stops with, without a prefix, and a dataset
## check_dataset() passes has none. Both are called with their defaults, so a
## case's own arguments aren't used
for (what in names(cases)) {
    test_that(paste0("qc_dataset() fails first where check_dataset() stops: ", what), {
        d <- cases[[what]]$d
        qc <- qc_dataset(d)
        expect_identical(qc$detail[qc$status == "FAIL"][1], outcome(check_dataset(d))$error)
    })
}

## ---- a mistake 6.0.3 didn't check doesn't take over its messages -------------
## The rules added since 6.0.3 are placed so that a dataset 6.0.3 refused with
## its own message keeps it: a rule for input 6.0.3 failed on with R's error
## goes where that error was, and a rule for input 6.0.3 took goes after
## 6.0.3's rules. Each mistake below is added to each case's dataset that holds
## the element, and wherever 6.0.3 still stops at the case's own message,
## check_dataset() must too. Each case is checked as coloc.abf() checks it, and
## as SuSiE does, with LD required, unless the case sets req itself.
mistakes <- list(
    "a character beta" = function(d)
        if (!is.null(d$beta)) { d$beta <- as.character(d$beta); d },
    "a character varbeta" = function(d)
        if (!is.null(d$varbeta)) { d$varbeta <- as.character(d$varbeta); d },
    "character p-values" = function(d)
        if (!is.null(d$pvalues)) { d$pvalues <- as.character(d$pvalues); d },
    "a negative varbeta" = function(d)
        if (!is.null(d$varbeta)) { d$varbeta[2] <- -abs(d$varbeta[2]); d },
    "two types" = function(d) { d$type <- rep(d$type, 2); d },
    "two values of s" = function(d) { d$s <- c(0.3, 0.3); d },
    "a vector LD" = function(d)
        if (!is.null(d$LD)) { d$LD <- as.vector(d$LD); d },
    "an N of 0" = function(d) { d$N <- 0; d },
    "an infinite N" = function(d) { d$N <- Inf; d },
    "an sdY of 0" = function(d) { d$sdY <- 0; d },
    "a z with a missing value" = function(d)
        if (!is.null(d$snp)) { d$z <- replace(rep(1, length(d$snp)), 1, NA); d },
    "a short z" = function(d)
        if (!is.null(d$snp)) { d$z <- rep(1, length(d$snp) - 1); d },
    "a character z" = function(d)
        if (!is.null(d$snp)) { d$z <- rep("1", length(d$snp)); d },
    "an infinite z" = function(d)
        if (!is.null(d$snp)) { d$z <- replace(rep(1, length(d$snp)), 1, Inf); d })

for (what in names(mistakes)) {
    test_that(paste0("a dataset 6.0.3 refused keeps its message with ", what, " added"), {
        for (use in c("single", "susie")) for (case in cases) {
            if (is.null(case$error) || !is.list(case$d))
                next
            d <- mistakes[[what]](case$d)
            if (is.null(d))
                next
            args <- c(list(d, 1), case$args)
            if (is.null(args$req))
                args$req <- .dataset_uses[[use]]
            then <- outcome(get_match_6_0_3(do.call(check_dataset, args)))$error
            if (identical(then, case$error))
                expect_identical(outcome(do.call(check_dataset, args))$error, today(case, then, 1))
        }
    })
}

## ---- what else check_dataset() and check_ld() must keep -----------------------
## The rules run in .gate() (R/check.R), but an error or warning shows the call
## the user made, as 6.0.3's did (call_of() is in helper-datasets.R)

test_that("check_dataset() errors and warns with its own call, as 6.0.3's did", {
    dup <- D1
    dup$snp[2] <- dup$snp[1]
    weak <- D1
    weak$beta <- weak$beta / 100
    per_snp_N <- D1
    per_snp_N$N <- rep(D1$N, length(D1$snp))
    expect_identical(call_of(check_dataset(dup, 1)), quote(check_dataset(dup, 1)))
    expect_identical(call_of(check_dataset(dup, 1)),
                     call_of(get_match_6_0_3(check_dataset(dup, 1))))
    expect_identical(call_of(check_dataset(weak, 1)),
                     call_of(get_match_6_0_3(check_dataset(weak, 1))))
    expect_identical(call_of(check_dataset(1:3, 1)),
                     call_of(get_match_6_0_3(check_dataset(1:3, 1))))
    expect_identical(call_of(check_dataset(per_snp_N, 1)),
                     call_of(get_match_6_0_3(check_dataset(per_snp_N, 1))))
    ## and through a function that calls it
    expect_identical(call_of(suppressMessages(coloc.abf(dup, D1))),
                     quote(check_dataset(d = dataset1, 1)))
})

test_that("check_ld(), today's against 6.0.3's", {
    names_differ <- D3$LD
    rownames(names_differ) <- paste0("r", seq_len(nrow(names_differ)))
    lds <- list(valid = D3$LD, missing = NULL, data_frame = as.data.frame(D3$LD),
                not_square = D3$LD[, -1], names_differ = names_differ,
                snp_missing = D3$LD[-1, -1])
    for (nm in names(lds)) {
        LD <- lds[[nm]]
        expect_identical(outcome(check_ld(D3, LD)),
                         outcome(get_match_6_0_3(check_ld(D3, LD))), label = nm)
        expect_identical(call_of(check_ld(D3, LD)),
                         call_of(get_match_6_0_3(check_ld(D3, LD))), label = nm)
    }
})
