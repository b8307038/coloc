## Values coloc can't use, refused at each way a dataset is checked, against
## 6.0.3's own copies: get_match_6_0_3(), in helper-v603.R, runs a call as
## 6.0.3 made it.

D1 <- coloc_test_data$D1
D2 <- coloc_test_data$D2
D3 <- coloc_test_data$D3
D4 <- coloc_test_data$D4

fit3 <- suppressMessages(runsusie(D3))

## ---- values coloc can't use, at each way a dataset is checked ---------------
## check_dataset() checks a dataset for coloc.abf(), finemap.abf(), the
## constructors, wakefield_bf(), runsusie() and coloc(). at_each_door() checks
## that each of them refuses d with the error given, or with error = NULL that
## each goes on. runsusie() is asked only about a dataset it could fine-map.
## 6.0.3's copies raise R's own errors, so a test reads only whether they
## stopped.

at_each_door <- function(d, error, other = D1, susie = FALSE) {
    ## each door, and the dataset number its error begins with: 2 where d is the
    ## second of two datasets, 1 for runsusie(), none where d is checked alone
    doors <- list(list(call = quote(check_dataset(d, 2)), dataset = 2),
                  list(call = quote(coloc.abf(other, d)), dataset = 2),
                  list(call = quote(finemap.abf(d)), dataset = NULL),
                  list(call = quote(as_coloc_dataset(d)), dataset = NULL),
                  list(call = quote(coloc(other, d)), dataset = 2),
                  list(call = quote(wakefield_bf(d)), dataset = NULL))
    if (susie)
        doors <- c(doors, list(list(call = quote(runsusie(d)), dataset = 1)))
    quietly <- function(expr) utils::capture.output(suppressWarnings(suppressMessages(expr)))
    ## !! puts the door's call in the message of a failed expectation
    for (door in doors) {
        if (!is.null(error)) {
            expected <- if (is.null(door$dataset)) error
                        else paste0("dataset ", door$dataset, ": ", error)
            expect_error(quietly(!!door$call), expected, fixed = TRUE)
        } else {
            expect_no_error(quietly(!!door$call))
        }
    }
}

## coloc.abf()'s summary, without its printing
abf_summary <- function(d1, d2) {
    utils::capture.output(res <- coloc.abf(d1, d2))
    res$summary
}

## checks that a call stops with the error given, before SuSiE said anything,
## as it does when it starts fine-mapping
refused_before_susie <- function(expr, error) {
    susie_ran <- FALSE
    withCallingHandlers(
        expect_error(expr, error, fixed = TRUE),
        message = function(m) {
            susie_ran <<- TRUE
            invokeRestart("muffleMessage")
        })
    expect_false(susie_ran)
}

test_that("beta, varbeta and p-values must be numeric where they are read, where 6.0.3 failed with R's errors", {
    chr_beta <- D1
    chr_beta$beta <- as.character(chr_beta$beta)
    chr_varbeta <- D1
    chr_varbeta$varbeta <- as.character(chr_varbeta$varbeta)
    chr_p <- to_pvalues(D1)
    chr_p$pvalues <- as.character(chr_p$pvalues)
    at_each_door(chr_beta, error = "beta must be numeric", susie = TRUE)
    at_each_door(chr_varbeta, error = "varbeta must be numeric", susie = TRUE)
    ## 6.0.3's check let character p-values through, and coloc.abf() failed after
    expect_error(get_match_6_0_3(coloc.abf(D1, chr_p)))
    at_each_door(chr_p, error = "pvalues must be numeric")
    ## 6.0.3 computed with a logical beta
    lgl_beta <- D1
    lgl_beta$beta <- lgl_beta$beta > 0
    at_each_door(lgl_beta, error = "beta must be numeric")
    ## where nothing reads them they pass, with 6.0.3's answer
    p_beside_beta <- D1
    p_beside_beta$pvalues <- chr_p$pvalues
    beta_beside_p <- to_pvalues(D1)
    beta_beside_p$beta <- chr_beta$beta
    varbeta_beside_p <- to_pvalues(D1)
    varbeta_beside_p$varbeta <- chr_varbeta$varbeta
    for (d in list(p_beside_beta, beta_beside_p, varbeta_beside_p)) {
        at_each_door(d, error = NULL)
        expect_equal(abf_summary(D1, d), get_match_6_0_3(coloc.abf(D1, d))$summary)
    }
    ## with a zero varbeta as well
    zero_too <- chr_beta
    zero_too$varbeta[1] <- 0
    at_each_door(zero_too, error = "varbeta has zeros, such as at s1; a variance must be above 0")
})

test_that("a negative varbeta is refused, where 6.0.3 failed with R's error", {
    neg <- D1
    neg$varbeta[5] <- -neg$varbeta[5]
    at_each_door(neg, error = "varbeta must be positive", susie = TRUE)
    ## with a zero as well
    zero_too <- neg
    zero_too$varbeta[1] <- 0
    at_each_door(zero_too, error = "varbeta has zeros, such as at s1; a variance must be above 0")
    ## beside p-values varbeta isn't read, and passes as in 6.0.3
    beside_p <- to_pvalues(D1)
    beside_p$varbeta <- -D1$varbeta
    at_each_door(beside_p, error = NULL)
    ## a coloc_dataset changed to hold one gave NA posteriors, as its minimum
    ## p-value comparison is skipped
    made <- as_coloc_dataset(D1)
    made$varbeta[5] <- -made$varbeta[5]
    expect_error(coloc(made, D2), "dataset 1: varbeta must be positive", fixed = TRUE)
    expect_error(coloc.abf(made, D2), "dataset 1: varbeta must be positive", fixed = TRUE)
    expect_error(wakefield_bf(made), "varbeta must be positive", fixed = TRUE)
})

test_that("type and s must be single values, where 6.0.3 failed with R's errors", {
    for (type in list(c("quant", "cc"), rep("quant", length(D1$snp)), character(0))) {
        d <- D1
        d$type <- type
        at_each_door(d, error = "type must be a single value, quant or cc")
    }
    two_s <- as_cc(D1)
    two_s$s <- c(0.5, 0.5)
    at_each_door(two_s, error = "s must be a single number between 0 and 1")
})

test_that("LD must be a matrix, and for SuSiE a numeric one, where 6.0.3 failed with R's errors", {
    vector_ld <- D1
    vector_ld$LD <- as.vector(D1$LD)
    at_each_door(vector_ld, error = "LD must be a matrix", susie = TRUE)
    ## 6.0.3 took an LD with a third dimension
    three_d <- D1
    three_d$LD <- array(D1$LD, dim = c(dim(D1$LD), 1), dimnames = c(dimnames(D1$LD), list(NULL)))
    at_each_door(three_d, error = "LD must be a matrix", susie = TRUE)
    ## and one that isn't square as well
    not_square <- D1
    not_square$LD <- array(D1$LD[, -1], c(500, 499, 1), dimnames = list(D1$snp, D1$snp[-1], NULL))
    expect_error(check_dataset(not_square, 2),
                 "dataset 2: LD must be square, but it has 500 rows and 499 columns", fixed = TRUE)
    ## a data frame LD has two dimensions: where LD isn't required it passes,
    ## with 6.0.3's answer
    df_ld <- D1
    df_ld$LD <- as.data.frame(D1$LD)
    at_each_door(df_ld, error = NULL)
    expect_equal(abf_summary(D1, df_ld), get_match_6_0_3(coloc.abf(D1, df_ld))$summary)
    ## SuSiE needs a numeric matrix. 6.0.3's check_ld() refused a data frame
    ## after check_dataset(), and a character matrix failed in R; now the check
    ## refuses both, before either trait is fine-mapped
    df_ld3 <- D3
    df_ld3$LD <- as.data.frame(D3$LD)
    chr_ld3 <- D3
    storage.mode(chr_ld3$LD) <- "character"
    expect_error(get_match_6_0_3(runsusie(df_ld3)))
    expect_error(suppressMessages(get_match_6_0_3(runsusie(chr_ld3))))
    for (d in list(df_ld3, chr_ld3)) {
        expect_error(runsusie(d), "dataset 1: LD must be a numeric matrix", fixed = TRUE)
        refused_before_susie(coloc(D4, d, method = "susie"),
                             "dataset 2: LD must be a numeric matrix")
    }
})

test_that("N must be a finite number above 0 wherever it is given, where 6.0.3 took any", {
    ## beside beta and sdY nothing reads it, and 6.0.3 gave its answer
    for (n in c(0, -5, Inf)) {
        d <- D1
        d$N <- n
        at_each_door(d, error = "N must be a finite number above 0", susie = TRUE)
    }
    ## without sdY, which is estimated from N, 6.0.3 went on all the same
    no_sdY <- D1
    no_sdY$sdY <- NULL
    no_sdY$N <- 0
    expect_no_error(suppressWarnings(get_match_6_0_3(coloc.abf(D1, no_sdY))))
    at_each_door(no_sdY, error = "N must be a finite number above 0")
    ## given as p-values, 0 is refused, as in 6.0.3, and so is Inf, which
    ## 6.0.3 took
    p_zero <- to_pvalues(D1)
    p_zero$N <- 0
    at_each_door(p_zero, error = "N must be a finite number above 0")
    p_inf <- to_pvalues(D1)
    p_inf$N <- Inf
    at_each_door(p_inf, error = "N must be a finite number above 0")
    ## for SuSiE, refused before either trait is fine-mapped; 6.0.3 left it to
    ## susieR, after the other trait was
    zero4 <- D4
    zero4$N <- 0
    expect_error(suppressMessages(get_match_6_0_3(runsusie(zero4))))
    refused_before_susie(coloc(D3, zero4, method = "susie"), "dataset 2: N must be a finite number above 0")
    ## and above 1, as susieR needs, which 6.0.3 also left to susieR. On the
    ## Wakefield route an N of 1 gives 6.0.3's answer
    for (n in c(1, 0.5)) {
        small4 <- D4
        small4$N <- n
        expect_error(suppressMessages(get_match_6_0_3(runsusie(small4))))
        refused_before_susie(coloc(D3, small4, method = "susie"),
                             "dataset 2: N must be above 1 for SuSiE")
    }
    one <- D1
    one$N <- 1
    at_each_door(one, error = NULL)
})

test_that("sdY must be a single finite number above 0 where it is read, where 6.0.3 took any", {
    ## a quantitative trait given as beta and varbeta
    for (sdY in list(0, -1, c(1, 1), Inf, "1")) {
        d <- D1
        d$sdY <- sdY
        at_each_door(d, error = "sdY must be a single finite number above 0")
    }
    ## where nothing reads it, it passes with 6.0.3's answer: on a case-control
    ## trait, and with p-values
    cc <- as_cc(D1)
    cc$sdY <- 0
    p <- to_pvalues(D1)
    p$sdY <- 0
    for (d in list(cc, p)) {
        at_each_door(d, error = NULL)
        expect_equal(abf_summary(D1, d), get_match_6_0_3(coloc.abf(D1, d))$summary)
    }
    ## and for SuSiE, which never reads it
    susie3 <- D3
    susie3$sdY <- 0
    expect_null(check_dataset(susie3, 1, req = .dataset_uses$susie))
})

test_that("for SuSiE LD must be a correlation matrix, where 6.0.3 took any numeric matrix", {
    above_one <- D3
    above_one$LD[1, 2] <- 1.5
    above_one$LD[2, 1] <- 1.5
    covariance <- D3
    covariance$LD <- D3$LD * 0.5
    asymmetric <- D3
    asymmetric$LD[1, 2] <- asymmetric$LD[1, 2] + 0.1
    for (d in list(above_one, covariance, asymmetric)) {
        ## 6.0.3's runsusie() checked its LD with check_ld(), which took it
        expect_null(get_match_6_0_3(check_ld(d, d$LD)))
        refused_before_susie(runsusie(d),
                             "dataset 1: LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal")
        refused_before_susie(coloc(D4, d, method = "susie"),
                             "dataset 2: LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal")
        ## where LD isn't required it isn't tested, and the answer is 6.0.3's
        expect_null(check_dataset(d, 2))
        expect_equal(abf_summary(D4, d), get_match_6_0_3(coloc.abf(D4, d))$summary)
    }
    ## to within 1e-6, on the diagonal and off it
    susie_check <- function(d) check_dataset(d, 1, req = .dataset_uses$susie)
    for (off in c(5e-7, 2e-6)) {
        diagonal <- D3
        diag(diagonal$LD) <- 1 + off
        pair <- D3
        pair$LD[1, 2] <- pair$LD[1, 2] + off
        for (d in list(diagonal, pair)) {
            if (off < 1e-6)
                expect_null(susie_check(d))
            else
                expect_error(susie_check(d),
                             "dataset 1: LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal",
                             fixed = TRUE)
        }
    }
    ## only the dataset's snps, the part runsusie() uses: an LD covering more
    ## snps passes with a fault among the others, and not with one among them
    fewer <- D3
    for (v in c("beta", "varbeta", "snp", "position", "MAF"))
        fewer[[v]] <- fewer[[v]][-(1:10)]
    outside <- fewer
    outside$LD[1, 2] <- 1.5
    outside$LD[2, 1] <- 1.5
    inside <- fewer
    inside$LD[11, 12] <- 1.5
    inside$LD[12, 11] <- 1.5
    expect_null(susie_check(outside))
    expect_error(susie_check(inside),
                 "dataset 1: LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal",
                 fixed = TRUE)
    ## and a non-square LD is refused for SuSiE too
    not_square <- D3
    not_square$LD <- D3$LD[, -1]
    expect_error(susie_check(not_square),
                 "dataset 1: LD must be square, but it has 500 rows and 499 columns", fixed = TRUE)
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
    ## each z, and whether 6.0.3's runsusie() failed on it with R's or susieR's
    ## error or fine-mapped it as given
    broken <- list("a missing value" = list(z = replace(z3, 5, NA), failed = FALSE),
                   "one value short" = list(z = z3[-1], failed = TRUE),
                   "one value too many" = list(z = c(z3, 1), failed = TRUE),
                   "character" = list(z = as.character(z3), failed = TRUE),
                   "logical" = list(z = z3 > 0, failed = FALSE),
                   "an infinite value" = list(z = replace(z3, 5, Inf), failed = FALSE))
    ## and the error each gets now, after its "dataset 1: " label
    errors <- c("a missing value" = "z contains missing values",
                "one value short" = "snp has 500 values but z has 499",
                "one value too many" = "snp has 500 values but z has 501",
                "character" = "z must be numeric",
                "logical" = "z must be numeric",
                "an infinite value" = "z has infinite values, such as at s5")
    for (what in names(broken)) {
        b <- broken[[what]]
        d <- D3
        d$z <- b$z
        failed <- tryCatch({
            suppressWarnings(suppressMessages(get_match_6_0_3(runsusie(d))))
            FALSE
        }, error = function(e) TRUE)
        expect_identical(failed, b$failed, label = what)
        expect_error(runsusie(d), paste0("dataset 1: ", errors[[what]]), fixed = TRUE)
        refused_before_susie(coloc(D4, d, method = "susie"), paste0("dataset 2: ", errors[[what]]))
        ## where nothing reads z it passes, with 6.0.3's answer
        d1 <- D1
        d1$z <- b$z
        at_each_door(d1, error = NULL)
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

test_that("per-snp elements must have the same length", {
    short <- D1
    short$beta <- short$beta[-1]
    at_each_door(short, error = "snp has 500 values but beta has 499", susie = TRUE)
    ## each per-snp element is compared with snp, position too
    short_position <- D1
    short_position$position <- short_position$position[-1]
    expect_error(check_dataset(short_position),
                 "snp has 500 values but position has 499", fixed = TRUE)
    short_snp <- D1
    short_snp$snp <- short_snp$snp[-1]
    expect_error(check_dataset(short_snp), "snp has 499 values but beta has 500", fixed = TRUE)
    ## without snp, with the first element given
    no_snp <- short
    no_snp$snp <- NULL
    expect_error(check_dataset(no_snp, req = "type"),
                 "beta has 499 values but varbeta has 500", fixed = TRUE)
})

test_that("N must be a single number, not one per snp", {
    for (n in list(rep(D1$N, length(D1$snp)), as.character(D1$N))) {
        d <- D1
        d$N <- n
        at_each_door(d, susie = TRUE,
                     error = "N must be a single number; coloc does not take per-snp sample sizes")
    }
})

test_that("check_alignment() refuses a dataset without LD, beta or varbeta, or whose LD isn't a correlation matrix, where 6.0.3 failed with R's error or gave a share", {
    no_ld <- D3
    no_ld$LD <- NULL
    p_only <- D3[c("snp", "position", "MAF", "N", "type", "LD")]
    p_only$pvalues <- 2 * pnorm(-abs(D3$beta / sqrt(D3$varbeta)))
    for (d in list(no_ld, p_only))
        expect_error(get_match_6_0_3(check_alignment(d, do_plot = FALSE)))
    expect_error(check_alignment(no_ld, do_plot = FALSE),
                 "LD is missing, and coloc needs it here", fixed = TRUE)
    expect_error(check_alignment(p_only, do_plot = FALSE),
                 "beta and varbeta are missing, and coloc needs both here", fixed = TRUE)
    scaled_ld <- D3
    scaled_ld$LD <- 1.5 * D3$LD
    expect_no_error(get_match_6_0_3(check_alignment(scaled_ld, do_plot = FALSE)))
    expect_error(check_alignment(scaled_ld, do_plot = FALSE),
                 "LD must be a correlation matrix: values in [-1, 1], symmetric, with 1 on the diagonal",
                 fixed = TRUE)
    expect_equal(check_alignment(D3, do_plot = FALSE), 0.997, tolerance = 1e-3)
})
