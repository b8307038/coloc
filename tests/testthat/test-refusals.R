## Values coloc can't use, refused at each way a dataset is checked, and the
## messages reworded since 6.0.3, against 6.0.3's own copies:
## get_match_6_0_3(), in helper-v603.R, runs a call as 6.0.3 made it.

D1 <- plain_test_data$D1
D2 <- plain_test_data$D2
D3 <- plain_test_data$D3
D4 <- plain_test_data$D4

fit3 <- suppressMessages(runsusie(D3))

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
    ## the reports take whatever a door refuses, never stopping or warning
    testthat::expect_no_condition(qc_dataset(d, engine = if (susie) "susie" else "wakefield"))
    testthat::expect_no_condition(qc_pair(other, d))
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
    ## refused after 6.0.3's LD rules, so one 6.0.3 refused keeps its message
    not_square <- D1
    not_square$LD <- array(D1$LD[, -1], c(500, 499, 1), dimnames = list(D1$snp, D1$snp[-1], NULL))
    expect_identical(outcome(get_match_6_0_3(check_dataset(not_square, 2)))$error, "LD not square")
    expect_identical(outcome(check_dataset(not_square, 2))$error, "dataset 2: LD not square")
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
    ## and above 1, as susieR needs, which 6.0.3 also left to susieR. On the
    ## Wakefield route an N of 1 gives 6.0.3's answer
    for (n in c(1, 0.5)) {
        small4 <- D4
        small4$N <- n
        expect_error(suppressMessages(get_match_6_0_3(runsusie(small4))), "n must be greater than 1")
        expect_identical(error_before_susie(coloc(D3, small4, method = "susie")),
                         "trait 2: N must be above 1 for SuSiE")
    }
    one <- D1
    one$N <- 1
    expect_identical(at_each_door(one), everywhere(NA))
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
