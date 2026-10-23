## qc_dataset(): every rule check_dataset() runs, then four of its own, listed
## without stopping. That its first failure is check_dataset()'s message is
## tested in test-check.R, over that file's cases.

D1 <- plain_test_data$D1
D3 <- plain_test_data$D3

## a report's statuses and details, named by rule id
status_of <- function(qc) setNames(qc$status, qc$check)
detail_of <- function(qc) setNames(qc$detail, qc$check)
## the rows a print shows
to_look_at <- function(qc) qc$check[qc$status %in% c("FAIL", "warn", "not checked")]

test_that("qc_dataset() lists every rule check_dataset() runs, then its own, in order", {
    qc <- qc_dataset(D1)
    expect_s3_class(qc, "coloc_qc")
    expect_identical(names(qc), c("check", "status", "detail"))
    expect_identical(qc$check, c(names(.dataset_rules), names(.report_rules)))
})

test_that("a coloc_dataset made from each packaged dataset has nothing to look at, with either engine", {
    for (name in c("D1", "D2", "D3", "D4")) for (engine in c("wakefield", "susie")) {
        qc <- qc_dataset(as_coloc_dataset(plain_test_data[[name]]), engine = engine)
        expect_identical(to_look_at(qc), character(0), label = paste(name, engine))
    }
    ## as a plain list, the one row to look at says so
    expect_identical(to_look_at(qc_dataset(D1)), "coloc_dataset")
})

test_that("standard errors given as varbeta show in the minimum p-value and sdY rows", {
    se <- D1
    se$varbeta <- sqrt(D1$varbeta)
    ## a coloc_dataset too, which check_dataset() doesn't warn about
    made <- as_coloc_dataset(D1)
    made$varbeta <- sqrt(made$varbeta)
    for (d in list(se, made)) {
        qc <- qc_dataset(d)
        expect_identical(status_of(qc)[c("min_p", "sdY_plausible")],
                         c(min_p = "warn", sdY_plausible = "warn"))
        expect_identical(detail_of(qc)[["min_p"]], "minimum p value is: 0.028794")
        expect_match(detail_of(qc)[["sdY_plausible"]], "does varbeta hold variances?", fixed = TRUE)
    }
    expect_null(check_dataset(made))
})

test_that("z-scores whose signs disagree with LD are caught, and a dataset with no signal isn't judged", {
    made <- as_coloc_dataset(D1)
    expect_identical(detail_of(qc_dataset(made))[["z_ld_signs"]], "100% of 11628 pairs agree")
    set.seed(1)
    flip <- sample(500, 50)
    made$beta[flip] <- -made$beta[flip]
    qc <- qc_dataset(made)
    expect_identical(to_look_at(qc), "z_ld_signs")
    expect_match(detail_of(qc)[["z_ld_signs"]], "84% of 11628 pairs agree in sign", fixed = TRUE)
    ## no snp near |z| = 3, so no pair to judge, unless every pair counts
    weak <- as_coloc_dataset(D1)
    weak$beta <- weak$beta / 100
    expect_identical(status_of(qc_dataset(weak))[["z_ld_signs"]], "not checked")
    expect_match(detail_of(qc_dataset(weak))[["z_ld_signs"]], "too few pairs to judge: 0 with", fixed = TRUE)
    expect_identical(status_of(qc_dataset(weak, z_thr = 0))[["z_ld_signs"]], "ok")
})

test_that("a z that differs from beta / sqrt(varbeta) is caught", {
    made <- as_coloc_dataset(D3)
    made$z <- made$beta / sqrt(made$varbeta)
    expect_identical(status_of(qc_dataset(made))[["z_matches"]], "ok")
    made$z[7] <- made$z[7] + 0.5
    qc <- qc_dataset(made)
    expect_identical(to_look_at(qc), "z_matches")
    expect_identical(detail_of(qc)[["z_matches"]],
                     "z differs from beta / sqrt(varbeta) at 1 snp, by up to 0.5")
})

test_that("qc_dataset() never stops or warns, and gives each of several mistakes a FAIL row", {
    na_beta <- D1
    na_beta$beta[3] <- NA
    na_varbeta <- D1
    na_varbeta$varbeta[4] <- NA
    na_p <- to_pvalues(D1)
    na_p$pvalues[1] <- NA
    several <- D1
    several$type <- "quantitative"
    several$snp[2] <- several$snp[1]
    several$MAF[1] <- 1.2
    several$beta[3] <- NA
    for (d in list(na_beta, na_varbeta, na_p, several)) {
        got <- outcome(qc_dataset(d))
        expect_identical(got$error, NA_character_)
        expect_length(got$warned, 0)
    }
    qc <- qc_dataset(several)
    expect_identical(qc$check[qc$status == "FAIL"], c("no_na", "type_valid", "snp_unique", "maf_range"))
    ## R can't work out the p-values with an NA in beta
    expect_identical(status_of(qc)[["min_p"]], "not checked")
})

test_that("the checks SuSiE needs run whatever the engine, and the engine decides what is required", {
    asymmetric <- as_coloc_dataset(D3)
    asymmetric$LD[1, 2] <- asymmetric$LD[1, 2] + 0.1
    expect_identical(status_of(qc_dataset(asymmetric))[["LD_correlation"]], "FAIL")
    expect_null(check_dataset(asymmetric))
    no_N <- as_coloc_dataset(D3)
    no_N$N <- NULL
    expect_identical(status_of(qc_dataset(no_N))[["required"]], "ok")
    qc <- qc_dataset(no_N, engine = "susie")
    expect_identical(status_of(qc)[["required"]], "FAIL")
    expect_identical(detail_of(qc)[["required"]], "missing required element(s) N")
    ## SuSiE doesn't read sdY
    zero <- as_coloc_dataset(D3)
    zero$sdY <- 0
    expect_identical(status_of(qc_dataset(zero))[["sdY_valid"]], "FAIL")
    expect_identical(status_of(qc_dataset(zero, engine = "susie"))[["sdY_valid"]], "ok")
})

test_that("something other than a list gets one row", {
    qc <- qc_dataset(1:3)
    expect_s3_class(qc, "coloc_qc")
    expect_identical(c(qc$check, qc$status, qc$detail), c("list", "FAIL", "the dataset is not a list"))
})

test_that("a coloc_qc prints a count, then only the rows to look at", {
    old <- options(width = 200)
    on.exit(options(old))
    se <- D1
    se$varbeta <- sqrt(D1$varbeta)
    qc <- qc_dataset(se)
    out <- capture.output(expect_invisible(print(qc)))
    expect_identical(out[1], "25 checks ok, 4 not given; 4 to look at:")
    expect_length(out, 6)       # the count, the column names and four rows
    expect_length(capture.output(print(qc, all = TRUE)), 2 + nrow(qc))
    expect_identical(capture.output(print(qc_dataset(as_coloc_dataset(D1)))),
                     "29 checks ok, 4 not given; nothing to look at")
})
