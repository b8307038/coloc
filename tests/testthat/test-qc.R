## qc_dataset(): every rule check_dataset() runs, then four of its own, listed
## without stopping. That its first failure is check_dataset()'s message is
## tested in test-check.R, over that file's cases. Then qc_pair(), two
## datasets against each other.

D1 <- plain_test_data$D1
D2 <- plain_test_data$D2
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
    ## where they pass, the rows say what they found
    qc <- qc_dataset(D1)
    expect_identical(detail_of(qc)[c("min_p", "sdY_plausible")],
                     c(min_p = "minimum p value is < 2.22e-16", sdY_plausible = "sdY 1.1, estimated 1.29"))
})

test_that("an sdY that isn't a single number leaves its plausibility not checked", {
    for (sdY in list(c(1.1, 1.1), numeric(0), rep(1.1, 500), "1.1")) {
        d <- D1
        d$sdY <- sdY
        got <- outcome(qc_dataset(d))
        expect_identical(got$error, NA_character_)
        expect_length(got$warned, 0)
        expect_identical(status_of(got$value)[c("sdY_valid", "sdY_plausible")],
                         c(sdY_valid = "FAIL", sdY_plausible = "not checked"))
        expect_identical(detail_of(got$value)[["sdY_plausible"]], "sdY isn't a single number")
    }
})

test_that("a check that gives anything but one sentence is not checked", {
    for (said in list(c("one", "two"), character(0), NA_character_, 1)) {
        rule <- .rule("odd", level = "report", function(d, ctx) said)
        expect_identical(.report_row(rule, list(), list(nd = character(0))),
                         c("not checked", "could not be checked: no single message"))
    }
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

test_that("a z that differs from beta / sqrt(varbeta) is caught, but not rounding", {
    made <- as_coloc_dataset(D3)
    made$z <- made$beta / sqrt(made$varbeta)
    expect_identical(status_of(qc_dataset(made))[["z_matches"]], "ok")
    ## beta and varbeta rounded to 3 significant figures, as summary
    ## statistics often are, against 2
    rounded <- made
    rounded$beta <- signif(made$beta, 3)
    rounded$varbeta <- signif(made$varbeta, 3)
    expect_identical(status_of(qc_dataset(rounded))[["z_matches"]], "ok")
    rounded$beta <- signif(made$beta, 2)
    rounded$varbeta <- signif(made$varbeta, 2)
    expect_identical(status_of(qc_dataset(rounded))[["z_matches"]], "warn")
    ## a z for the other allele at one snp
    made$z[7] <- -made$z[7]
    qc <- qc_dataset(made)
    expect_identical(to_look_at(qc), "z_matches")
    expect_match(detail_of(qc)[["z_matches"]],
                 "z differs from beta / sqrt(varbeta) by more than 1% at 1 snp, by up to", fixed = TRUE)
})

test_that("z_valid compares z alone with snp, as lengths compares the rest", {
    d <- D3
    d$z <- d$beta / sqrt(d$varbeta)
    d$beta <- d$beta[-1]
    qc <- qc_dataset(d, engine = "susie")
    expect_identical(qc$check[qc$status == "FAIL"], "lengths")
    d$z <- d$z[-1]
    expect_identical(detail_of(qc_dataset(d, engine = "susie"))[["z_valid"]],
                     "snp has 500 values but z has 499")
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

test_that("the engine decides which checks run, and what is required", {
    ## SuSiE's checks on LD and z, only for SuSiE, as check_dataset() runs them
    asymmetric <- as_coloc_dataset(D3)
    asymmetric$LD[1, 2] <- asymmetric$LD[1, 2] + 0.1
    expect_identical(status_of(qc_dataset(asymmetric))[["LD_correlation"]], "ok")
    expect_null(check_dataset(asymmetric))
    expect_identical(status_of(qc_dataset(asymmetric, engine = "susie"))[["LD_correlation"]], "FAIL")
    expect_error(check_dataset(asymmetric, req = .dataset_uses$susie), "LD must be a correlation matrix")
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

test_that("a part of a coloc_qc is a plain data frame", {
    qc <- qc_dataset(D1)
    expect_s3_class(qc, "coloc_qc")
    fails <- qc[qc$status == "FAIL", ]
    expect_identical(class(fails), "data.frame")
    expect_identical(class(qc[1:3, c("check", "status")]), "data.frame")
    expect_identical(qc[2, "check"], "no_na")
})

test_that("something other than a list gets one row", {
    qc <- qc_dataset(1:3)
    expect_s3_class(qc, "coloc_qc")
    expect_identical(c(qc$check, qc$status, qc$detail), c("list", "FAIL", "the dataset is not a list"))
})

test_that("a coloc_qc prints a count, then only the rows to look at, one line each", {
    ## R's default width, at which a data frame's print splits a long column
    old <- options(width = 80)
    on.exit(options(old))
    se <- D1
    se$varbeta <- sqrt(D1$varbeta)
    qc <- qc_dataset(se)
    out <- capture.output(expect_invisible(print(qc)))
    expect_identical(out[1], "26 checks ok, 4 not given; 4 to look at:")
    expect_length(out, 6)       # the count, the column names and four rows
    expect_identical(out[3], " min_p         warn        minimum p value is: 0.028794")
    expect_length(capture.output(print(qc, all = TRUE)), 2 + nrow(qc))
    expect_identical(capture.output(print(qc_dataset(as_coloc_dataset(D1)))),
                     "30 checks ok, 4 not given; nothing to look at")
})

## ---- qc_pair() -----------------------------------------------------------------

## D2 with its positions moved, as on another genome build, and 20 MAFs given
## for the other allele
moved <- D2
moved$position <- moved$position + 1000L
moved$MAF[1:20] <- 1 - moved$MAF[1:20]

test_that("qc_pair() finds nothing to look at in two packaged datasets, in either order", {
    qc <- qc_pair(D1, D2)
    expect_s3_class(qc, "coloc_qc")
    expect_identical(qc$check, names(.pair_rules))
    expect_identical(to_look_at(qc), character(0))
    expect_identical(detail_of(qc)[["snps_shared"]], "500 of 500 and 500")
    ## matched by snp id
    reversed <- D2
    for (v in c("beta", "varbeta", "MAF", "snp", "position"))
        reversed[[v]] <- rev(reversed[[v]])
    expect_identical(to_look_at(qc_pair(D1, reversed)), character(0))
})

test_that("qc_pair() warns about positions and MAFs that differ for the same snp", {
    qc <- qc_pair(D1, moved)
    expect_identical(to_look_at(qc), c("positions_agree", "maf_agree"))
    expect_identical(detail_of(qc)[["positions_agree"]],
                     "500 shared snps differ in position: another genome build?")
    expect_match(detail_of(qc)[["maf_agree"]], "20 shared snps differ in MAF by more than 0.1", fixed = TRUE)
    expect_identical(status_of(qc_pair(D1, moved, maf_diff = 0.99))[["maf_agree"]], "ok")
    one <- D2
    one$position[1] <- one$position[1] + 1L
    expect_identical(detail_of(qc_pair(D1, one))[["positions_agree"]],
                     "1 shared snp differs in position: another genome build?")
})

test_that("qc_pair() fails when no snp is shared, and leaves out what either dataset lacks", {
    apart <- D2
    apart$snp <- paste0(apart$snp, "x")
    qc <- qc_pair(D1, apart)
    expect_identical(status_of(qc), c(snps_shared = "FAIL", positions_agree = "not checked",
                                      maf_agree = "not checked"))
    part <- D2
    for (v in c("beta", "varbeta", "MAF", "snp", "position"))
        part[[v]] <- part[[v]][1:100]
    expect_identical(detail_of(qc_pair(D1, part))[["snps_shared"]], "100 of 500 and 100")
    expect_identical(status_of(qc_pair(D1, no_position(D2)))[["positions_agree"]], "not given")
    expect_identical(status_of(qc_pair(no_position(D2), D1))[["positions_agree"]], "not given")
})

test_that("qc_pair() takes only lists", {
    qc <- qc_pair(D1, 1:3)
    expect_identical(c(qc$check, qc$status, qc$detail), c("lists", "FAIL", "each dataset must be a list"))
})

test_that("a coloc_qc with one check ok says so in the singular", {
    expect_identical(capture.output(print(qc_pair(D1, moved)))[1],
                     "1 check ok, 0 not given; 2 to look at:")
})
