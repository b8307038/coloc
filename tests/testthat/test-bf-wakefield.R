## wakefield_bf() against 6.0.3's process.dataset(), copied into helper-v603.R

data(coloc_test_data)
D1 <- coloc_test_data$D1

## D1 on each route wakefield_bf() takes
as_cc <- function(d) {
    d$type <- "cc"
    d$s <- 0.3
    d$sdY <- NULL
    d
}
to_pvalues <- function(d) {
    d$pvalues <- 2 * pnorm(-abs(d$beta / sqrt(d$varbeta)))
    d[c("beta", "varbeta", "sdY")] <- NULL
    d
}
routes <- list(
    "beta and varbeta, sdY given" = D1,
    "beta and varbeta, sdY estimated" = local({ d <- D1; d$sdY <- NULL; d }),
    "beta and varbeta, case-control" = as_cc(D1),
    "p-values" = to_pvalues(D1),
    "p-values, case-control" = to_pvalues(as_cc(D1)),
    "beta and varbeta, no position" = local({ d <- D1; d$position <- NULL; d }),
    "p-values, no position" = local({ d <- to_pvalues(D1); d$position <- NULL; d }))

for (what in names(routes)) {
    test_that(paste0("wakefield_bf() gives 6.0.3's numbers: ", what), {
        d <- routes[[what]]
        got <- suppressWarnings(wakefield_bf(d))
        expected <- suppressWarnings(get_match_6_0_3(process.dataset(d, suffix = "df1")))
        expect_identical(got$snp, expected$snp)
        for (v in c("z", "V", "r", "lABF"))
            expect_identical(got[[v]], expected[[paste0(v, ".df1")]], label = v)
        if ("position" %in% names(d)) {
            expect_identical(got$position, expected$position)
        } else {
            expect_identical(got$position, rep(NA_integer_, nrow(got)))
        }
    })
}

test_that("wakefield_bf() warns when it estimates sdY, as 6.0.3 did", {
    d <- routes[["beta and varbeta, sdY estimated"]]
    msg <- "estimating sdY from maf and varbeta"
    expect_warning(wakefield_bf(d), msg, fixed = TRUE)
    expect_warning(get_match_6_0_3(process.dataset(d, suffix = "df1")), msg, fixed = TRUE)
})

test_that("wakefield_bf() returns one row per snp, with no suffixes", {
    got <- wakefield_bf(D1)
    expect_named(got, c("snp", "position", "z", "V", "r", "lABF"))
    expect_identical(got$snp, D1$snp)
    expect_identical(got$position, D1$position)
})

test_that("effect_priors reaches both routes", {
    ## 6.0.3's approx.bf.estimates() takes effect_priors, though
    ## process.dataset() never passed it on (AL06) and the p-value route had
    ## 0.15 and 0.2 written in. Its formula is the p-value route's when sdY is 1
    wide <- c(quant = 0.3, cc = 0.4)
    for (what in c("beta and varbeta, sdY given", "beta and varbeta, case-control",
                   "p-values", "p-values, case-control")) {
        d <- routes[[what]]
        default <- wakefield_bf(d)
        got <- wakefield_bf(d, effect_priors = wide)
        sdY <- if (is.null(d$sdY)) 1 else d$sdY
        expected <- get_match_6_0_3(approx.bf.estimates(z = default$z, V = default$V,
                                                        type = d$type, sdY = sdY,
                                                        effect_priors = wide))
        expect_identical(got$lABF, expected$lABF, label = what)
        expect_false(identical(got$lABF, default$lABF), label = what)
        ## the prior for the other type of trait changes nothing
        other <- if (d$type == "quant") c(quant = 0.15, cc = 0.4) else c(quant = 0.3, cc = 0.2)
        expect_identical(wakefield_bf(d, effect_priors = other), default, label = what)
    }
})

test_that("wakefield_bf() checks its dataset, and takes a coloc_dataset", {
    d <- D1
    d$beta[1] <- NA
    expect_error(wakefield_bf(d), "beta contains missing values", fixed = TRUE)
    expect_identical(wakefield_bf(as_coloc_dataset(D1)), wakefield_bf(D1))
})

test_that("wakefield_bf() refuses effect_priors without a positive number for each type", {
    msg <- "effect_priors should be positive numbers named quant and cc"
    expect_error(wakefield_bf(D1, effect_priors = c(quant = 0.15)), msg, fixed = TRUE)
    expect_error(wakefield_bf(D1, effect_priors = c(0.15, 0.2)), msg, fixed = TRUE)
    expect_error(wakefield_bf(D1, effect_priors = c(quant = 0, cc = 0.2)), msg, fixed = TRUE)
    expect_error(wakefield_bf(D1, effect_priors = c(quant = "0.15", cc = "0.2")), msg, fixed = TRUE)
})
