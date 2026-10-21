## wakefield_bf() against 6.0.3's process.dataset(), copied into helper-v603.R

data(coloc_test_data)
D1 <- coloc_test_data$D1

## D1 on each route wakefield_bf() takes
routes <- list(
    "beta and varbeta, sdY given" = D1,
    "beta and varbeta, sdY estimated" = local({ d <- D1; d$sdY <- NULL; d }),
    "beta and varbeta, case-control" = as_cc(D1),
    "p-values" = to_pvalues(D1),
    "p-values, case-control" = to_pvalues(as_cc(D1)),
    "beta and varbeta, no position" = no_position(D1),
    "p-values, no position" = no_position(to_pvalues(D1)),
    "beta and varbeta, type a factor" = type_factor(D1, c("cc", "quant")),
    "p-values, case-control, type a factor" = type_factor(to_pvalues(as_cc(D1)), "cc"))

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
    ## the row names are row numbers, whatever the input's vectors are named
    named <- D1
    names(named$beta) <- names(named$position) <- D1$snp
    expect_identical(rownames(wakefield_bf(named)), as.character(seq_along(D1$snp)))
    ## no snps and no position: no rows, where a single NA position would fail
    none <- no_position(D1)
    for (v in c("beta", "varbeta", "MAF", "snp"))
        none[[v]] <- none[[v]][0]
    expect_identical(nrow(suppressWarnings(wakefield_bf(none))), 0L)
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

test_that("wakefield_bf() makes a list a coloc_dataset, and checks a coloc_dataset again", {
    ## a list goes through as_coloc_dataset(), which refuses a data.frame,
    ## where check_dataset() alone would report "missing required element(s) type"
    expect_error(wakefield_bf(as.data.frame(D1[c("snp", "beta", "varbeta")])),
                 "expected a coloc dataset", fixed = TRUE)
    ## a coloc_dataset changed since it was made is checked again
    changed <- as_coloc_dataset(D1)
    changed$beta[1] <- NA
    expect_error(wakefield_bf(changed), "beta contains missing values", fixed = TRUE)
    ## a weak list gets the minimum p-value warning once
    weak <- D1
    weak$beta <- weak$beta / 100
    warned <- character(0)
    withCallingHandlers(wakefield_bf(weak), warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    expect_length(grep("minimum p value", warned, fixed = TRUE), 1)
    ## for now an error shows the call of the check inside, as_coloc_dataset()'s
    ## for a list and check_dataset()'s for a coloc_dataset; pinned, so that a
    ## change to it is made on purpose
    bad <- D1
    bad$beta[1] <- NA
    expect_identical(tryCatch(wakefield_bf(bad), error = conditionCall),
                     quote(as_coloc_dataset(d)))
    expect_identical(tryCatch(wakefield_bf(changed), error = conditionCall),
                     quote(check_dataset(d)))
})

test_that("wakefield_bf() refuses effect_priors without a positive number for each type", {
    msg <- "effect_priors should be positive numbers named quant and cc"
    expect_error(wakefield_bf(D1, effect_priors = c(quant = 0.15)), msg, fixed = TRUE)
    expect_error(wakefield_bf(D1, effect_priors = c(0.15, 0.2)), msg, fixed = TRUE)
    expect_error(wakefield_bf(D1, effect_priors = c(quant = 0, cc = 0.2)), msg, fixed = TRUE)
    expect_error(wakefield_bf(D1, effect_priors = c(quant = "0.15", cc = "0.2")), msg, fixed = TRUE)
})
