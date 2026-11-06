## The plain-list note: coloc() and wakefield_bf() convert a plain list with
## as_coloc_dataset() and say so once per session, one note for both.
## setup-notes.R marks it as given for the other tests; each test here runs
## inside forgetting_notes(), as a session that hasn't had it. Each call's
## warnings are worked out before they are checked, as some expectations
## evaluate their object twice, and the note comes only once.

L1 <- coloc_test_data$D1
L2 <- coloc_test_data$D2
D1 <- as_coloc_dataset(L1)
D2 <- as_coloc_dataset(L2)

## each warning an expression gives, by its own class, the first of its classes
warning_classes <- function(expr) {
    got <- character(0)
    withCallingHandlers(expr, warning = function(w) {
        got <<- c(got, class(w)[1])
        invokeRestart("muffleWarning")
    })
    got
}

test_that("coloc() notes a plain list once per session", {
    forgetting_notes({
        first <- warning_classes(coloc(L1, D2))
        second <- warning_classes(coloc(L1, D2))
        both <- warning_classes(coloc(L1, L2))
    })
    expect_identical(first, "coloc_plain_list")
    expect_identical(second, character(0))
    expect_identical(both, character(0))
})

test_that("one note comes whichever trait, or both, is a plain list", {
    both <- forgetting_notes(warning_classes(coloc(L1, L2)))
    second <- forgetting_notes(warning_classes(coloc(D1, L2)))
    expect_identical(both, "coloc_plain_list")
    expect_identical(second, "coloc_plain_list")
})

test_that("coloc_datasets get no note", {
    expect_no_warning(forgetting_notes(coloc(D1, D2)))
})

test_that("wakefield_bf() notes a plain list, and shares the note with coloc()", {
    forgetting_notes({
        from_bf <- warning_classes(wakefield_bf(L1))
        then_coloc <- warning_classes(coloc(L1, D2))
    })
    expect_identical(from_bf, "coloc_plain_list")
    expect_identical(then_coloc, character(0))
    forgetting_notes({
        from_coloc <- warning_classes(coloc(L1, D2))
        then_bf <- warning_classes(wakefield_bf(L1))
    })
    expect_identical(from_coloc, "coloc_plain_list")
    expect_identical(then_bf, character(0))
})

test_that("the older functions take plain lists without a note, and leave it to come", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    L3 <- coloc_test_data$D3
    L4 <- coloc_test_data$D4
    fit4 <- suppressMessages(runsusie(L4))
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
    forgetting_notes({
        older <- list(coloc.susie = warning_classes(suppressMessages(coloc.susie(L3, L4))),
                      coloc.susie_bf = warning_classes(suppressMessages(coloc.susie_bf(L3, bf4))),
                      coloc.abf = warning_classes(utils::capture.output(coloc.abf(L1, L2))),
                      runsusie = warning_classes(suppressMessages(runsusie(L3))),
                      check_dataset = warning_classes(check_dataset(L1)))
        then_coloc <- warning_classes(coloc(L1, D2))
    })
    for (fn in names(older))
        expect_identical(older[[fn]], character(0), label = fn)
    expect_identical(then_coloc, "coloc_plain_list")
})

test_that("a broken plain list gets only its error, and the note is still to come", {
    bad <- L1
    bad$beta[1] <- NA
    forgetting_notes({
        expect_no_warning(expect_error(coloc(bad, D2),
                                       "dataset 1: beta contains missing values", fixed = TRUE))
        then_good <- warning_classes(coloc(L1, D2))
    })
    expect_identical(then_good, "coloc_plain_list")
    forgetting_notes({
        expect_no_warning(e <- expect_error(wakefield_bf(bad)))
        then_good <- warning_classes(wakefield_bf(L1))
    })
    expect_identical(conditionMessage(e), "beta contains missing values")
    expect_identical(then_good, "coloc_plain_list")
})

test_that("the note counts as given when its warning stops the call, as under options(warn = 2)", {
    ## options(warn = 2) can't be tried here, as testthat takes warnings
    ## itself; stopping at the warning leaves warning() the same way
    forgetting_notes({
        stopped <- tryCatch(coloc(L1, D2), warning = identity)
        again <- warning_classes(coloc(L1, D2))
    })
    expect_s3_class(stopped, "coloc_plain_list")
    expect_identical(again, character(0))
})

test_that("coloc.susie() and coloc.susie_bf() warn once per weak dataset, and leave the note to the user's own coloc()", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    weak <- function(d) {
        d$beta <- d$beta / 100
        d
    }
    w3 <- weak(coloc_test_data$D3)
    w4 <- weak(coloc_test_data$D4)
    fit4 <- suppressMessages(runsusie(coloc_test_data$D4))
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
    forgetting_notes({
        got <- warning_classes(suppressMessages(coloc.susie(w3, w4)))
        expect_equal(sum(got == "coloc_min_p"), 2)
        expect_false("coloc_plain_list" %in% got)
        got <- warning_classes(suppressMessages(coloc.susie_bf(w3, bf4)))
        expect_equal(sum(got == "coloc_min_p"), 1)
        expect_false(any(c("coloc_plain_list", "coloc_method_ignored") %in% got))
        ## the old functions didn't use up the note
        expect_true("coloc_plain_list" %in% warning_classes(coloc(L1, L2)))
    })
})
