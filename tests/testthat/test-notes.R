## The plain-list note: coloc() and wakefield_bf() convert a plain list with
## as_coloc_dataset() and say so once per session, one note for both.
## setup-notes.R marks it as given for the other tests; each test here runs
## inside forgetting_notes(), as a session that hasn't had it. Each call's
## warnings are worked out before they are checked, as some expectations
## evaluate their object twice, and the note comes only once.

L1 <- plain_test_data$D1
L2 <- plain_test_data$D2
D1 <- as_coloc_dataset(L1)
D2 <- as_coloc_dataset(L2)

## the warnings an expression gives
warned <- function(expr) outcome(expr)$warned

test_that("coloc() notes a plain list once per session, naming its trait", {
    forgetting_notes({
        first <- warned(coloc(L1, D2))
        second <- warned(coloc(L1, D2))
        both <- warned(coloc(L1, L2))
    })
    expect_identical(first, paste(
        "coloc() was given a plain list for trait 1 and converted it with as_coloc_dataset().",
        "Make datasets with coloc_dataset() or as_coloc_dataset().",
        "(shown once per session)", sep = "\n"))
    expect_identical(second, character(0))
    expect_identical(both, character(0))
})

test_that("the note names every trait given as a plain list", {
    both <- forgetting_notes(warned(coloc(L1, L2)))
    second <- forgetting_notes(warned(coloc(D1, L2)))
    expect_match(both, "coloc() was given plain lists for traits 1 and 2 and converted them", fixed = TRUE)
    expect_match(second, "coloc() was given a plain list for trait 2 and converted it", fixed = TRUE)
})

test_that("coloc_datasets get no note", {
    expect_identical(forgetting_notes(warned(coloc(D1, D2))), character(0))
})

test_that("wakefield_bf() notes a plain list, and shares the note with coloc()", {
    forgetting_notes({
        from_bf <- warned(wakefield_bf(L1))
        then_coloc <- warned(coloc(L1, D2))
    })
    expect_identical(from_bf, paste(
        "wakefield_bf() was given a plain list and converted it with as_coloc_dataset().",
        "Make datasets with coloc_dataset() or as_coloc_dataset().",
        "(shown once per session)", sep = "\n"))
    expect_identical(then_coloc, character(0))
    forgetting_notes({
        from_coloc <- warned(coloc(L1, D2))
        then_bf <- warned(wakefield_bf(L1))
    })
    expect_length(from_coloc, 1)
    expect_identical(then_bf, character(0))
})

test_that("the older functions take plain lists without a note, and leave it to come", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    L3 <- plain_test_data$D3
    L4 <- plain_test_data$D4
    fit4 <- suppressMessages(runsusie(L4))
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
    forgetting_notes({
        older <- list(coloc.susie = warned(suppressMessages(coloc.susie(L3, L4))),
                      coloc.susie_bf = warned(suppressMessages(coloc.susie_bf(L3, bf4))),
                      coloc.abf = warned(utils::capture.output(coloc.abf(L1, L2))),
                      runsusie = warned(suppressMessages(runsusie(L3))),
                      check_dataset = warned(check_dataset(L1)))
        then_coloc <- warned(coloc(L1, D2))
    })
    for (fn in names(older))
        expect_identical(older[[fn]], character(0), label = fn)
    expect_length(then_coloc, 1)
})

test_that("a broken plain list gets only its error, and the note is still to come", {
    bad <- L1
    bad$beta[1] <- NA
    forgetting_notes({
        broken <- outcome(coloc(bad, D2))
        then_good <- warned(coloc(L1, D2))
    })
    expect_identical(broken$error, "trait 1: beta contains missing values")
    expect_length(broken$warned, 0)
    expect_length(then_good, 1)
    forgetting_notes({
        broken <- outcome(wakefield_bf(bad))
        then_good <- warned(wakefield_bf(L1))
    })
    expect_identical(broken$error, "beta contains missing values")
    expect_length(broken$warned, 0)
    expect_length(then_good, 1)
})

test_that("the note counts as given when its warning stops the call, as under options(warn = 2)", {
    ## options(warn = 2) can't be tried here, as testthat takes warnings
    ## itself; stopping at the warning leaves warning() the same way
    forgetting_notes({
        stopped <- tryCatch(coloc(L1, D2), warning = conditionMessage)
        again <- warned(coloc(L1, D2))
    })
    expect_match(stopped, "coloc() was given a plain list", fixed = TRUE)
    expect_identical(again, character(0))
})
