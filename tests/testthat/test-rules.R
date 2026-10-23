## .dataset_rules, the list of rules check_dataset() runs in order (R/check.R).
## What each rule does is tested in test-check.R, against 6.0.3; these check
## the list itself.

test_that("each rule has a unique id, which names it in the list", {
    ids <- vapply(.dataset_rules, `[[`, "", "id", USE.NAMES = FALSE)
    expect_false(anyDuplicated(ids) > 0)
    expect_identical(names(.dataset_rules), ids)
})

test_that("each rule's elements are elements of a dataset", {
    for (rule in .dataset_rules)
        expect_true(all(rule$about %in% .dataset_items), label = rule$id)
})

test_that("the minimum p-value step is the one warning", {
    levels <- vapply(.dataset_rules, `[[`, "", "level")
    expect_true(all(levels %in% c("error", "warning")))
    expect_identical(names(levels)[levels == "warning"], "min_p")
})

test_that("the rules only qc_dataset() runs have ids of their own, and their elements are a dataset's", {
    ids <- vapply(.report_rules, `[[`, "", "id", USE.NAMES = FALSE)
    expect_identical(names(.report_rules), ids)
    expect_false(anyDuplicated(c(names(.dataset_rules), ids)) > 0)
    for (rule in .report_rules) {
        expect_identical(rule$level, "report", label = rule$id)
        expect_true(all(rule$about %in% .dataset_items), label = rule$id)
    }
})

test_that("qc_pair()'s rules have unique ids, and their elements are a dataset's", {
    ids <- vapply(.pair_rules, `[[`, "", "id", USE.NAMES = FALSE)
    expect_identical(names(.pair_rules), ids)
    expect_false(anyDuplicated(ids) > 0)
    for (rule in .pair_rules)
        expect_true(all(rule$about %in% .dataset_items), label = rule$id)
})

test_that("check_ld() runs the LD rules, in order", {
    expect_identical(names(.rules_about("LD")),
                     c("LD_dims", "LD_square", "LD_dimnames", "LD_covers", "LD_numeric",
                       "LD_correlation"))
})
