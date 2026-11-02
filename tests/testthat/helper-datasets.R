## Shared by the test files: a dataset in other forms, and what a call gives.

## code run in a session that hasn't had the plain-list note, which
## setup-notes.R marks as given; the notes given are put back after
forgetting_notes <- function(code) {
    old <- .coloc_state$noted
    on.exit(.coloc_state$noted <- old)
    .coloc_state$noted <- character(0)
    code
}

## coloc_test_data's datasets as plain lists, which is what most tests are
## written for: a coloc_dataset skips the minimum p-value warning, and keeps
## its class when an element is changed
plain_test_data <- lapply(coloc_test_data, unclass)

## the same dataset as case-control, as p-values, and without position
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
no_position <- function(d) {
    d$position <- NULL
    d
}
## type as a factor: with levels cc and quant, as a data frame column holding
## both types gives, or with its own level alone
type_factor <- function(d, levels) {
    d$type <- factor(d$type, levels = levels)
    d
}

## what a call gives: its value, its error message (NA when there is none) and
## its warnings
outcome <- function(expr) {
    value <- NULL
    warned <- character(0)
    error <- tryCatch(withCallingHandlers({
        value <- expr
        NA_character_
    }, warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = conditionMessage)
    list(value = value, error = error, warned = warned)
}

## the call in a call's error or warning; NULL when it gives neither
call_of <- function(expr) tryCatch({
    expr
    NULL
}, error = conditionCall, warning = conditionCall)
