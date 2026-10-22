## Shared by the test files: a dataset in other forms, and what a call gives.

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

## 6.0.3's messages with the prefix as check_dataset() gives it now. "dataset N: "
## stays, "dataset : " goes, and the four messages that named no dataset, and
## the minimum p-value warning, gain "dataset N: ". With word "trait", as
## coloc() gives them, "trait N: " in place of "dataset N: ". suffix is recycled
## over msg, and NA stays NA. A message it doesn't know stops it, so a mistake here, or
## one of R's own messages reaching it, can't pass unnoticed
reprefix <- function(msg, suffix, word = "dataset") {
    bare <- c("N should be a single numeric", "LD not square", "LD rownames != colnames",
              "colnames in LD do not contain all SNPs", "minimum p value is: ")
    suffix <- rep_len(as.character(suffix), length(msg))
    vapply(seq_along(msg), function(i) {
        if (is.na(msg[i]))
            return(NA_character_)
        prefix <- if (suffix[i] == "") "" else paste0(word, " ", suffix[i], ": ")
        old <- paste0("dataset ", suffix[i], ": ")
        if (startsWith(msg[i], old))
            return(paste0(prefix, substring(msg[i], nchar(old) + 1)))
        if (any(startsWith(msg[i], bare)))
            return(paste0(prefix, msg[i]))
        stop("reprefix() doesn't know this message: ", msg[i])
    }, "")
}
