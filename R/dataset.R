## A coloc dataset constructor

## The elements with one value per snp, in the order the lengths rule compares
## them and print() lists them
.per_snp_items <- c("beta", "varbeta", "pvalues", "MAF", "snp", "position", "z")

## The elements a coloc dataset may have. check_dataset() checks these and
## ignores anything else.
.dataset_items <- c(.per_snp_items, "N", "type", "s", "sdY", "LD")

## The elements each of coloc()'s methods needs, as check_dataset()'s req.
## "single" is also what a dataset is checked for when it is made,
## check_dataset()'s default. A req naming LD, as "susie"'s does, also brings
## the checks on LD and z that SuSiE needs.
.dataset_uses <- list(single = c("type", "snp"),
                      susie  = c("beta", "varbeta", "LD", "snp", "N"))

## The notes already given this session, by key: "plain_list" for the one
## coloc() and wakefield_bf() give when they convert a plain list
.coloc_state <- new.env(parent = emptyenv())
.coloc_state$noted <- character(0)

## Warn once per session for key. The key is recorded first, so a note turned
## into an error with options(warn = 2) still counts as given
.note_once <- function(key, message) {
    if (key %in% .coloc_state$noted)
        return(invisible(FALSE))
    .coloc_state$noted <- c(.coloc_state$noted, key)
    warning(message, call. = FALSE)
    invisible(TRUE)
}

## the plain-list note from fn, naming the traits given as plain lists, if any
.plain_list_note <- function(fn, traits = NULL) {
    given <- if (length(traits) > 1L)
                 paste("plain lists for traits", paste(traits, collapse = " and "),
                       "and converted them")
             else paste0("a plain list", if (length(traits)) paste(" for trait", traits),
                         " and converted it")
    paste0(fn, "() was given ", given, " with as_coloc_dataset().\n",
           "Make datasets with coloc_dataset() or as_coloc_dataset().\n",
           "(shown once per session)")
}

##' Make a coloc dataset
##'
##' `coloc_dataset()` builds a dataset from its parts and checks it with
##' [check_dataset()], so a mistake is reported where the dataset is made
##' rather than later, inside `coloc()`. The result has class
##' `"coloc_dataset"`, but is still a list with the same elements, so code that
##' reads `D$beta` works unchanged.
##'
##' `as_coloc_dataset()` does the same for a dataset that is already a list, as
##' made for earlier versions of coloc. Elements coloc does not use are kept.
##'
##' Which elements are needed depends on the data, and [check_dataset()] gives
##' the details. In brief: always give `snp` and `type`; give `beta` and
##' `varbeta` if you have them, or `pvalues`, `MAF` and `N` if not; give `s`
##' for a case-control trait, and `sdY` (or `MAF` and `N`, to estimate it) for
##' a quantitative one.
##'
##' A dataset is checked when it is made. If you change its elements
##' afterwards, coloc's functions check it again when they use it, but don't
##' do that. `as_coloc_dataset()` returns a coloc_dataset as it is:
##' [check_dataset()] checks one again for errors, and
##' `as_coloc_dataset(unclass(x))` makes it again.
##'
##' A `type` given as a factor, or with a name, comes back as a plain string.
##' 
##' @param snp character vector of snp ids
##' @param type `"quant"` for a quantitative trait, `"cc"` for
##'     case-control
##' @param beta,varbeta per-snp effect estimates, and their variances
##'     (std error squared)
##' @param pvalues per-snp p values
##' @param MAF per-snp minor allele frequencies
##' @param position per-snp positions, used for plotting
##' @param N number of samples
##' @param s for a case-control trait, the proportion of samples that
##'     are cases
##' @param sdY for a quantitative trait, the standard deviation of the
##'     trait
##' @param LD matrix of LD between snps, with snp ids as row and
##'     column names
##' @param z per-snp z-scores, optional: [runsusie()] fine-maps from them in
##'     place of `beta / sqrt(varbeta)`
##' @param warn.minp warn if no p value is smaller than this
##' @return a list of class `"coloc_dataset"`
##' @seealso [check_dataset()] for what each element must look like
##' @export
##' @examples
##' D <- coloc_dataset(snp = paste0("s", 1:5),
##'                    type = "quant",
##'                    beta = c(0.1, 0.5, 0.2, 0, 0.05),
##'                    varbeta = rep(0.01, 5),
##'                    sdY = 1)
##' D
##'
##' ## a dataset made as a list for an earlier version of coloc
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##'
##' ## mistakes are reported when the dataset is made
##' try(coloc_dataset(snp = paste0("s", 1:5),
##'                   type = "quant",
##'                   beta = c(0.1, 0.5, 0.2, 0, 0.05)))
coloc_dataset <- function(snp = NULL, type = NULL, beta = NULL, varbeta = NULL,
                          pvalues = NULL, MAF = NULL,
                          position = NULL, N = NULL, s = NULL, sdY = NULL,
                          LD = NULL, z = NULL, warn.minp = 1e-6) {
    ## every argument defaults to NULL so that a missing snp or type gets
    ## check_dataset()'s message, not R's "argument is missing"
    d <- list(beta = beta, varbeta = varbeta, pvalues = pvalues, MAF = MAF,
              snp = snp, position = position,
              N = N, type = type, s = s, sdY = sdY, LD = LD, z = z)
    d <- d[!sapply(d, is.null)]
    .new_coloc_dataset(d, warn.minp = warn.minp, call = sys.call())
}

##' @rdname coloc_dataset
##' @param x a list with the elements of a coloc dataset, or a coloc_dataset
##' @export
as_coloc_dataset <- function(x, warn.minp = 1e-6) {
    if (inherits(x, "coloc_dataset")) return(x)
    .new_coloc_dataset(x, warn.minp = warn.minp, call = sys.call())
}

## The one place a dataset gets the coloc_dataset class, so every
## coloc_dataset has been checked: when coloc_dataset() or as_coloc_dataset()
## makes it, and each time coloc() is given it. `call` is the call R shows in
## an error or warning: coloc_dataset(...), as_coloc_dataset(x) or coloc(...)
## as the user typed it, or as_coloc_dataset(d) when wakefield_bf() converts a
## list for them. coloc() also gives the trait's number, which its messages
## begin with as "trait N: ", and its method, which says what the dataset
## needs. A coloc_dataset keeps its class through the check, so it doesn't get
## the minimum p-value warning again: that was given when it was made.
.new_coloc_dataset <- function(d, warn.minp = 1e-6, call, trait = NULL,
                               use = "single", word = "trait") {
    ## a data.frame and a susie fit are lists too, but neither is a dataset
    if (!is.list(d) || is.data.frame(d) || inherits(d, "susie"))
        stop(simpleError(paste0(.prefix(trait, word),
                                "expected a coloc dataset: a list, or the result of coloc_dataset(), not a ",
                                class(d)[1]), call))
    if (!inherits(d, "coloc_dataset"))
        d <- unclass(d)
    .gate(d, suffix = trait, req = .dataset_uses[[use]], warn.minp = warn.minp,
          call = call, word = word)
    ## read by its label everywhere, so a factor or a name only gets in the way
    if (!is.null(d[["type"]]))
        d[["type"]] <- as.character(d[["type"]])
    structure(d, class = c("coloc_dataset", "list"))
}

##' @export
print.coloc_dataset <- function(x, ...) {
    ## the type by its label, as check_dataset() reads it: identical() would
    ## also compare a factor's or a named type's attributes
    cat("coloc dataset:",
        if (isTRUE(x$type == "cc")) "case-control" else "quantitative",
        "trait,", length(x$snp), "snps\n")
    cat("  per-snp:", paste(intersect(setdiff(.per_snp_items, "snp"), names(x)), collapse = ", "), "\n")
    sc <- unlist(x[intersect(c("N", "s", "sdY"), names(x))])
    if (length(sc))
        cat("  ", paste(names(sc), signif(sc, 4), sep = " = ", collapse = ", "),
            "\n", sep = "")
    ## [[ matches the name exactly, where $ would take an element like LD_ref
    if (!is.null(x[["LD"]])) cat("  LD:", nrow(x[["LD"]]), "x", ncol(x[["LD"]]), "\n")
    other <- setdiff(names(x), .dataset_items)
    if (length(other)) cat("  other elements:", paste(other, collapse = ", "), "\n")
    invisible(x)
}
