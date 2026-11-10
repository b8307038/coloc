## What x is, for a message such as "bf2 is a character vector": "a
## data.frame", "a numeric vector", "a SuSiE fit" or "NULL"
.what_is <- function(x) {
    if (is.null(x))
        return("NULL")
    what <- if (inherits(x, "susie")) "SuSiE fit"
            else if (is.object(x) || !is.atomic(x)) class(x)[1]
            else paste(mode(x), if (is.matrix(x)) "matrix" else "vector")
    paste(if (grepl("^[aeiou]", what)) "an" else "a", what)
}
