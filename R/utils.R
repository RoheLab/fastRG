# based on https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
is_invertible <- function(X) {
  !inherits(try(solve(X), silent = TRUE), "try-error")
}

sort_by_all_columns <- function(X) {
  args <- as.list(as.data.frame(as.matrix(X)))
  args$decreasing <- TRUE

  X <- X[do.call(order, args), ]
  X
}

sort_by_all_columns_idx <- function(X) {
  args <- as.list(as.data.frame(as.matrix(X)))
  args$decreasing <- TRUE

  do.call(order, args)
}

l1_normalize <- function(x) x / sum(x)

# `x` and `block` should be vectors of the same length. `block` should be
# discrete and `x` numeric. normalizes `x` to sum to one within each distinct
# block. ideally only call on positive x to avoid divide by zero errors
l1_normalize_within <- function(x, block) {
  do.call(c, lapply(split(x, block), l1_normalize))
}

check_multigraph <- function(multigraph, call = rlang::caller_env()) {
  if (!rlang::is_bool(multigraph)) {
    rlang::abort(
      paste0(
        "`multigraph` must be a single `TRUE` or `FALSE`, not ",
        obj_type(multigraph),
        "."
      ),
      call = call
    )
  }
  invisible(multigraph)
}

obj_type <- function(x) {
  if (is.null(x)) {
    return("`NULL`")
  }
  if (length(x) != 1) {
    return(paste0("a ", class(x)[1], " vector of length ", length(x)))
  }
  paste0("`", deparse(x), "`")
}

# collapse repeated (from, to) rows of an edgelist into a single row with
# an integer `weight` column counting the repeats. for undirected
# edgelists, callers must ensure from <= to beforehand so that (i, j) and
# (j, i) are collapsed together
collapse_multiedges <- function(edgelist) {
  edgelist <- tibble::as_tibble(edgelist)
  dplyr::count(edgelist, .data$from, .data$to, name = "weight")
}

utils::globalVariables(
  c("i", "j", "value", "x")
)
