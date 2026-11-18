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
# block. blocks whose elements are all zero are left as zeros.
# the result is in the same order as `x` (i.e. node order, not block order)
l1_normalize_within <- function(x, block) {
  totals <- stats::ave(x, block, FUN = sum)
  # blocks whose x are all zero stay zero rather than becoming NaN
  as.numeric(ifelse(totals == 0, 0, x / totals))
}

# put node-level weights `theta` into the non-zero entries of a block
# membership matrix `X` with exactly one non-zero entry per row. the
# entries of X@x are stored in column-major order, so they are not in node
# order unless the nodes are sorted by block. use the row indices X@i to
# put each weight in the right row
set_membership_weights <- function(X, theta) {
  if (methods::is(X, "CsparseMatrix")) {
    X@x <- theta[X@i + 1L]
  } else {
    X@x <- theta
  }
  X
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

check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(glue("`{arg}` must be `TRUE` or `FALSE`."), call. = FALSE)
  }

  invisible(x)
}

# block sizes for `n` nodes in `k` blocks that are as balanced as possible:
# every block gets floor(n / k) nodes, and the remainder is spread over the
# first blocks, one node each. this is the same rule `dcsbm()` uses by default
balanced_block_sizes <- function(n, k) {
  base_value <- floor(n / k)
  remainder <- n %% k
  c(rep(base_value + 1, remainder), rep(base_value, k - remainder))
}

utils::globalVariables(
  c("i", "j", "value", "x")
)
