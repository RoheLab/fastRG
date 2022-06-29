# based on https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
is_invertible <- function(X) {
  !inherits(try(solve(X), silent = TRUE), "try-error")
}

sort_by_all_columns <- function(X) {
  args <- as.list(as.data.frame(as.matrix(X)))
  args$decreasing <- TRUE

  X <- X[do.call(order, args),]
  X
}

l1_normalize <- function(x) x / sum(x)

# `x` and `block` should be vectors of the same length. `block` should be
# discrete and `x` numeric. normalizes `x` to sum to one within each distinct
# block. ideally only call on positive x to avoid divide by zero errors
l1_normalize_within <- function(x, block) {
  do.call(c, lapply(split(x, block), l1_normalize))
}
