#' Convert level set flat index (lsfi) to multi-index (lsmi)
#'
#' @param lsfi Level set flat index.
#' @param num_intervals Number of intervals.
#' @return A multi-index corresponding to the flat index.
#' @export
to_lsmi <- function(lsfi, num_intervals) {
  # level set flat index (lsfi)
  j <- c(1, num_intervals) # put 1 in front to make indexing easier in the product prod(j[1:k])
  f <- c()
  for (k in 1:length(num_intervals)) {
    # use lsfi-1 to shift from 1-based indexing to 0-based indexing
    f[k] <- floor((lsfi-1) / prod(j[1:k])) %% num_intervals[k]
  }
  # lsmi = f+1 = level set multi index
  return(f+1) # shift from 0-based indexing back to 1-based indexing
}

#' Convert level set multi-index (lsmi) to flat index (lsfi)
#'
#' @param lsmi Level set multi-index.
#' @param num_intervals Number of intervals.
#' @return A flat index corresponding to the multi-index.
#' @export
to_lsfi <- function(lsmi, num_intervals) {
  # level set multi index (lsmi)
  lsfi <- lsmi[1]
  if (length(num_intervals) > 1) {
    for (i in 2:length(num_intervals)) {
      lsfi <- lsfi + prod(num_intervals[1:(i-1)]) * (lsmi[i]-1)
    }
  }
  return(lsfi)
}