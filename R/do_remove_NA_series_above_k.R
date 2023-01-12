#' @title do_remove_NA_series_above_k
#'
#' @description Remove series of more than k succeeding NA's in x, y, and x_values.
#' This function is recommended before performing a permutation distancing test.
#'
#' @param x factor vector to indicate conditions or phases (e.g., "A" and "B")
#' @param y numeric vector with the observed y-values
#' @param k maximum allowed number of NA's
#' @param x_values numerical vector with distance (time markers) between observations
#'
#' @return List with the modified x, y, x_values:
#'   x = factor vector with conditions (e.g., "A" and "B").
#'   y = vector with observed values.
#'   x_values = vector with distance (time markers) between observations x,y.
#'
#' @examples
#' pdt::do_remove_NA_series_above_k(as.factor(c("A","A","A","B","B","B")),
#'   c(1.1,NA,NA,7.1,8.3,9.8), 1, c(1,2,4,5,6,8))
#'
#' @export


do_remove_NA_series_above_k <- function(x, y, k, x_values=NULL) {

  if (k > 0) {
    number_NA <- 0
    start_remove_indices <- c()
    end_remove_indices <- c()
    started <- FALSE
    for (m in 1:length(y)) {
      if (is.na(y[m]) && m<length(y)) {
        number_NA <- number_NA+1
        if (number_NA>k && !started) {
          start_remove_indices <- c(start_remove_indices, m)
          started <- TRUE
        }
      } else if (started && is.na(y[m]) && m==length(y)) {
        end_remove_indices <- c(end_remove_indices, m)
      } else {
        number_NA <- 0
        if (started) end_remove_indices <- c(end_remove_indices, m-1)
        started <- FALSE
      }
    }
  }

  if (k > 0 && length(start_remove_indices) > 0 ) {
    indices_to_remove <- c()
    for (m in 1:length(start_remove_indices)) indices_to_remove <- c(indices_to_remove, c(start_remove_indices[m]:end_remove_indices[m]))
    x_NA_series_removed <- x[-indices_to_remove]
    y_NA_series_removed <- y[-indices_to_remove]
    if (!is.null(x_values)) x_values_NA_series_removed <- x_values[-indices_to_remove]
  } else if (k == 0) {
    x_NA_series_removed <- x[!is.na(y)]
    y_NA_series_removed <- y[!is.na(y)]
    if (!is.null(x_values)) x_values_NA_series_removed <- x_values[!is.na(y)]
  } else {
    x_NA_series_removed <- x
    y_NA_series_removed <- y
    if (!is.null(x_values)) x_values_NA_series_removed <- x_values
  }

  if (is.null(x_values)) x_values_NA_series_removed <- NA

  list(
    x = x_NA_series_removed,
    y = y_NA_series_removed,
    x_values = x_values_NA_series_removed
  )
}

