#' @title detrend_A
#'
#' @description Detrends the phase A part of time series y.
#' detrend A is optional and not validated.
#' The mean of the detrended signal will be set to the predicted value based on detrend_A_position:
#' detrend_A_position = "first" : take predicted value for first valid observation
#' detrend_A_position = "center" : take predicted value for center observation
#' detrend_A_position = "last" : take predicted value for last valid observation.
#'
#' @param x factor vector to indicate conditions or phases (e.g., "A" and "B")
#' @param x_values numerical vector with distance (time markers) between observations
#' @param y numeric vector with the observed y-values
#' @param detrend_A_position character to indicate the mean
#'
#' @return List with the trend and the detrended y-values:
#'   x_values_A_trend = vector with distance (time markers) between A-detrended signal,
#'   y_A_trend = vector with computed A-trend,
#'   y_detrended = vector with computed A-detrended y values.
#'
#' @examples
#' pdt::detrend_A(as.factor(c(rep("A",20), rep("B",20))), 1:40,
#'   c(rnorm(20), rnorm(20)+2), detrend_A_position="center")
#'
#' @export


detrend_A <- function(x, x_values, y, detrend_A_position="center") {

  if (is.null(x_values)) x_values <- 1:length(x)

  length_A <- length(which(x==x[1]))
  x_A <- 1:length_A
  y_A <- y[1:length_A]
  A_df <- data.frame(x_A, y_A)

  A_lm <- stats::lm(y_A ~ x_A, na.action=stats::na.exclude, data=A_df)
  A_residuals <- as.vector(stats::residuals(A_lm))
  y_A_trend <- as.vector(stats::predict(A_lm))
  x_values_A_trend <- x_values[1:length_A]

  valid_y_A <- y_A[!is.na(y_A)]
  valid_x_A <- x_A[!is.na(y_A)]
  if (detrend_A_position=="first") {
    first_valid_x <- valid_x_A[1]
    first_A_predicted <- y_A_trend[first_valid_x]
    y_A_detrended <- A_residuals+first_A_predicted
  } else if (detrend_A_position=="last") {
    last_valid_x <- valid_x_A[length(valid_x_A)]
    last_A_predicted <- y_A_trend[last_valid_x]
    y_A_detrended <- A_residuals+last_A_predicted
  } else {
    center_valid_x <- valid_x_A[trunc(length(valid_x_A)/2)]
    center_A_predicted <- y_A_trend[center_valid_x]
    y_A_detrended <- A_residuals+center_A_predicted
  }

  y[1:length_A] <- y_A_detrended

  #plot(y, col="blue", type="p")
  #lines(y, col="red", type="p")
  list(
    x_values_A_trend = x_values_A_trend,
    y_A_trend = y_A_trend,
    y_detrended = y
  )

}

