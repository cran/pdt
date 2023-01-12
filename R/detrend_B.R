#' @title detrend_B
#'
#' @description Detrends the phase B part of time series y.
#' detrend B is optional and not validated.
#' The mean of the detrended signal will be set to the predicted value based on detrend_B_position:
#' detrend_B_position = "first" : take predicted value for first valid observation
#' detrend_B_position = "center" : take predicted value for center observation
#' detrend_B_position = "last" : take predicted value for last valid observation.
#'
#' @param x factor vector to indicate conditions or phases (e.g., "A" and "B")
#' @param x_values numerical vector with distance (time markers) between observations
#' @param y numeric vector with the observed y-values
#' @param detrend_B_position character to indicate the mean
#'
#' @return List with the trend and the detrended y values:
#'   x_values_B_trend = vector with distance (time markers) between B-detrended signal,
#'   y_B_trend = vector with computed B-trend,
#'   y_detrended = vector with computed B-detrended y values.
#'
#' @examples
#' pdt::detrend_B(as.factor(c(rep("A",20), rep("B",20))), 1:40,
#'   c(rnorm(20), rnorm(20)+2), detrend_B_position="center")
#'
#' @export


detrend_B <- function(x, x_values, y, detrend_B_position="center") {

  if (is.null(x_values)) x_values <- 1:length(x)

  length_B <- length(which(x==x[length(x)]))
  x_B <- 1:length_B
  y_B <- y[which(x==x[length(x)])[1]:length(y)]
  B_df <- data.frame(x_B, y_B)

  B_lm <- stats::lm(y_B ~ x_B, na.action=stats::na.exclude, data=B_df)
  B_residuals <- as.vector(stats::residuals(B_lm))
  y_B_trend <- as.vector(stats::predict(B_lm))
  x_values_B_trend <- x_values[which(x==x[length(x)])[1]:length(x_values)]

  valid_y_B <- y_B[!is.na(y_B)]
  valid_x_B <- x_B[!is.na(y_B)]
  if (detrend_B_position=="first") {
    first_valid_x <- valid_x_B[1]
    first_B_predicted <- y_B_trend[first_valid_x]
    y_B_detrended <- B_residuals+first_B_predicted
  } else if (detrend_B_position=="last") {
    last_valid_x <- valid_x_B[length(valid_x_B)]
    last_B_predicted <- y_B_trend[last_valid_x]
    y_B_detrended <- B_residuals+last_B_predicted
  } else {
    center_valid_x <- valid_x_B[trunc(length(valid_x_B)/2)]
    center_B_predicted <- y_B_trend[center_valid_x]
    y_B_detrended <- B_residuals+center_B_predicted
  }


  #plot(y, col="blue", type="p")

  y[which(x==x[length(x)])[1]:length(y)] <- y_B_detrended

  #lines(y, col="red", type="p")
  list(
    x_values_B_trend = x_values_B_trend,
    y_B_trend = y_B_trend,
    y_detrended = y
  )

}


