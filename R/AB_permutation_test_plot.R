#' @title AB_permutation_test_plot
#'
#' @description Creates a permutation distancing test plot.
#' Several plot options are available, e.g., to show both the observed and detrended lines.
#'
#' @param x factor vector to indicate conditions or phases (e.g., "A" and "B")
#' @param x_values numerical vector with distance (time markers) between observations
#' @param y numerical vector with the observed y-values
#' @param test_statistic_function character compute and compare "mean" or "median" for phase A and B
#' @param de_A_trend boolean de-trend A (optional)
#' @param detrend_A_position character c("first", "center", "last"), see detrend_A
#' @param show_de_A_trended boolean show de-trend A line (optional)
#' @param de_B_trend boolean de-trend B (optional)
#' @param detrend_B_position character c("first", "center", "last"), see detrend_B
#' @param show_de_B_trended boolean show de-trend B line (optional)
#' @param show_plot_header character header
#' @param xlab character x-axis label
#' @param ylab character y-axis label
#' @param ylim numerical vector of y-axis limits
#' @param labels character vector of labels
#' @param line_colors character vector with colors of the succeeding lines c("blue", "red", "blue", "red", "blue", "cyan"),
#' @param show_legend boolean show legend
#'
#' @return NULL.
#'
#' @examples
#' pdt::AB_permutation_test_plot(
#'   as.factor(c(rep("A",20), rep("B",20))),
#'   1:40,
#'   c(rnorm(20), rnorm(20)+2),
#'   test_statistic_function="mean",
#'   de_A_trend=TRUE,
#'   detrend_A_position="center",
#'   show_de_A_trended=TRUE,
#'   de_B_trend=TRUE,
#'   detrend_B_position="center",
#'   show_de_B_trended=TRUE,
#'   show_plot_header="",
#'   xlab="",
#'   ylab="",
#'   ylim=NULL,
#'   labels=NULL,
#'   line_colors=c("blue", "red", "blue", "red", "blue", "cyan"),
#'   show_legend=TRUE)
#'
#' @export


AB_permutation_test_plot <- function(x, x_values, y, test_statistic_function="mean",
                                     de_A_trend=FALSE, detrend_A_position="center", show_de_A_trended=FALSE,
                                     de_B_trend=FALSE, detrend_B_position="center", show_de_B_trended=FALSE,
                                     show_plot_header="", xlab="", ylab="",
                                     ylim=NULL, labels=NULL, line_colors=c("blue", "red", "blue", "red", "blue", "cyan"), show_legend=TRUE)  {

  # x= treatment condition as factor ("A" vs. "B"),
  # x_values = x-axis observation time values
  # y=observed values

  if (de_A_trend) {
    dt_A <- detrend_A(x, x_values, y, detrend_A_position=detrend_A_position)
    x_values_A_trend = dt_A$x_values_A_trend
    y_A_trend <- dt_A$y_A_trend
    y_de_A_trend <- dt_A$y_detrended
  }
  if (de_B_trend) {
    dt_B <- detrend_B(x, x_values, y, detrend_B_position=detrend_B_position)
    x_values_B_trend = dt_B$x_values_B_trend
    y_B_trend <- dt_B$y_B_trend
    y_de_B_trend <- dt_B$y_detrended
  }

  #remove NA-series of more than 1 NA's in y
  # this will improve the test-statistic lines
  NA_reduced <- do_remove_NA_series_above_k(x, y, 1, x_values)
  x_NA_reduced <- NA_reduced$x
  y_NA_reduced <- NA_reduced$y
  x_values_de_trend <- NA_reduced$x_values
  if (de_A_trend) {
    NA_reduced_de_A_trend <- do_remove_NA_series_above_k(x, y_de_A_trend, 1, x_values)
    y_de_A_trend <- NA_reduced_de_A_trend$y
  }
  if (de_B_trend) {
    NA_reduced_de_B_trend <- do_remove_NA_series_above_k(x, y_de_B_trend, 1, x_values)
    y_de_B_trend <- NA_reduced_de_B_trend$y
  }

  #just to make sure 1,2,3, etc. are used for levels
  for (i in 1:length(levels(x_NA_reduced))) levels(x_NA_reduced)[levels(x_NA_reduced)==levels(x_NA_reduced)[i]] <- as.character(i)

  # prepare to plot lines and points
  col_data <- line_colors
  if (de_A_trend | de_B_trend) col_detrend <- line_colors[length(line_colors)]

  i_col <- 1
  y_tmp <- y_NA_reduced
  if (de_A_trend) y_tmp_de_A_trend <- y_de_A_trend
  if (de_B_trend) y_tmp_de_B_trend <- y_de_B_trend
  x_tmp <- x_NA_reduced
  if (is.null(x_values)[1] || is.na(x_values)[1]) {
    x_values_tmp <- 1:length(x_NA_reduced)
  } else {
    x_values_tmp <- x_values_de_trend
  }
  length_first_level <- length(which(x_tmp==1))
  if (is.null(ylim)) ylim=c(min(y_tmp, na.rm=TRUE)-1, max(y_tmp, na.rm=TRUE)+1)
  xlim <- c(min(x_values_tmp, na.rm=TRUE), max(x_values_tmp, na.rm=TRUE))

  if (!de_A_trend | (de_A_trend & !show_de_A_trended)) {
    plot(x_values_tmp[1:length_first_level], y_tmp[1:length_first_level], type = "l", col=col_data[i_col], xlab=xlab, ylab=ylab, main=show_plot_header, xlim=xlim, ylim=ylim)
    graphics::lines(x_values_tmp[1:length_first_level], y_tmp[1:length_first_level], type = "p", col=col_data[i_col], cex=.5)
  } else {
    plot(x_values_tmp[1:length_first_level], y_tmp_de_A_trend[1:length_first_level], type = "l", col=col_detrend, xlab=xlab, ylab=ylab, main=show_plot_header, xlim=xlim, ylim=ylim)
    graphics::lines(x_values_tmp[1:length_first_level], y_tmp_de_A_trend[1:length_first_level], type = "p", col=col_detrend, cex=.5)
    graphics::lines(x_values_tmp[1:length_first_level], y_tmp[1:length_first_level], type = "l", col=col_data[i_col])
    graphics::lines(x_values_tmp[1:length_first_level], y_tmp[1:length_first_level], type = "p", col=col_data[i_col], cex=.5)
  }

  if (length(levels(x_NA_reduced))>1) {
    if (!de_B_trend | (de_B_trend & !show_de_B_trended)) {

      for (i in 2:length(levels(x_NA_reduced))) { # start with the second, the first is already plotted
        i_col <- i_col+1
        y_changepoint <- which(x_tmp==i)[1] # a changepoint is the first element with the next level

        y_tmp <- c(y_tmp, NA)
        y_tmp[(y_changepoint+1):length(y_tmp)] <- y_tmp[y_changepoint:(length(y_tmp)-1)]
        y_tmp[y_changepoint] <- NA
        x_tmp <- c(as.character(x_tmp), "NA")
        x_tmp[(y_changepoint+1):length(x_tmp)] <- as.character(x_tmp[y_changepoint:(length(x_tmp)-1)])
        x_tmp[y_changepoint] <- "NA"
        x_values_tmp <- c(x_values_tmp, NA)
        x_values_tmp[(y_changepoint+1):length(x_values_tmp)] <- x_values_tmp[y_changepoint:(length(x_values_tmp)-1)]
        x_values_tmp[y_changepoint] <- NA

        graphics::lines(x_values_tmp[(y_changepoint+1):length(y_tmp)], y_tmp[(y_changepoint+1):length(y_tmp)], type = "l", col=col_data[i_col])
        graphics::lines(x_values_tmp[(y_changepoint+1):length(y_tmp)], y_tmp[(y_changepoint+1):length(y_tmp)], type = "p", col=col_data[i_col], cex=.5)
      }
    } else {

      for (i in 2:length(levels(x_NA_reduced))) { # start with the second, the first is already plotted
        i_col <- i_col+1
        y_changepoint <- which(x_tmp==i)[1] # a changepoint is the first element with the next level

        y_tmp <- c(y_tmp, NA)
        y_tmp[(y_changepoint+1):length(y_tmp)] <- y_tmp[y_changepoint:(length(y_tmp)-1)]
        y_tmp[y_changepoint] <- NA
        y_tmp_de_B_trend <- c(y_tmp_de_B_trend, NA)
        y_tmp_de_B_trend[(y_changepoint+1):length(y_tmp_de_B_trend)] <- y_tmp_de_B_trend[y_changepoint:(length(y_tmp_de_B_trend)-1)]
        y_tmp_de_B_trend[y_changepoint] <- NA
        x_tmp <- c(as.character(x_tmp), "NA")
        x_tmp[(y_changepoint+1):length(x_tmp)] <- as.character(x_tmp[y_changepoint:(length(x_tmp)-1)])
        x_tmp[y_changepoint] <- "NA"
        x_values_tmp <- c(x_values_tmp, NA)
        x_values_tmp[(y_changepoint+1):length(x_values_tmp)] <- x_values_tmp[y_changepoint:(length(x_values_tmp)-1)]
        x_values_tmp[y_changepoint] <- NA

        graphics::lines(x_values_tmp[(y_changepoint+1):length(y_tmp_de_B_trend)], y_tmp_de_B_trend[(y_changepoint+1):length(y_tmp_de_B_trend)], type = "l", col=col_detrend)
        graphics::lines(x_values_tmp[(y_changepoint+1):length(y_tmp_de_B_trend)], y_tmp_de_B_trend[(y_changepoint+1):length(y_tmp_de_B_trend)], type = "p", col=col_detrend, cex=.5)

        graphics::lines(x_values_tmp[(y_changepoint+1):length(y_tmp)], y_tmp[(y_changepoint+1):length(y_tmp)], type = "l", col=col_data[i_col])
        graphics::lines(x_values_tmp[(y_changepoint+1):length(y_tmp)], y_tmp[(y_changepoint+1):length(y_tmp)], type = "p", col=col_data[i_col], cex=.5)
      }


    }

  }

  # if y-values with x=NA exist: they are not part of the AB dataset, plot them in light gray
  if (length(x_NA_reduced[is.na(x_NA_reduced)])>0) {
    y_with_nax <- y_NA_reduced[is.na(x_NA_reduced)]
    x_values_with_nax <- x_values[is.na(x_NA_reduced)]
    graphics::lines(x_values_with_nax, y_with_nax, type = "l", col="gray80")
    graphics::lines(x_values_with_nax, y_with_nax, type = "p", col="gray80", cex=.5)
  }

  i_col <- 1
  for (i in 1:length(levels(x_NA_reduced))) { # start with the first
    if (de_A_trend && i==1) {
      y2_tmp_de_A_trend <- rep(NA, length(y_tmp_de_A_trend))
      if (test_statistic_function == "median") {
        y2_tmp_de_A_trend[x_tmp==as.character(i)] <- stats::median(y_tmp_de_A_trend[x_tmp==as.character(i)], na.rm=TRUE)
      } else {
        y2_tmp_de_A_trend[x_tmp==as.character(i)] <- mean(y_tmp_de_A_trend[x_tmp==as.character(i)], na.rm=TRUE)
      }
      #graphics::lines(x_values_A_trend, y_A_trend, type = "l", col=col_data[i_col], lty=2)
      graphics::lines(x_values_A_trend, y_A_trend, type = "p", col=col_data[i_col], cex=.35)
      if (show_de_A_trended) {
        graphics::lines(x_values_tmp, y2_tmp_de_A_trend, type = "l", col=col_detrend, lty=2)
      }
    } else if (de_B_trend && i==2) {
      y2_tmp_de_B_trend <- rep(NA, length(y_tmp_de_B_trend))
      if (test_statistic_function == "median") {
        y2_tmp_de_B_trend[x_tmp==as.character(i)] <- stats::median(y_tmp_de_B_trend[x_tmp==as.character(i)], na.rm=TRUE)
      } else {
        y2_tmp_de_B_trend[x_tmp==as.character(i)] <- mean(y_tmp_de_B_trend[x_tmp==as.character(i)], na.rm=TRUE)
      }

      #graphics::lines(x_values_B_trend, y_B_trend, type = "l", col=col_data[i_col], lty=2)
      graphics::lines(x_values_B_trend, y_B_trend, type = "p", col=col_data[i_col], cex=.35)
      if (show_de_B_trended) {
        graphics::lines(x_values_tmp, y2_tmp_de_B_trend, type = "l", col=col_detrend, lty=2)
      }
    } else {
      y2_tmp <- rep(NA, length(y_tmp))
      if (test_statistic_function == "median") {
        y2_tmp[x_tmp==as.character(i)] <- stats::median(y_tmp[x_tmp==as.character(i)], na.rm=TRUE)
      } else {
        y2_tmp[x_tmp==as.character(i)] <- mean(y_tmp[x_tmp==as.character(i)], na.rm=TRUE)
      }
      graphics::lines(x_values_tmp, y2_tmp, type = "l", col=col_data[i_col], lty=3)
    }
    i_col <- i_col+1
  }

  if (show_legend) {
    if (!de_A_trend || !show_de_A_trended) {
      if (!is.null(labels)) graphics::legend("topright", "(x,y)", legend=labels, col=col_data[1:length(levels(x))], lty = 1:1, bty = "n", cex=0.8)
    } else {
      if (!is.null(labels)) graphics::legend("topright", "(x,y)", legend=c(labels, "detrended"), col=c(col_data[1:length(levels(x))], col_detrend), lty = 1:1, bty = "n", cex=0.8)
    }

  }

}

