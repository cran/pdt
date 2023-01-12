#' @title AB_permutation_distancing_test
#'
#' @description Performs a randomisation test for two phases (A and B) that corrects for dependency.
#' The correcting is done through stepwise resampling the time series while varying
#' the distance between observations. The required distance 0,1,2,3.. is determined based
#' on repeated dependency testing while stepwise increasing the distance.
#' The input x and y values should be equidistant (with NA's included) using insert_NA_and_try_to_shift.
#' The distance per cycles = k-1.
#' k_max should be max 25% of (non missing) observations.
#' If de_A_trend=TRUE, phase A will be first de_A_trended.
#' If de_B_trend=TRUE, phase B will be first de_B_trended.
#' If detrend_x_position = "first" : take predicted value for first valid observation.
#' If detrend_x_position = "center": take predicted value for center observation.
#' If detrend_x_position = "last"  : take predicted value for last valid observation.
#' The p-value returned corresponds with the lowest Ljung-Box test (minimal) p-value found.
#' The statistic returned correspond with median chi-square k with p_box larger than alpha_p_box_test
#' or statistic_box smaller than max_statistic_p_box_test (i.e., not dependent).
#'
#' @param x factor vector to indicate conditions or phases (e.g., "A" and "B")
#' @param y numerical vector with the observed y-values
#' @param test_statistic character how to compute the test statistic c("A-B", "B-A", "*") *=two-sided
#' @param test_statistic_function character compute and compare "mean" or "median" for A and B
#' @param reps_max numerical maximum number of permutation replications (the theoretical number= n!)
#' @param k_max numerical maximum k value
#' @param alpha_p_box_test numerical see above
#' @param max_statistic_p_box_test numerical see above
#' @param no_duplicates boolean do a permutation test without duplicates (makes it much slower)
#' @param remove_NA_series_above_k boolean first clean the data by skipping repeated NA's
#' @param de_A_trend boolean de-trend A first (optional)
#' @param detrend_A_position character c("first", "center", "last"), see detrend_A
#' @param de_B_trend boolean de-trend B first (optional)
#' @param detrend_B_position character c("first", "center", "last"), see detrend_B
#' @param show_plot boolean show test plot of statistical test
#' @param show_plot_header character header of test plot
#'
#' @return List with the permutation distancing test results:
#'     de_A_trend setting in call,
#'     detrend_A_position in call,
#'     de_B_trend setting in call,
#'     detrend_B_position in call,
#'     ar1 = vector of computed ar1 values per distancing step (0,1,2, etc),
#'     p_box = vector of computed box-test p-values per distancing step,
#'     statistic_box = vector of box-test statistics per distancing step,
#'     observed_test_statistic = computed overall AB test statistic (before distancing),
#'     effect_size_overall = computed overall effect size (before distancing),
#'     p = vector of computed permutation test p-values per distancing step,
#'     effect_size vector of computed permutation test effect-sizes per distancing step,
#'     p_fitted = vector of lm-fitted line p-values through p_box,
#'     k_max = k_max setting in call or computed based on the number of observations,
#'     k_selected_based_on_Box_test = selected k values,
#'     p_selected_based_on_Box_test = selected p-value,
#'     effect_size_selected_based_on_Box_test = selected effect-size values.
#'
#' @examples
#' pdt::AB_permutation_distancing_test(
#'   as.factor(c(rep("A",20), rep("B",20))),
#'   c(rnorm(20), rnorm(20)+2),
#'   test_statistic="B-A",
#'   test_statistic_function="mean",
#'   reps_max=1000,
#'   k_max=NULL,
#'   alpha_p_box_test=0.1,
#'   max_statistic_p_box_test=2.7,
#'   no_duplicates=FALSE,
#'   remove_NA_series_above_k=TRUE,
#'   de_A_trend=FALSE,
#'   detrend_A_position="center",
#'   de_B_trend=FALSE,
#'   detrend_B_position="center",
#'   show_plot=FALSE,
#'   show_plot_header="")
#'
#' @export



AB_permutation_distancing_test <- function(x, y, test_statistic="*", test_statistic_function="mean", reps_max=2000,
                                           k_max=NULL, alpha_p_box_test=0.1, max_statistic_p_box_test=2.7, no_duplicates=FALSE, remove_NA_series_above_k=TRUE,
                                           de_A_trend=FALSE, detrend_A_position="center", de_B_trend=FALSE, detrend_B_position="center",
                                           show_plot=FALSE, show_plot_header="")  {


  #if (de_A_trend && length(y[x==levels(x)[1]][!is.na(y[x==levels(x)[1]])]) < 20 ||
  #    length(y[x==levels(x)[2]][!is.na(y[x==levels(x)[2]])]) < 20) warning("de_A_trend produces unreliable results in case of small (< 20) A or B number of observations")

  if (is.null(k_max)) { # auto_break when Ljung-Box test not significant
    k_max <- trunc(length(y[!is.na(y)])/4)
    auto_break <- TRUE
  } else {
    #test for too large setting of k_max; use minimum number of observations of at least 4 * k
    if (length(y[!is.na(y)])/k_max < 4){
      k_max <- trunc(length(y[!is.na(y)])/4)
      warning(paste("reducing k_max to 25% of (non missing) observations: ", as.character(k_max), sep=''))
    }
    auto_break <- FALSE
  }


  #de_A_trend and de_B_trend if requested; do it separately
  #plot(c(1:length(x)), y, type = "l", col="blue")
  if (de_A_trend) y <- detrend_A(x, c(1:length(x)), y, detrend_A_position=detrend_A_position)$y_detrended
  if (de_B_trend) y <- detrend_B(x, c(1:length(x)), y, detrend_B_position=detrend_B_position)$y_detrended
  #graphics::lines(c(1:length(x)), y, type = "l", col="red")

  test_results <- AB_permutation_test(x, y, test_statistic=test_statistic, test_statistic_function=test_statistic_function, reps_max=reps_max, no_duplicates=FALSE, show_plot=FALSE)
  observed_test_statistic <- test_results$observed_test_statistic
  effect_size_overall <- test_results$effect_size

  if (is.na(observed_test_statistic)) {

    ar1 <- NA
    p_box <- NA
    statistic_box <- NA
    observed_test_statistic <- NA
    effect_size_overall <- NA
    p <- NA
    effect_size <- NA
    p_fitted <- NA
    k_max <- NA
    k_selected_based_on_Box_test <- NA
    p_selected_based_on_Box_test <- NA
    effect_size_selected_based_on_Box_test <- NA

  } else {

    indices <- 1:length(x) # indices based on equidistant x,y

    p <- c()
    ar1 <- c()
    p_box <- c()
    statistic_box <- c()
    effect_size <- c()

    for (k in 1:k_max) {

      # k-1 is number of distancing NA's between observations; k=1 : no distancing
      # remove NA-series of more than k NA's in y
      # this will protect against inflation of p as a result of too many NA combinations during distancing
      NA_reduced <- do_remove_NA_series_above_k(x, y, k-1, x_values=indices)
      x_tmp <- NA_reduced$x
      y_tmp <- NA_reduced$y
      indices_tmp <- NA_reduced$x_values  # these indices are based on original equidistant y

      #print(k)
      #print(y_tmp)
      # n_k <- sample_size_overall/(k+1)
      # for each k, create k (l=1 to k) new distancing time series by NA intermediate values,
      # and then perform the random assignments

      random_assignments <- c()
      ar1_l <- c()
      p_box_l <- c()
      statistic_box_l <- c()
      effect_size_l <- c()

      for (l in 1:k) {

        indices_distancing <- indices[seq(l, length(y), by=k)] # these indices are based on original equidistant y
        indices_distancing_in_indices_tmp <- indices_distancing[indices_distancing %in% indices_tmp] # these indices are based on original equidistant y
        #print(indices_distancing)
        #print(indices_distancing_in_indices_tmp)

        #method without remove_NA_series_above_k
        if (!remove_NA_series_above_k) {
          x_distancing <- x
          y_distancing <- rep(NA, length(y))
          y_distancing[indices_distancing] <- y[indices_distancing]
        } else {
          #method with remove_NA_series_above_k
          x_distancing <- x_tmp
          y_distancing <- rep(NA, length(y_tmp))
          for (m in 1:length(indices_distancing_in_indices_tmp)) y_distancing[which(indices_tmp %in% indices_distancing_in_indices_tmp[m])] <- y[indices_distancing_in_indices_tmp[m]]
          #onderstaande is sneller, maar gaat het ook goed?
          #y_distancing <- rep(NA, length(y_tmp))
          #y_distancing[indices_tmp %in% indices_distancing_in_indices_tmp] <- y[indices_distancing_in_indices_tmp]
        }

        #print(y_distancing)
        plot_header <- paste("distancing k: ", as.character(k), "  series l: ", as.character(l), sep='')
        test_results_l <- AB_permutation_test(x_distancing, y_distancing, test_statistic=test_statistic, test_statistic_function=test_statistic_function, reps_max=reps_max, no_duplicates=no_duplicates,
                                              show_plot=FALSE)
        random_assignments <- c(random_assignments, test_results_l$random_assignments)
        ar1_l <- c(ar1_l, stats::acf(y_distancing[!is.na(y_distancing)], lag.max=1, type="correlation", plot=FALSE, na.action=stats::na.pass, demean=TRUE)$acf[2])

        p_box_l <- c(p_box_l, stats::Box.test(y_distancing[!is.na(y_distancing)], type = "Ljung-Box")$p.value)
        statistic_box_l <- c(statistic_box_l, stats::Box.test(y_distancing[!is.na(y_distancing)], type = "Ljung-Box")$statistic)

        effect_size_l <- c(effect_size_l, test_results_l$effect_size)

      }

      #print(ar1_l)
      ar1 <- c(ar1, stats::median(ar1_l, na.rm=TRUE)) # median because ar1 values are skewed
      p_box <- c(p_box, min(p_box_l, na.rm=TRUE)) #take smallest: one of them is significant
      statistic_box <- c(statistic_box, stats::median(statistic_box_l, na.rm=TRUE)) #median as well
      effect_size <- c(effect_size, mean(effect_size_l, na.rm=TRUE)) #mean

      # correct length(random_assignments) by number of NA's
      p_randomization_AB_k <- sum(random_assignments>observed_test_statistic, na.rm=TRUE) / length(random_assignments[!is.na(random_assignments)])
      p <- c(p, p_randomization_AB_k)

      if (auto_break && (p_box[length(p_box)]>alpha_p_box_test || statistic_box[length(statistic_box)]<max_statistic_p_box_test) ) break

    }

    p[is.nan(p)] <- NA

    if (!all(is.na(p))) {
      indices <- which(!is.na(p))
      p <- p[indices]
      ar1 <- ar1[indices]
      p_box <- p_box[indices]
      statistic_box <- statistic_box[indices]
      effect_size <- effect_size[indices]
      # get second order fitted regression line for p
      k <- 1:length(p) # length of vectors with computed pdt values
      df_tmp <- data.frame(k, ar1, p, effect_size)
      p_fitted <- as.numeric(stats::lm(formula = p ~ k + I(k^2), data=df_tmp)$fitted.values)
      effect_size_fitted <- as.numeric(stats::lm(formula = effect_size ~ k + I(k^2), data=df_tmp)$fitted.values)
      #get first k with Box test not significant or with small chi-square
      k_selected_based_on_Box_test <- which(p_box>alpha_p_box_test | statistic_box<max_statistic_p_box_test)[1]
    } else {
      p_fitted <- NA
      effect_size_fitted <- NA
      k_selected_based_on_Box_test <- NA
    }

    if (!is.na(k_selected_based_on_Box_test)) {
      p_selected_based_on_Box_test <- p_fitted[k_selected_based_on_Box_test]
      effect_size_selected_based_on_Box_test <- effect_size_fitted[k_selected_based_on_Box_test]
    } else {
      p_selected_based_on_Box_test <- NA
      effect_size_selected_based_on_Box_test <- NA
    }

    if (show_plot) {
      plot(k, ar1, type="l", col="blue", xlab="distancing k", ylab="ar1", main=show_plot_header)
      graphics::lines(k_selected_based_on_Box_test, ar1[k_selected_based_on_Box_test], type="p", col="red")
      plot(k, p_box, type="l", col="blue", xlab="distancing k", ylab="Ljung Box test p-value", main=show_plot_header)
      graphics::lines(k_selected_based_on_Box_test, p_box[k_selected_based_on_Box_test], type="p", col="red")
      plot(k, statistic_box, type="l", col="blue", xlab="distancing k", ylab="Ljung Box test statistic", main=show_plot_header)
      graphics::lines(k_selected_based_on_Box_test, statistic_box[k_selected_based_on_Box_test], type="p", col="red")
      plot(k, p, type="l", col="blue", xlab="distancing k", ylab="one-sided p-value", main=show_plot_header)
      graphics::lines(k, p_fitted, type="l", col="green")
      graphics::lines(k_selected_based_on_Box_test, p_selected_based_on_Box_test, type="p", col="red")
      plot(k, effect_size, type="l", col="blue", xlab="distancing k", ylab="effect size", main=show_plot_header)
      graphics::lines(k, effect_size_fitted, type="l", col="green")
      graphics::lines(k_selected_based_on_Box_test, effect_size_selected_based_on_Box_test, type="p", col="red")
    }

  }

  list(
    de_A_trend = de_A_trend,
    detrend_A_position = detrend_A_position,
    de_B_trend = de_B_trend,
    detrend_B_position = detrend_B_position,
    ar1 = ar1,
    p_box = p_box,
    statistic_box = statistic_box,
    observed_test_statistic = observed_test_statistic,
    effect_size_overall = effect_size_overall,
    p = p,
    effect_size = effect_size,
    p_fitted = p_fitted,
    k_max = k_max,
    k_selected_based_on_Box_test = k_selected_based_on_Box_test,
    p_selected_based_on_Box_test = p_selected_based_on_Box_test,
    effect_size_selected_based_on_Box_test = effect_size_selected_based_on_Box_test
  )

}

