% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AB_permutation_distancing_test.R
\name{AB_permutation_distancing_test}
\alias{AB_permutation_distancing_test}
\title{AB_permutation_distancing_test}
\usage{
AB_permutation_distancing_test(
  x,
  y,
  test_statistic = "*",
  test_statistic_function = "mean",
  reps_max = 2000,
  k_max = NULL,
  alpha_p_box_test = 0.1,
  max_statistic_p_box_test = 2.7,
  no_duplicates = FALSE,
  remove_NA_series_above_k = TRUE,
  de_A_trend = FALSE,
  detrend_A_position = "center",
  de_B_trend = FALSE,
  detrend_B_position = "center",
  show_plot = FALSE,
  show_plot_header = ""
)
}
\arguments{
\item{x}{factor vector to indicate conditions or phases (e.g., "A" and "B")}

\item{y}{numerical vector with the observed y-values}

\item{test_statistic}{character how to compute the test statistic c("A-B", "B-A", "*") *=two-sided}

\item{test_statistic_function}{character compute and compare "mean" or "median" for A and B}

\item{reps_max}{numerical maximum number of permutation replications (the theoretical number= n!)}

\item{k_max}{numerical maximum k value}

\item{alpha_p_box_test}{numerical see above}

\item{max_statistic_p_box_test}{numerical see above}

\item{no_duplicates}{boolean do a permutation test without duplicates (makes it much slower)}

\item{remove_NA_series_above_k}{boolean first clean the data by skipping repeated NA's}

\item{de_A_trend}{boolean de-trend A first (optional)}

\item{detrend_A_position}{character c("first", "center", "last"), see detrend_A}

\item{de_B_trend}{boolean de-trend B first (optional)}

\item{detrend_B_position}{character c("first", "center", "last"), see detrend_B}

\item{show_plot}{boolean show test plot of statistical test}

\item{show_plot_header}{character header of test plot}
}
\value{
List with the permutation distancing test results:
    de_A_trend setting in call,
    detrend_A_position in call,
    de_B_trend setting in call,
    detrend_B_position in call,
    ar1 = vector of computed ar1 values per distancing step (0,1,2, etc),
    p_box = vector of computed box-test p-values per distancing step,
    statistic_box = vector of box-test statistics per distancing step,
    observed_test_statistic = computed overall AB test statistic (before distancing),
    effect_size_overall = computed overall effect size (before distancing),
    p = vector of computed permutation test p-values per distancing step,
    effect_size vector of computed permutation test effect-sizes per distancing step,
    p_fitted = vector of lm-fitted line p-values through p_box,
    k_max = k_max setting in call or computed based on the number of observations,
    k_selected_based_on_Box_test = selected k values,
    p_selected_based_on_Box_test = selected p-value,
    effect_size_selected_based_on_Box_test = selected effect-size values.
}
\description{
Performs a randomisation test for two phases (A and B) that corrects for dependency.
The correcting is done through stepwise resampling the time series while varying
the distance between observations. The required distance 0,1,2,3.. is determined based
on repeated dependency testing while stepwise increasing the distance.
The input x and y values should be equidistant (with NA's included) using insert_NA_and_try_to_shift.
The distance per cycles = k-1.
k_max should be max 25% of (non missing) observations.
If de_A_trend=TRUE, phase A will be first de_A_trended.
If de_B_trend=TRUE, phase B will be first de_B_trended.
If detrend_x_position = "first" : take predicted value for first valid observation.
If detrend_x_position = "center": take predicted value for center observation.
If detrend_x_position = "last"  : take predicted value for last valid observation.
The p-value returned corresponds with the lowest Ljung-Box test (minimal) p-value found.
The statistic returned correspond with median chi-square k with p_box larger than alpha_p_box_test
or statistic_box smaller than max_statistic_p_box_test (i.e., not dependent).
}
\examples{
pdt::AB_permutation_distancing_test(
  as.factor(c(rep("A",20), rep("B",20))),
  c(rnorm(20), rnorm(20)+2),
  test_statistic="B-A",
  test_statistic_function="mean",
  reps_max=1000,
  k_max=NULL,
  alpha_p_box_test=0.1,
  max_statistic_p_box_test=2.7,
  no_duplicates=FALSE,
  remove_NA_series_above_k=TRUE,
  de_A_trend=FALSE,
  detrend_A_position="center",
  de_B_trend=FALSE,
  detrend_B_position="center",
  show_plot=FALSE,
  show_plot_header="")

}
