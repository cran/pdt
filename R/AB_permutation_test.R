#' @title AB_permutation_test
#'
#' @description Performs a regular permutations test for two conditions or phases (A and B).
#'
#' @param x factor vector to indicate conditions or phases (e.g., "A" and "B")
#' @param y numerical vector with the observed y-values
#' @param test_statistic character how to compute the test statistic c("A-B", "B-A", "*") *=two-sided
#' @param test_statistic_function character compute and compare "mean" or "median" for A and B
#' @param reps_max numerical maximum number of permutation replications (the theoretical number= n!)
#' @param no_duplicates boolean do a permutation test without duplicates (makes it much slower)
#' @param show_plot boolean show test plot of statistical test
#' @param show_plot_header character header of test plot
#'
#' @return List with the permutation test results:
#'   observed_test_statistic = computed test statistic,
#'   effect_size = computed effect size (similar to Cohen's d),
#'   random_assignments,
#'   p_randomization_AB = p value randomization AB test,
#'   one_sided_p = one-sided p-value in case of B-A or A-B.
#'
#' @examples
#' pdt::AB_permutation_test(
#'   as.factor(c(rep("A",20), rep("B",20))),
#'   c(rnorm(20), rnorm(20)+2),
#'   test_statistic="B-A",
#'   test_statistic_function="mean",
#'   reps_max=1000,
#'   no_duplicates=FALSE,
#'   show_plot=FALSE,
#'   show_plot_header="")
#'
#' @export


AB_permutation_test <- function(x, y, test_statistic="*", test_statistic_function="mean", reps_max=2000, no_duplicates=FALSE, show_plot=FALSE, show_plot_header="")  {

  # do a permutation test without duplicates
  # theoretical number of permutations = n! = factorial(length(y))

  #return NA below this threshold
  MIN_OBSERVATIONS <- 2
  MAX_COUNT <- reps_max*100

  rperm <- function(m, size, max_count=1000) { # Obtain m unique permutations of 1:size
    # Returns a `size` by `m` matrix; each column is a permutation of 1:size.

    # Function to add a new permutation of 1:size as a new column to the permutations matrix
    newperm <- function(permutations, i, size, max_count=1000) {

      count <- 0  # Protects against infinite loops
      p <- sample(1:size)
      is.new <- FALSE

      # Generate another permutation column and check against previous ones.
      while (count < max_count && !is.new) {

        for (j in 1:i) {
          if (!all(is.na(permutations[, j])) && all(permutations[, j] == p)) {
            is.new <- FALSE
            #message("not new")
            p <- sample(1:size)
            break
          } else {
            count <- count+1
          }
        }
        if (j==i) {
          is.new <- TRUE
          permutations[, i] <- p
        }

      }
      if (count>=max_count) warning("maxcount reached in generating permutations without duplicates")
      permutations

    }

    # Obtain m unique permutations.
    permutations <- matrix(nrow=size, ncol=m) #start with empty matrix filled with NA
    for (i in 1:m) permutations <- newperm(permutations, i, size, max_count=max_count)

    permutations

  }


  if (length(levels(x)) > 2) stop("More than two levels found, only two allowed")
  if (!all(x == x[order(x)])) stop("No AB order found, only AB designs allowed")

  A <- y[x==levels(x)[1]]
  B <- y[x==levels(x)[2]]
  sample_size_overall <- length(x)

  if (all(is.na(A)) || all(is.na(B)) || length(A[!is.na(A)])<MIN_OBSERVATIONS || length(B[!is.na(B)])<MIN_OBSERVATIONS ||
      (stats::sd(A, na.rm=TRUE)+stats::sd(B, na.rm=TRUE))==0) {

    observed_test_statistic <- NA
    effect_size <- NA
    random_assignments <- c()
    p_randomization_AB <- NA

  } else {

    #observed_test_statistic_sign is used for unknown test_statistic
    if (test_statistic_function=="median") {
      if (test_statistic=="B-A") {
        observed_test_statistic <- stats::median(B, na.rm=TRUE)-stats::median(A, na.rm=TRUE)
      } else if (test_statistic=="A-B") {
        observed_test_statistic <- stats::median(A, na.rm=TRUE)-stats::median(B, na.rm=TRUE)
      } else {
        observed_test_statistic <- stats::median(A, na.rm=TRUE)-stats::median(B, na.rm=TRUE)
        if (observed_test_statistic < 0) observed_test_statistic_sign <- -1 else observed_test_statistic_sign <- 1
        observed_test_statistic <- observed_test_statistic_sign*observed_test_statistic
      }
    } else {
      if (test_statistic=="B-A") {
        observed_test_statistic <- mean(B, na.rm=TRUE)-mean(A, na.rm=TRUE)
      } else if (test_statistic=="A-B") {
        observed_test_statistic <- mean(A, na.rm=TRUE)-mean(B, na.rm=TRUE)
      } else {
        observed_test_statistic <- mean(A, na.rm=TRUE)-mean(B, na.rm=TRUE)
        if (observed_test_statistic < 0) observed_test_statistic_sign <- -1 else observed_test_statistic_sign <- 1
        observed_test_statistic <- observed_test_statistic_sign*observed_test_statistic
      }
    }
    sd_pooled <- (stats::sd(A, na.rm=TRUE)+stats::sd(B, na.rm=TRUE))/2
    effect_size <- observed_test_statistic/sd_pooled

    if (factorial(sample_size_overall) < reps_max) reps <- factorial(sample_size_overall) else reps <- reps_max

    #set.seed(seed = 14412)
    if (no_duplicates) {
      perm_matrix <- rperm(reps, size=sample_size_overall, max_count=MAX_COUNT) # Obtain m unique permutations; each column is a permutation of 1:size
      random_assignments <- c()
      for (i in 1:ncol(perm_matrix)) {
        perm_indices_i <- perm_matrix[,i]
        perm_indices_i_A <- perm_indices_i[1:length(A)]
        perm_indices_i_B <- perm_indices_i[(length(A)+1):(length(A)+length(B))]
        if (test_statistic_function=="median") {
          if (test_statistic=="B-A") {
            random_assignments_i <- stats::median(y[perm_indices_i_B], na.rm=TRUE)-stats::median(y[perm_indices_i_A], na.rm=TRUE)
          } else if (test_statistic=="A-B") {
            random_assignments_i <- stats::median(y[perm_indices_i_A], na.rm=TRUE)-stats::median(y[perm_indices_i_B], na.rm=TRUE)
          } else {
            random_assignments_i <- observed_test_statistic_sign*(stats::median(y[perm_indices_i_A], na.rm=TRUE)-stats::median(y[perm_indices_i_B], na.rm=TRUE))
          }
        } else {
          if (test_statistic=="B-A") {
            random_assignments_i <- mean(y[perm_indices_i_B], na.rm=TRUE)-mean(y[perm_indices_i_A], na.rm=TRUE)
          } else if (test_statistic=="A-B") {
            random_assignments_i <- mean(y[perm_indices_i_A], na.rm=TRUE)-mean(y[perm_indices_i_B], na.rm=TRUE)
          } else {
            random_assignments_i <- observed_test_statistic_sign*(mean(y[perm_indices_i_A], na.rm=TRUE)-mean(y[perm_indices_i_B], na.rm=TRUE))
          }
        }
        random_assignments <- c(random_assignments, random_assignments_i)
      }
    } else {
      random_assignments <- numeric(reps)
      for (i in 1:reps) {
        temp <- sample(y)
        if (test_statistic_function=="median") {
          if (test_statistic=="B-A") {
            random_assignments[i] <- stats::median(temp[(length(A)+1):length(y)], na.rm=TRUE)-stats::median(temp[1:length(A)], na.rm=TRUE)
          } else if (test_statistic=="A-B") {
            random_assignments[i] <- stats::median(temp[1:length(A)], na.rm=TRUE)-stats::median(temp[(length(A)+1):length(y)], na.rm=TRUE)
          } else {
            random_assignments[i] <- observed_test_statistic_sign*(stats::median(temp[1:length(A)], na.rm=TRUE)-stats::median(temp[(length(A)+1):length(y)], na.rm=TRUE))
          }
        } else {
          if (test_statistic=="B-A") {
            random_assignments[i] <- mean(temp[(length(A)+1):length(y)], na.rm=TRUE)-mean(temp[1:length(A)], na.rm=TRUE)
          } else if (test_statistic=="A-B") {
            random_assignments[i] <- mean(temp[1:length(A)], na.rm=TRUE)-mean(temp[(length(A)+1):length(y)], na.rm=TRUE)
          } else
            random_assignments[i] <- observed_test_statistic_sign*(mean(temp[1:length(A)], na.rm=TRUE)-mean(temp[(length(A)+1):length(y)], na.rm=TRUE))
        }
      }
    }

    random_assignments[is.na(random_assignments)] <- NA

    if (show_plot) {
      graphics::hist(random_assignments, main=show_plot_header)
      graphics::abline(v=observed_test_statistic, col="red")
    }

    # correct reps by number of NA's
    p_randomization_AB <- sum(random_assignments>observed_test_statistic, na.rm=TRUE)/(reps-length(which(is.na(random_assignments))))

  }

  one_sided_p <- NA
  if (test_statistic=="B-A" || test_statistic=="A-B") one_sided_p <- p_randomization_AB

  list(
    observed_test_statistic = observed_test_statistic,
    effect_size = effect_size,
    random_assignments = random_assignments,
    p_randomization_AB = p_randomization_AB,
    one_sided_p = one_sided_p
  )

}

