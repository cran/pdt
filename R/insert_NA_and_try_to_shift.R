#' @title insert_NA_and_try_to_shift
#'
#' @description Makes the input time series equidistant.
#' This is recommended before performing a permutation distancing test.
#' This function first inserts NA's for missing x_values,
#' then it tries to shift double value to previous or next NA's,
#' finally it aggregates the remaining identical x_values.
#'
#' @param x factor vector to indicate conditions or phases (e.g., "A" and "B")
#' @param x_values numerical vector with distance (time markers) between observations
#' @param y numeric vector with the observed y-values
#'
#' @return List with the modified x, x_values, y:
#'   x = factor vector with conditions or phases (e.g., "A" and "B").
#'   x_values = (optional) vector with distance (time markers) between observations.
#'   y = vector with observed values.
#'
#' @examples
#' pdt::insert_NA_and_try_to_shift(as.factor(c("A","A","A","B","B","B")),
#'   c(1,2,4,5,6,8), c(1.1,3.2,5.3,7.1,8.3,9.8))
#'
#' @export


insert_NA_and_try_to_shift <- function(x, x_values, y)  {

  if (is.null(x_values)) {

    x_values_new <- NULL
    x_new <- x
    y_new <- y

  } else if (!all(trunc(x_values)==x_values)) {

    # insert NA's for missing x_values
    x_values_all <- min(trunc(x_values)):max(trunc(x_values))
    x_values_new <- c()
    x_new <- c()
    y_new <- c()
    for (n in x_values_all) {
      indices <- which(trunc(x_values)==n)
      if (length(indices)==0) { # no values, insert NA
        x_values_new <- c(x_values_new, n)
        x_new <- c(x_new, as.character(x_new[length(x_new)]))
        y_new <- c(y_new, NA)
      } else {
        for (m in 1:length(indices)) {
          x_values_new <- c(x_values_new, x_values[indices[m]])
          x_new <- c(x_new, as.character(x[indices[m]]))
          y_new <- c(y_new, y[indices[m]])

        }
      }
    }

  } else {

    # insert NA's for missing x_values
    x_values_new <- min(x_values):max(x_values)
    x_new <- c()
    y_new <- c()
    indices_with_doubles <- c()
    number_doubles <- c()
    for (i in x_values_new) {
      indices <- which(x_values==i)
      if (length(indices)==0) { # no values, insert NA
        x_new <- c(x_new, as.character(x_new[length(x_new)]))
        y_new <- c(y_new, NA)
      } else if (length(indices)==1) { # one single values, use this
        x_new <- c(x_new, as.character(x[indices]))
        y_new <- c(y_new, y[indices])
      } else { # more values, take the mean mean
        x_new <- c(x_new, as.character(x[indices][1]))
        y_new <- c(y_new, mean(y[indices]))
        number_doubles <- c(number_doubles, length(indices))
        indices_with_doubles <- c(indices_with_doubles, indices)
      }
    }
    # plot(x_values, y, type ="p", col="blue")
    # lines(x_values_new, y_new, type ="p", col="red")

    # try to shift double value to previous or next NA
    if (length(number_doubles)>0) {
      cum_number <- 0
      for (i in 1:length(number_doubles)) {
        number_doubles_i <- number_doubles[i]
        doubles_i <- indices_with_doubles[(cum_number+1):(cum_number+number_doubles_i)]
        x_value_doubles_i <- x_values[doubles_i][1] #always the same
        y_doubles_i <- y[doubles_i]
        remaining_doubles <- c()
        cum_number <- cum_number+number_doubles_i
        #print(x_value_doubles_i)
        #print(y_doubles_i)
        indices_new <- which(x_values_new==x_value_doubles_i)
        for (j in 1:length(y_doubles_i)) {
          if (number_doubles_i>1 && j==1 && indices_new>1 && is.na(y_new[indices_new-1]) && x_new[indices_new-1]==x_new[indices_new]) {
            y_new[indices_new-1] <- y_doubles_i[j]
            number_doubles_i <- number_doubles_i-1
          } else if (number_doubles_i>1 && j==length(y_doubles_i) && indices_new<length(y_new) && is.na(y_new[indices_new+1]) && x_new[indices_new+1]==x_new[indices_new]) {
            y_new[indices_new+1] <- y_doubles_i[j]
            number_doubles_i <- number_doubles_i-1
          } else {
            remaining_doubles <- c(remaining_doubles, y_doubles_i[j])
          }
        }
        y_new[indices_new] <- mean(remaining_doubles) # aggregate remaining identical x_values
      }
    }

  }


  # plot(x_values, y, type ="p", col="blue")
  # lines(x_values_new, y_new, type ="p", col="green")

  list(
    x = factor(x_new, levels=levels(x)),
    x_values = x_values_new,
    y = y_new
  )

}


