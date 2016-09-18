draw_max_entropy_conditional_uniform_corr <- function(
  x,
  i,
  j,
  n = 1000,
  increment = .001) {
  # x is a valid correlation matrix
  # i j are the row and column indices
  # increment is the increment by which we search for bounds
  # note, the larger "increment" the faster, but the bounds are less precise
  current_value <- x[i,j]
  upper <- current_value
  still_positive_definite <- matrixcalc::is.positive.definite(x)
  x_new <- x
  while (still_positive_definite & upper <= 1+increment) {
    upper <- upper + increment
    x_new <- x
    x_new[i,j] <- upper
    x_new[j,i] <- upper
    still_positive_definite <- matrixcalc::is.positive.definite(x_new)
  }
  upper <- upper - increment
  lower <- current_value
  still_positive_definite <- matrixcalc::is.positive.definite(x)
  x_new <- x
  while(still_positive_definite & lower >= -1-increment){
    lower <- lower - increment
    x_new <- x
    x_new[i,j] <-  lower
    x_new[j,i] <- lower
    still_positive_definite <- matrixcalc::is.positive.definite(x_new)
  }
  lower <- lower + increment
  # Returns n simulated correlation coefficients that adhere to the positive definiteness constraint
  cat("Max Ent Lower Bound:",lower,"Upper Bound:",upper,"\n")
  runif(n,lower,upper)
}
