calculate_standard_errors <- function(hessian,
                                      GERGM_Object){

  # calculate standard errors inside of try statement to better handle cases
  # where we end up with degeneracy
  std_errors <- try(sqrt(diag(solve(-hessian))))

  # if we get an error, set all standard errors to NaN.
  if (class(std_errors) == "try-error") {
    std_errors <- matrix(NaN,nrow = nrow(hessian), ncol = nrow(hessian))
  }

  if (is.null(GERGM_Object@hessians)) {
    GERGM_Object@hessians <- list(hessian)
  } else {
    GERGM_Object@hessians <- append(GERGM_Object@hessians,hessian)
  }

  if (is.null(GERGM_Object@standard_errors)) {
    GERGM_Object@standard_errors <- list(std_errors)
  } else {
    GERGM_Object@standard_errors <- append(GERGM_Object@standard_errors,
                                           std_errors)
  }

  return(list(GERGM_Object = GERGM_Object,
              std_errors = std_errors))
}
