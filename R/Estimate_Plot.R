#' Generate parameter estimate plot with 95 percent CI's from a GERGM object.
#'
#' @param GERGM_Object The object returned by the estimation procedure using the
#' GERGM function.
#' @param normalize_coefficients Defaults to FALSE, if TRUE then parameter
#' estimates will be converted be deivided by their standard deviations with
#' and displayed with 95 percent confidence intervals. These coefficients will
#' no longer be comparable, but make graphical interpretation of significance
#' and sign easier.
#' @param coefficients_to_plot An optional argument indicating which kind of
#' parameters to plot. Can be one of "both","covariate", or "structural". Useful
#' for creating separate parameter plots for covariates and structural
#' parameters when these parameters are on very different scales.
#' @param coefficient_names Defualts to NULL. Can be a string vector of names
#' for coefficients to be used in making publication quality plots.
#' @param leave_out_coefficients Defaults to NULL. Can be a string vector of
#' coefficient names as they appear in the plot. These coefficients will be
#' removed from the final plot. Useful if the intercept term is much larger in
#' magnitude than other estimates, and the user wishes to clarify the other
#' parameter estimates without normalizing.
#' given
#' @return A parameter estimate plot.
#' @export
Estimate_Plot <- function(
  GERGM_Object,
  normalize_coefficients = FALSE,
  coefficients_to_plot = c("both","covariate","structural"),
  coefficient_names = NULL,
  leave_out_coefficients = NULL
  ){

  coefficients_to_plot <- coefficients_to_plot[1]

  #define colors
  UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
  UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)
  Model <- Variable <- Coefficient <- SE <- NULL

  if (coefficients_to_plot == "both") {
    # make sure that we use rows as estimates and se are in a two row matrix
    modelFrame <- data.frame(Variable = colnames(GERGM_Object@theta.coef) ,
                             Coefficient = as.numeric(GERGM_Object@theta.coef[1,]),
                             SE = as.numeric(GERGM_Object@theta.coef[2,]),
                             Model = "Theta Estimates"
    )
    data <- data.frame(modelFrame)

    #now add in lambda estimates
    if(length(GERGM_Object@lambda.coef) > 1){
      temp1 <- as.numeric(GERGM_Object@lambda.coef[1,])
      temp1 <- temp1[1:(length(temp1)-1)]
      temp2 <- as.numeric(GERGM_Object@lambda.coef[2,])
      temp2 <- temp2[1:(length(temp2)-1)]
      modelFrame2 <- data.frame(Variable = dimnames(GERGM_Object@data_transformation)[[3]] ,
                                Coefficient = temp1,
                                SE = temp2,
                                Model = "Lambda Estimates"
      )
      data2 <- data.frame(modelFrame2)
      data <- rbind(data,data2)
    }
  } else if (coefficients_to_plot == "covariate") {
    #now add in lambda estimates
    if(length(GERGM_Object@lambda.coef) > 1){
      temp1 <- as.numeric(GERGM_Object@lambda.coef[1,])
      temp1 <- temp1[1:(length(temp1)-1)]
      temp2 <- as.numeric(GERGM_Object@lambda.coef[2,])
      temp2 <- temp2[1:(length(temp2)-1)]
      modelFrame2 <- data.frame(Variable = dimnames(GERGM_Object@data_transformation)[[3]] ,
                                Coefficient = temp1,
                                SE = temp2,
                                Model = "Lambda Estimates"
      )
      data <- data.frame(modelFrame2)
    }
  } else if (coefficients_to_plot == "structural") {
    # make sure that we use rows as estimates and se are in a two row matrix
    modelFrame <- data.frame(Variable = colnames(GERGM_Object@theta.coef) ,
                             Coefficient = as.numeric(GERGM_Object@theta.coef[1,]),
                             SE = as.numeric(GERGM_Object@theta.coef[2,]),
                             Model = "Theta Estimates"
    )
    data <- data.frame(modelFrame)
  }

  # rename coefficients if necessary
  if (!is.null(coefficient_names)) {
    if (length(coefficient_names) != nrow(data)) {
      stop("coefficient_names must be the same length as the number of covariates in the plot. Try setting coefficient_names = NULL and counting the number of coefficients to check that you have provided the right number.")
    }
    cat("Replacing:\n")
    print(data$Variable)
    with("With..\n")
    print(coefficient_names)
    data$Variable <- coefficient_names
  }

  #remove variables
  if (!is.null(leave_out_coefficients)) {
    for (j in 1:length(leave_out_coefficients)) {
      remove <- which(data$Variable == leave_out_coefficients[j])
      if (length(remove) == 1) {
        data <- data[-remove,]
        cat("Removing variable",leave_out_coefficients[j],"\n")
      } else if (length(remove) > 1) {
        cat("Your argument matches more than one variable. Please respecify.\n")
      } else {
        cat ("There is no variable",leave_out_coefficients[j],
             ".Please respecify.\n")
      }
    }
  }


  # standardize coefficients
  if (normalize_coefficients) {
    for (i in 1:nrow(data)) {
      data$Coefficient[i] <- data$Coefficient[i]/data$SE[i]
      data$SE[i] <- 1
    }
  }

  # Plot
  if (length(GERGM_Object@lambda.coef[,1]) > 0 &
     (coefficients_to_plot == "covariate" | coefficients_to_plot == "both")) {
    zp1 <- ggplot2::ggplot(data, ggplot2::aes(colour = Model)) +
      ggplot2::scale_color_manual(values = c(UMASS_BLUE,UMASS_RED))
  } else {
    zp1 <- ggplot2::ggplot(data, ggplot2::aes(colour = Model)) +
      ggplot2::scale_color_manual(values = UMASS_BLUE)
  }
  zp1 <- zp1 + ggplot2::geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
  zp1 <- zp1 + ggplot2::geom_linerange( ggplot2::aes(x = Variable,
                  ymin = Coefficient - SE*(-qnorm((1-0.9)/2)),
                  ymax = Coefficient + SE*(-qnorm((1-0.9)/2))),
                  lwd = 1,
                  position = ggplot2::position_dodge(width = 1/2))
  zp1 <- zp1 + ggplot2::geom_pointrange(ggplot2::aes(x = Variable,
                  y = Coefficient,
                  ymin = Coefficient - SE*(-qnorm((1-0.95)/2)),
                  ymax = Coefficient + SE*(-qnorm((1-0.95)/2))),
                  lwd = 1/2,
                  position = ggplot2::position_dodge(width = 1/2),
                  shape = 21, fill = "WHITE")
  if(normalize_coefficients){
    zp1 <- zp1  + ggplot2::theme_bw() +
      ggplot2::coord_flip() +
      ggplot2::theme(legend.position="none") +
      ggplot2::ylab("Normalized Coefficient")
  }else{
    zp1 <- zp1  + ggplot2::theme_bw() +
      ggplot2::coord_flip() +
      ggplot2::theme(legend.position="none")
  }
  print(zp1)
}
