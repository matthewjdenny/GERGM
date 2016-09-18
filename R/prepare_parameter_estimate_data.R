prepare_parameter_estimate_data <- function(GERGM_Object,
                                            normalize_coefficients,
                                            coefficients_to_plot,
                                            coefficient_names,
                                            leave_out_coefficients,
                                            model){

  if (coefficients_to_plot == "both") {
    # make sure that we use rows as estimates and se are in a two row matrix
    modelFrame <- data.frame(Variable = colnames(GERGM_Object@theta.coef) ,
                             Coefficient = as.numeric(GERGM_Object@theta.coef[1,]),
                             SE = as.numeric(GERGM_Object@theta.coef[2,]),
                             Coefficient_Type = "Theta Estimates",
                             Model = model
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
                                Coefficient_Type = "Lambda Estimates",
                                Model = model
      )
      data2 <- data.frame(modelFrame2)
      data <- rbind(data,data2)

      # add in separate dispersion parameter estimate
      temp1 <- as.numeric(GERGM_Object@lambda.coef[1,])
      temp1 <- temp1[length(temp1)]
      temp2 <- as.numeric(GERGM_Object@lambda.coef[2,])
      temp2 <- temp2[length(temp2)]
      modelFrame2 <- data.frame(Variable = "Dispersion Parameter" ,
                                Coefficient = temp1,
                                SE = temp2,
                                Coefficient_Type = "Dispersion Parameter",
                                Model = model
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
                                Coefficient_Type = "Lambda Estimates",
                                Model = model
      )
      data <- data.frame(modelFrame2)

      # add in separate dispersion parameter estimate
      temp1 <- as.numeric(GERGM_Object@lambda.coef[1,])
      temp1 <- temp1[length(temp1)]
      temp2 <- as.numeric(GERGM_Object@lambda.coef[2,])
      temp2 <- temp2[length(temp2)]
      modelFrame2 <- data.frame(Variable = "Dispersion Parameter",
                                Coefficient = temp1,
                                SE = temp2,
                                Coefficient_Type = "Dispersion Parameter",
                                Model = model
      )
      data2 <- data.frame(modelFrame2)
      data <- rbind(data,data2)
    }
  } else if (coefficients_to_plot == "structural") {
    # make sure that we use rows as estimates and se are in a two row matrix
    modelFrame <- data.frame(Variable = colnames(GERGM_Object@theta.coef) ,
                             Coefficient = as.numeric(GERGM_Object@theta.coef[1,]),
                             SE = as.numeric(GERGM_Object@theta.coef[2,]),
                             Coefficient_Type = "Theta Estimates",
                             Model = model
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
    cat("With..\n")
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

  return(data)
}
