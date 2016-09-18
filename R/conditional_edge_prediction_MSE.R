#' @title A Function to calcualte Mean Edgewise MSE to evaluate edge predictions.
#' @description Calculates mean edgewise MSE for predicted edge values.
#'
#' @param edge_prediction_results A list object returned by the
#' `conditional_edge_prediction()` function.
#' @return A list of MSE's.
#' @export
conditional_edge_prediction_MSE <- function(
  edge_prediction_results) {

  cur <- edge_prediction_results$max_ent_edge_predictions[,,1]
  # calculate average MSE
  dimensions <- ncol(cur)
  mses <- matrix(NA, ncol = dimensions, nrow = dimensions)
  mses_max_ent <- mses
  obv <- edge_prediction_results$observed_network

  for (i in 1:dimensions) {
    for (j in 1:dimensions) {
      if (i != j) {
        cur <- edge_prediction_results$edge_predictions[i,j,]
        cur2 <- edge_prediction_results$max_ent_edge_predictions[i,j,]
        mses[i,j] <- mean((cur - obv[i,j])^2)
        mses_max_ent[i,j] <- mean((cur2 - obv[i,j])^2)
      }
    }
  }

  mse_predicted <- mean(c(mses), na.rm = TRUE)
  mse_max_ent <- mean(c(mses_max_ent), na.rm = TRUE)

  mse1 <- round(mse_predicted,4)
  mse2 <- round(mse_max_ent,4)
  cat("Mean MSE for Predicted Edge Values:", mse1,"\n")
  cat("Mean MSE for Max Ent Predicted Edge Values:", mse2,"\n")

  reduction <- (1 - mse1/mse2)
  cat("This represents a",round(reduction*100,2),
      "percent reduction in the average edgewise MSE when using the GERGM model.\n")

  return(list(model_prediction_average_MSE = mse1,
              max_ent_prediction_average_MSE = mse2))
}

