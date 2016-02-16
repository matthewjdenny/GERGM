#' Generate trace plot of network density from a GERGM object.
#'
#' @param GERGM_Object The object returned by the estimation procedure using the
#' GERGM function.
#' @return A trace plot of network density.
#' @export
Trace_Plot <- function(GERGM_Object){
  UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
  edges <- NULL
  stats <- GERGM_Object@simulated_statistics_for_GOF

  # get the total number of indices
  indexes <- 1:length(stats[,1])

  iteration_number <- seq(1,GERGM_Object@number_of_simulations,
                          length.out = length(indexes))
  nr <- nrow(GERGM_Object@network)
  normalizer <- nr * (nr - 1)
  stats <- cbind(stats/normalizer,iteration_number)
  actual_density <- GERGM_Object@stats[2,6]/normalizer

  p <- ggplot2::ggplot(stats, ggplot2::aes(x = iteration_number, y = edges))
  p  <- p + ggplot2::geom_line(color = UMASS_BLUE) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab("Constrained Network Density") +
    ggplot2::ggtitle("Trace Plot of Density for Simulated Networks") +
    ggplot2::geom_abline(intercept = actual_density, slope = 0) +
    ggplot2::scale_x_continuous(labels = scales::comma)
  print(p)
}
