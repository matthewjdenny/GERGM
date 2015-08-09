#' Generate trace plot of network density from a GERGM object.
#'
#' @param GERGM_Object The object returned by the estimation procedure using the GERGM function.
#' @return A trace plot of network density.
#' @export
Trace_Plot <- function(GERGM_Object){
  UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
  UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)
  stats <- GERGM_Object@MCMC_output$Statistics
  indexes <- 1:length(stats[,1])
  stats <- cbind(stats,indexes)
  p <- ggplot2::ggplot(stats, ggplot2::aes(x=indexes, y=edgeweight ))
  p  <- p + ggplot2::geom_line(color = UMASS_BLUE) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab("Network Density") +
    ggplot2::ggtitle("Trace Plot of Density for Simulated Networks")
  print(p)
}
