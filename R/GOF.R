#' Generate Goodness Of Fit plot from a GERGM object.
#'
#' @param GERGM_Object The object returned by the estimation procedure using the GERGM function.
#' @return A GOF plot.
#' @export
GOF <- function(GERGM_Object){
  #define colors
  UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
  UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)

  if(GERGM_Object@simulation_only){
    temp <- GERGM_Object@MCMC_output$Statistics

    boxplot(temp, medcol = UMASS_BLUE,
            xlab = "Network Statistic",
            ylab = "Statistic Values",
            main = "Simulated Network Statistics")
  }else{
    #input is the GERGM_Object, we are going to normalize all statistics
    temp <- GERGM_Object@MCMC_output$Statistics

    temp2 <- apply(GERGM_Object@MCMC_output$Statistics,2,sd)
    for(i in 1:ncol(temp)){
      temp[,i] <- temp[,i] - GERGM_Object@stats[2,i]
      temp[,i] <- temp[,i]/temp2[i]
    }

    boxplot(temp, medcol = UMASS_RED,
            xlab = "Network Statistic",
            ylab = "Standardized, Normalized Values",
            main = "Blue = Observed Statistic, Red = Simulated Mean")
    zero_line <- rep(0,length(GERGM_Object@stats[2, ]))
    zero_plot <- rbind(zero_line,zero_line)
    boxplot(zero_plot, add =T, medcol=UMASS_BLUE, names = F)
  }
}
