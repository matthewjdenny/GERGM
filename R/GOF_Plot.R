# Function for plotting goodness of fit statistics in a boxplot
Gof_Plot <- function(GERGM_Object){
  #define colors
  UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
  UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)

  #input is the GERGM_Object, we are going to normalize all statistics
  temp <- (GERGM_Object@MCMC_output$Statistics - GERGM_Object@stats[2, ])/apply(GERGM_Object@MCMC_output$Statistics,1,sd)

  boxplot(temp, medcol = UMASS_RED,
          xlab = "Network Statistic",
          ylab = "Standardized, Normalized Values",
          main = "Blue = Observed Statistic, Red = Simulated Mean")
  zero_line <- rep(0,length(GERGM_Object@stats[2, ]))
  zero_plot <- rbind(zero_line,zero_line)
  boxplot(zero_plot, add =T, medcol=UMASS_BLUE, names = F)
}

Comparison_GOF_Plot <- function(MH.obj, Gibbs.obj, gergm.obj){
  par(mfrow = c(2,3))
  violins(data.frame(Gibbs.obj$Statistics[,1], MH.obj$Statistics[,1]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Out-2-Stars")
  abline(h = gergm.obj@stats[2,1], col = "red", lty = 2, lwd = 2)

  violins(data.frame(Gibbs.obj$Statistics[,2], MH.obj$Statistics[,2]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "In-2-Stars")
  abline(h = gergm.obj@stats[2,2], col = "red", lty = 2, lwd = 2)

  violins(data.frame(Gibbs.obj$Statistics[,3], MH.obj$Statistics[,3]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Cyclic Triads")
  abline(h = gergm.obj@stats[2,3], col = "red", lty = 2, lwd = 2)

  violins(data.frame(Gibbs.obj$Statistics[,4], MH.obj$Statistics[,4]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Reciprocity")
  abline(h = gergm.obj@stats[2,4], col = "red", lty = 2, lwd = 2)

  violins(data.frame(Gibbs.obj$Statistics[,5], MH.obj$Statistics[,5]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Transitive Triads")
  abline(h = gergm.obj@stats[2,5], col = "red", lty = 2, lwd = 2)

  violins(data.frame(Gibbs.obj$Statistics[,6], MH.obj$Statistics[,6]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Total Edge Weight")
  abline(h = gergm.obj@stats[2,6], col = "red", lty = 2, lwd = 2)
}
