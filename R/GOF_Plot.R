# Function for plotting goodness of fit statistics in a boxplot
Gof_Plot <- function(GERGM_Object){
  # simulation.obj is a list object resulting from simulations
  # gergm.obj is an object of class "gergm"
  boxplot(GERGM_Object@MCMC_output$Statistics, medcol = "red")
  boxplot(rbind(GERGM_Object@stats[2, ], GERGM_Object@stats[2, ]), add =T,
          medcol="blue", names = F)
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
