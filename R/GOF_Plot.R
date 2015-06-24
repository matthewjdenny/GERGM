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

Estimate_Plot <- function(GERGM_Object){
  #define colors
  UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
  UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)
  GERGM_Object@theta.coef
  modelFrame <- data.frame(Variable = colnames(GERGM_Object@theta.coef) ,
                            Coefficient = GERGM_Object@theta.coef[,1],
                            SE = GERGM_Object@theta.coef[,2],
                            Model = "Theta Estimates"
  )
  data <- data.frame(modelFrame)

  # Plot
  zp1 <- ggplot2::ggplot(data, ggplot2::aes(colour = Model)) +
    ggplot2::scale_color_manual(values = rgb(51,51,153,255,maxColorValue = 255))
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
  zp1 <- zp1  + ggplot2::theme_bw() +
                ggplot2::coord_flip() +
                ggplot2::theme(legend.position="none")
  print(zp1)
}

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

