#' Generate Goodness Of Fit plot from a GERGM object.
#'
#' @param GERGM_Object The object returned by the estimation procedure using the
#' GERGM function.
#' @param column_names Optional argument allowing the user to specify names for
#' statistics to be plotted.
#' @param modularity_group_memberships Optional numeric vector of node group
#' memberships indexed from 1, which will be used to calculate network
#' modularities.
#' @param return_GERGM_Object Optional argument to return the GERGM Object that
#' was passed in, but now with additional GOF statistics such as modularity
#' included in the `@simulated_statistics_for_GOF` and `@additional_stats`
#' fields
#' @param observed_support logical indicating whether GOF plots should use
#' observed support. Defaults to FALSE.
#' @param ... Additional Arguments can be passed in. Included for eventual
#' compatibility with XERGM package.
#' @return A set of box plots where of simulated network statistics centered at
#' the observed value for those statistics and normalized by their standard
#' deviation. This aids in interpretation as the y-axis can be interpreted as
#' the number of simulated-sample standard deviations above or below the
#' observed statistic. Optionally also returns a GERGM object with updated
#' statistics.
#' @export
GOF <- function(GERGM_Object,
                column_names = NULL,
                modularity_group_memberships = NULL,
                return_GERGM_Object = FALSE,
                observed_support = FALSE,
                ...){
  #define colors
  UMASS_BLUE <- rgb(51,51,153,155,maxColorValue = 255)
  UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)

  GERGM_Object <- calculate_additional_GOF_statistics(
    GERGM_Object,
    modularity_group_memberships)

  # deal with case where wer are not including the diagonal
  if (!GERGM_Object@include_diagonal) {
    rem <- which(colnames(GERGM_Object@simulated_statistics_for_GOF) == "diagonal")
    if (length(rem) > 0) {
      GERGM_Object@simulated_statistics_for_GOF <- GERGM_Object@simulated_statistics_for_GOF[,-rem]
      GERGM_Object@stats <- GERGM_Object@stats[,-rem]
    }
  }

  if (observed_support) {
    temp <- GERGM_Object@simulated_statistics_for_GOF
    temp3 <- GERGM_Object@MCMC_output$Networks
    temp4 <- GERGM_Object@MCMC_output$Statistics
    # recalculate statistics
    for(l in 1:nrow(temp4)) {
      temp4[l,] <- calculate_h_statistics(
        GERGM_Object,
        GERGM_Object@statistic_auxiliary_data,
        all_weights_are_one = FALSE,
        calculate_all_statistics = TRUE,
        use_constrained_network = FALSE,
        network = temp3[,,l])
    }
    temp[,1:ncol(temp4)] <- temp4
    stats <- GERGM_Object@stats[1,]

  } else {
    temp <- GERGM_Object@simulated_statistics_for_GOF
    stats <- GERGM_Object@stats[2,]
  }
  # thin chain for accurate t-test stats
  temp <- Thin_Statistic_Samples(temp)

  cat("\nStatistics used in comparison:\n")
  print(stats)
  cat("\nMean statistic values from simulated networks:\n")
  print(colMeans(temp))
  cat("\n")
  # calculate t-test statistics
  t_stats <- stats
  for (i in 1:length(stats)) {
    if (length(unique(temp[,i])) > 1) {
      t_stats[i] <- as.numeric(t.test(temp[,i], mu = stats[i])$statistic)
    }
  }



  cat("t-statistics for difference of simulated mean from observed statistic...\n")
  print(t_stats)

  for (i in 1:length(stats)) {
    mn <- mean(temp[,i])
    ste <- sd(temp[,i])
    t_stats[i] <- (stats[i] - mn)/ste
  }
  cat("\n")
  cat("t-statistics for test of whether observed statistic is an outlier with respect to simulated statistics...\n")
  print(t_stats)
  cat("\n")

  # check to see if user provided column names, if not, then generate them
  if (is.null(column_names)) {
    column_names <- GOF_statistic_names(GERGM_Object)
  } else {
    print("Original statistic names:")
    print(GOF_statistic_names(GERGM_Object))
    print("User specified statistic names:")
    print(column_names)
  }

  if (GERGM_Object@simulation_only) {
    if (!GERGM_Object@directed_network) {
      #input is the GERGM_Object, we are going to normalize all statistics

      par(mar = c(2.5,4,2.5,2))
      layout(matrix(c(1,2), 1, 2, byrow = TRUE),
             widths=c(3,1), heights=c(1))
      boxplot(temp, medcol = UMASS_BLUE,
              xlab = "Network Statistic",
              ylab = "Statistic Values",
              main = "Simulated Network Statistics")

      par(mar = c(2.5,2,2.5,2))
      sim <- GERGM_Object@additional_stats$simulated_in_degrees
      plot (density(sim),
            main = "Degree Distribution",
            col = UMASS_RED,
            lwd = 3,
            xlim = c(min(c(sim)), max(c(sim))),
            ylim = c(0, max(density(sim)$y)))

    } else {
      #input is the GERGM_Object, we are going to normalize all statistics

      par(mar = c(2.5,4,2.5,2))
      layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE),
             widths = c(3,1), heights=c(1,1))
      boxplot(temp, medcol = UMASS_BLUE,
              xlab = "Network Statistic",
              ylab = "Statistic Values",
              main = "Simulated Network Statistics")

      par(mar = c(2.5,2,2.5,2))
      sim <- GERGM_Object@additional_stats$simulated_in_degrees
      plot (density(sim),
            main = "In-Degree Distribution",
            col = UMASS_BLUE,
            lwd = 3,
            xlim = c(min(c(sim)), max(c(sim))),
            ylim = c(0, max(density(sim)$y)))


      sim <- GERGM_Object@additional_stats$simulated_out_degrees
      plot (density(sim),
            main = "Out-Degree Distribution",
            col = UMASS_BLUE,
            lwd = 3,
            xlim = c(min(c(sim)), max(c(sim))),
            ylim = c(0, max(density(sim)$y)))
    }


  } else {
    # if we are not just working with simulated networks, but with a GERGM fit:
    if (!GERGM_Object@directed_network) {
      #input is the GERGM_Object, we are going to normalize all statistics

      temp2 <- apply(temp,2,sd)
      for (i in 1:ncol(temp)) {
        temp[,i] <- temp[,i] - stats[i]
        temp[,i] <- temp[,i]/temp2[i]
      }

      #now we are only dealing with ttriads and twostars since this is an undirected network
      colnames(temp) <- column_names

      par(mar = c(2.5,4,2.5,2))
      layout(matrix(c(1,2), 1, 2, byrow = TRUE),
             widths=c(3,1), heights=c(1))
      boxplot(temp, medcol = UMASS_RED,
              xlab = "",
              ylab = "Normalized Statistic Values",
              main = "Blue = Observed Statistic, Red = Simulated Mean",
              xaxt = "n",
              las = 1)

      zero_line <- rep(0,length(stats))
      zero_plot <- rbind(zero_line,zero_line)

      axis(side = 1,at = 1:ncol(temp),colnames(temp),tick = FALSE, line = 1)
      boxplot(zero_plot, add = T, medcol = UMASS_BLUE, names = F, axes = F)

      par(mar = c(2.5,2,2.5,2))
      sim <- GERGM_Object@additional_stats$simulated_in_degrees
      obs <- GERGM_Object@additional_stats$observed_in_degrees
      plot (density(sim),
            main = "Degree Distribution",
            col = UMASS_RED,
            lwd = 3,
            xlim = c(min(c(sim,obs)), max(c(sim,obs))),
            ylim = c(0, max(max(density(sim)$y),max(density(obs)$y))),
            xlab = "")
      lines (density(obs),
             col = UMASS_BLUE,
             lwd = 3)


    } else {
      #input is the GERGM_Object, we are going to normalize all statistics

      colnames(temp) <- column_names

      temp2 <- apply(temp,2,sd)
      for(i in 1:ncol(temp)){
        temp[,i] <- temp[,i] - stats[i]
        temp[,i] <- temp[,i]/temp2[i]
      }
      par(mar = c(2.5,4,2.5,2))
      layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE),
             widths=c(3,1), heights=c(1,1))
      boxplot(temp, medcol = UMASS_RED,
              xlab = "",
              ylab = "Normalized Statistic Values",
              main = "Blue = Observed Statistic, Red = Simulated Mean",
              xaxt = "n",
              las = 1)
      zero_line <- rep(0,length(stats))
      zero_plot <- rbind(zero_line,zero_line)
      axis(side = 1,at = 1:ncol(temp),colnames(temp),tick = FALSE, line = 1)
      boxplot(zero_plot, add = T, medcol = UMASS_BLUE, names = F, axes = F)

      par(mar = c(2.5,2,2.5,2))
      sim <- GERGM_Object@additional_stats$simulated_in_degrees
      obs <- GERGM_Object@additional_stats$observed_in_degrees
      plot (density(sim),
            main = "In-Degree Distribution",
            col = UMASS_RED,
            lwd = 3,
            xlim = c(min(c(sim,obs)), max(c(sim,obs))),
            ylim = c(0, max(max(density(sim)$y),max(density(obs)$y))),
            xlab = "")
      lines (density(obs),
             col = UMASS_BLUE,
             lwd = 3)

      sim <- GERGM_Object@additional_stats$simulated_out_degrees
      obs <- GERGM_Object@additional_stats$observed_out_degrees
      plot (density(sim),
            main = "Out-Degree Distribution",
            col = UMASS_RED,
            lwd = 3,
            xlim = c(min(c(sim,obs)), max(c(sim,obs))),
            ylim = c(0, max(max(density(sim)$y),max(density(obs)$y))),
            xlab = "")
      lines (density(obs),
             col = UMASS_BLUE,
             lwd = 3)
    }
  }
  par(mfrow = c(1,1))
  if (return_GERGM_Object) {
    return(GERGM_Object)
  }
}
