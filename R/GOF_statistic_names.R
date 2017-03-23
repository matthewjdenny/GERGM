GOF_statistic_names <- function(GERGM_Object){

  if (GERGM_Object@include_diagonal) {
    if (!GERGM_Object@directed_network) {
      column_names <- c("Two\nStars",
                        "Transitive\nTriads",
                        "Network\nDensity",
                        "Diagonal\n",
                        "Average\nIntensity")
    } else {
      if (length(which(colnames(GERGM_Object@stats)== "modularity")) == 1) {
        column_names <- c("Out\n2-Stars",
                          "In\n2-Stars",
                          "Cyclic\nTriads",
                          "Mutual\nDyads",
                          "Transitive\nTriads",
                          "Network\nDensity",
                          "Diagonal\n",
                          "Average\nIntensity",
                          "Modularity\n")
      } else {
        column_names <- c("Out\n2-Stars",
                          "In\n2-Stars",
                          "Cyclic\nTriads",
                          "Mutual\nDyads",
                          "Transitive\nTriads",
                          "Network\nDensity",
                          "Diagonal\n",
                          "Average\nIntensity")
      }
    }
  } else {
    if (!GERGM_Object@directed_network) {
      column_names <- c("Two\nStars",
                        "Transitive\nTriads",
                        "Network\nDensity",
                        "Average\nIntensity")
    } else {
      if (length(which(colnames(GERGM_Object@stats)== "modularity")) == 1) {
        column_names <- c("Out\n2-Stars",
                          "In\n2-Stars",
                          "Cyclic\nTriads",
                          "Mutual\nDyads",
                          "Transitive\nTriads",
                          "Network\nDensity",
                          "Average\nIntensity",
                          "Modularity\n")
      } else {
        column_names <- c("Out\n2-Stars",
                          "In\n2-Stars",
                          "Cyclic\nTriads",
                          "Mutual\nDyads",
                          "Transitive\nTriads",
                          "Network\nDensity",
                          "Average\nIntensity")
      }
    }
  }


  # store the number of base stats
  base_length <- length(column_names)

  # determine if additional non-base covariates were provided
  non_base_stats <- which(GERGM_Object@non_base_statistic_indicator == 1)
  if (length(non_base_stats) > 0) {

    # get the full theta names that are non-base
    subset_inds <- GERGM_Object@statistic_auxiliary_data$specified_statistic_indexes_in_full_statistics
    subset_inds <- subset_inds[non_base_stats]
    names <- GERGM_Object@full_theta_names[subset_inds]
    types <- column_names[GERGM_Object@stats_to_use[non_base_stats]]

    # now append on the subset to each type
    for (i in 1:length(types)) {
      tpe <- types[i]
      category <- stringr::str_split(names[i], "_")[[1]]
      column_names <- c(column_names,
                        paste(tpe,"\n",category[2]," = ",category[3],sep = ""))
    }
  }

  # add on additional row to each of the base stats to make things line up
  # nicely
  for (j in 1:base_length) {
    column_names[j] <- paste(column_names[j],"\n",sep = "")
  }

  return(column_names)
}
