# convert a formula object to statistics and network of interest (formerly weighted.ergm.data)
Parse_Formula_Object <- function(formula,
                                 possible.stats,
                                 theta = NULL,
                                 alpha = NULL) {

  num_stats <- length(possible.stats)

  # parse the formula
  if (class(formula) != "formula") {
    stop("'formula' must be a formula object.")
  }
  lhs <- deparse(formula[[2]])  # name of the response variable
  net <- eval(parse(text = lhs))  # get the actual response data

  # check that the response is a matrix
  if (class(net) != "matrix") {
    stop("Response must be a matrix object.")
  }

  rhs <- paste0(deparse(formula[[3]]), collapse = "")  # rhs of formula
  rhs <- gsub("\\s+", " ", rhs)  # get rid of redundant spaces
  rhs <- strsplit(rhs, " \\+ ")[[1]]  # parse separate formula elements

  # check that the names of the statistics match those that are possible
  possible <- 1:length(rhs)
  actual <- which(rhs %in% possible.stats)
  if (length(possible) != length(actual)) {
    stop(paste("the specified variable", rhs[setdiff(possible, actual)],
               "is not an available statistic.", sep = " "))
  }

  # if alpha is NULL, assume all ones
  if (is.null(alpha) == TRUE) {
    alpha <- rep(1, length(rhs))
  }

  # if theta is NULL, assume all ones
  if (is.null(theta) == TRUE) {
    theta <- rep(1, length(rhs))
  }
  # check that alpha and the number of statistics are equal
  if (length(rhs) != length(alpha)) {
    stop("'alpha' must be the same length as the number of statistics")
  }

  # check that theta and the number of statistics are equal
  if (length(rhs) != length(theta)) {
    stop("'theta' must be the same length as the number of statistics")
  }
  stat.indx <- which(possible.stats %in% rhs)
  statistics <- rep(0, num_stats)
  statistics[stat.indx] <- 1
  alphas <- rep(1, num_stats)
  thetas <- rep(0, num_stats)
  for (i in 1:length(rhs)) {
    alphas[which(rhs[i] == possible.stats)] <- alpha[i]
    thetas[which(rhs[i] == possible.stats)] <- theta[i]
  }
  return(list(net = net, statistics = statistics,
              alphas = alphas, thetas = thetas))
}
