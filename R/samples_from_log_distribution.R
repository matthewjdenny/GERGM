samples_from_log_distribution <- function(log_sampling_weights,
                                          values_to_be_sampled,
                                          num_samples) {

  # allocate a vector to hold indices of values we actually want to sample.
  indices <- rep(0,num_samples)

  # populate indices with log space sampler
  for(i in 1:num_samples) {
    indices[i] <- log_space_multinomial_sampler(log_sampling_weights,
                                             runif(1))
  }

  sampled_values <- values_to_be_sampled[indices]

  return(sampled_values)
}

