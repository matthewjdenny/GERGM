// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]

#include <RcppParallel.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using std::pow;
using std::exp;
using std::log;
using std::sqrt;

namespace wobj {

  using std::pow;
  using std::exp;
  using std::log;
  using std::sqrt;

  // Function to calculate the number of out 2-stars
  double Out2Star(arma::mat net,
                  arma::mat triples,
                  double alpha,
                  int together) {

      int number_of_triples = triples.n_rows;
      double st1 = 0;
      double st2 = 0;
      double st3 = 0;

      double to_return = 0;
      if (together == 1) {
          for (int i = 0; i < number_of_triples; ++i) {
              st1 += net(triples(i, 0), triples(i, 1)) * net(triples(i, 0),
                         triples(i, 2));
              st2 += net(triples(i, 1), triples(i, 0)) * net(triples(i, 1),
                         triples(i, 2));
              st3 += net(triples(i, 2), triples(i, 0)) * net(triples(i, 2),
                         triples(i, 1));
          }
          to_return = pow((st1 + st2 + st3), alpha);
      } else {
          for (int i = 0; i < number_of_triples; ++i) {
              st1 += pow(net(triples(i, 0), triples(i, 1)), alpha) *
                  pow(net(triples(i, 0), triples(i, 2)), alpha);
              st2 += pow(net(triples(i, 1), triples(i, 0)), alpha) *
                  pow(net(triples(i, 1), triples(i, 2)), alpha);
              st3 += pow(net(triples(i, 2), triples(i, 0)), alpha) *
                  pow(net(triples(i, 2),triples(i, 1)), alpha);
          }
          to_return = st1 + st2 + st3;
      }
      return to_return;
  };

  // Function to calculate the number of in 2-stars
  double In2Star(arma::mat net,
                 arma::mat triples,
                 double alpha,
                 int together) {

      int number_of_triples = triples.n_rows;
      double st1 = 0;
      double st2 = 0;
      double st3 = 0;

      double to_return = 0;
      if (together == 1) {
          for (int i = 0; i < number_of_triples; ++i) {
              st1 += net(triples(i, 2), triples(i, 0)) * net(triples(i, 1),
                         triples(i, 0));
              st2 += net(triples(i, 2), triples(i, 1)) * net(triples(i, 0),
                         triples(i, 1));
              st3 += net(triples(i, 0), triples(i, 2)) * net(triples(i, 1),
                         triples(i, 2));
          }
          to_return = pow((st1 + st2 + st3), alpha);
      } else {
          for (int i = 0; i < number_of_triples; ++i) {
              st1 += pow(net(triples(i, 2), triples(i, 0)), alpha) *
                  pow(net(triples(i, 1), triples(i, 0)), alpha);
              st2 += pow(net(triples(i, 2), triples(i, 1)), alpha) *
                  pow(net(triples(i, 0), triples(i, 1)), alpha);
              st3 += pow(net(triples(i, 0), triples(i, 2)), alpha) *
                  pow(net(triples(i, 1),triples(i, 2)), alpha);
          }
          to_return = st1 + st2 + st3;
      }
      return to_return;
  };

  // Function to calculate the number of transitive triads
  double TTriads(arma::mat net,
                 arma::mat triples,
                 double alpha,
                 int together) {

      int number_of_triples = triples.n_rows;
      double st1 = 0;
      double st2 = 0;
      double st3 = 0;
      double st4 = 0;
      double st5 = 0;
      double st6 = 0;

      double to_return = 0;
      if (together == 1) {
          for (int i = 0; i < number_of_triples; ++i) {
              st1 += net(triples(i, 0), triples(i, 1)) * net(triples(i, 1),
                         triples(i, 2)) * net(triples(i, 0), triples(i, 2));
              st2 += net(triples(i, 0), triples(i, 1)) * net(triples(i, 2),
                         triples(i, 1)) * net(triples(i, 2), triples(i, 0));
              st3 += net(triples(i, 0), triples(i, 1)) * net(triples(i, 2),
                         triples(i, 1)) * net(triples(i, 0), triples(i, 2));
              st4 += net(triples(i, 1), triples(i, 0)) * net(triples(i, 1),
                         triples(i, 2)) * net(triples(i, 2), triples(i, 0));
              st5 += net(triples(i, 1), triples(i, 0)) * net(triples(i, 1),
                         triples(i, 2)) * net(triples(i, 0), triples(i, 2));
              st6 += net(triples(i, 1), triples(i, 0)) * net(triples(i, 2),
                         triples(i, 1)) * net(triples(i, 2), triples(i, 0));
          }
          to_return = pow((st1 + st2 + st3 + st4 + st5 + st6), alpha);
      } else {
          for (int i = 0; i < number_of_triples; ++i) {
              st1 += pow(net(triples(i, 0), triples(i, 1)), alpha) *
                  pow(net(triples(i, 1), triples(i, 2)), alpha) *
                  pow(net(triples(i, 0), triples(i, 2)), alpha);
              st2 += pow(net(triples(i, 0), triples(i, 1)), alpha) *
                  pow(net(triples(i, 2), triples(i, 1)), alpha) *
                  pow(net(triples(i, 2), triples(i, 0)), alpha);
              st3 += pow(net(triples(i, 0), triples(i, 1)), alpha) *
                  pow(net(triples(i, 2),triples(i, 1)), alpha) *
                  pow(net(triples(i, 0), triples(i, 2)), alpha);
              st4 += pow(net(triples(i, 1), triples(i, 0)), alpha) *
                  pow(net(triples(i, 1), triples(i, 2)), alpha) *
                  pow(net(triples(i, 2), triples(i, 0)), alpha);
              st5 += pow(net(triples(i, 1), triples(i, 0)), alpha) *
                  pow(net(triples(i, 1), triples(i, 2)), alpha) *
                  pow(net(triples(i, 0), triples(i, 2)), alpha);
              st6 += pow(net(triples(i, 1), triples(i, 0)), alpha) *
                  pow(net(triples(i, 2), triples(i, 1)), alpha) *
                  pow(net(triples(i, 2), triples(i, 0)), alpha);
          }
          to_return = st1 + st2 + st3 + st4 + st5 + st6;
      }
      return to_return;
  };

  // Function to calculate the number of closed triads
  double CTriads(arma::mat net,
                 arma::mat triples,
                 double alpha,
                 int together){

      int number_of_triples = triples.n_rows;
      double st1 = 0;
      double st2 = 0;

      double to_return = 0;
      if (together == 1) {
          for (int i = 0; i < number_of_triples; ++i) {
              st1 += net(triples(i, 0), triples(i, 1)) * net(triples(i, 1),
                         triples(i, 2)) * net(triples(i, 2), triples(i, 0));
              st2 += net(triples(i, 1), triples(i, 0)) * net(triples(i, 2),
                         triples(i,1)) * net(triples(i, 0), triples(i, 2));
          }
          to_return = pow((st1 + st2), alpha);
      } else {
          for (int i = 0; i < number_of_triples; ++i) {
              st1 += pow(net(triples(i, 0), triples(i, 1)), alpha) *
                  pow(net(triples(i, 1), triples(i, 2)), alpha) *
                  pow(net(triples(i, 2), triples(i, 0)), alpha);
              st2 += pow(net(triples(i, 1), triples(i, 0)), alpha) *
                  pow(net(triples(i, 2), triples(i,1)), alpha) *
                  pow(net(triples(i, 0), triples(i, 2)), alpha);
          }
          to_return = st1 + st2;
      }
      return to_return;
  };

  // Function to calculate the number of reciprocated edges
  double Recip(arma::mat net,
               arma::mat pairs,
               double alpha,
               int together) {

      int number_of_pairs = pairs.n_rows;
      double st1 = 0;

      double to_return = 0;
      if(together == 1) {
          for (int i = 0; i < number_of_pairs; ++i) {
              st1 += net(pairs(i, 0), pairs(i, 1)) * net(pairs(i, 1), pairs(i, 0));
          }
          to_return = pow(st1, alpha);
      } else {
          for (int i = 0; i < number_of_pairs; ++i) {
              st1 += pow(net(pairs(i, 0), pairs(i, 1)), alpha) *
                  pow(net(pairs(i, 1), pairs(i, 0)), alpha);
          }
          to_return = st1;
      }


      return to_return;
  };

  // Function to calculate the density of the network
  double EdgeDensity(arma::mat net,
                     arma::mat pairs,
                     double alpha,
                     int together) {

      int number_of_pairs = pairs.n_rows;
      double st1 = 0;
      double st2 = 0;
      double to_return = 0;

      for (int i = 0; i < number_of_pairs; ++i) {
          st1 += net(pairs(i, 0), pairs(i, 1));
          st1 += net(pairs(i, 1), pairs(i, 0));
          st2 += pow(net(pairs(i, 0), pairs(i, 1)), alpha);
          st2 += pow(net(pairs(i, 1), pairs(i, 0)), alpha);
      }

      if (together == 1) {
          to_return = pow(st1, alpha);
      } else {
          to_return = st2;
      }
      return to_return;
  };

  // Function that will calculate h statistics
  double CalculateNetworkStatistics(arma::mat current_network,
                                    arma::vec statistics_to_use,
                                    arma::vec thetas,
                                    arma::mat triples,
                                    arma::mat pairs,
                                    arma::vec alphas,
                                    int together) {

    double to_return = 0;

    if (statistics_to_use[0] == 1) {
      to_return += thetas[0] * Out2Star(current_network, triples, alphas[0],
                                        together);
    }
    if (statistics_to_use[1] == 1) {
      to_return += thetas[1] * In2Star(current_network, triples, alphas[1],
                                       together);
    }
    if (statistics_to_use[2] == 1) {
      to_return += thetas[2] * CTriads(current_network, triples, alphas[2],
                                       together);
    }
    if (statistics_to_use[3] == 1) {
      to_return += thetas[3] * Recip(current_network, pairs, alphas[3],
                                     together);
    }
    if (statistics_to_use[4] == 1) {
      to_return += thetas[4] * TTriads(current_network, triples, alphas[4],
                                       together);
    }
    if (statistics_to_use[5] == 1) {
      to_return += thetas[5] * EdgeDensity(current_network, pairs,
                                           alphas[5], together);
    }
    return to_return;
  };


  // Function that will calculate h statistics
  double integrand(arma::mat current_network,
                  arma::vec statistics_to_use,
                  arma::vec thetas,
                  arma::mat triples,
                  arma::mat pairs,
                  arma::vec alphas,
                  int together,
                  int sender,
                  int recipient,
                  double edge_value) {

    if (edge_value != -1) {
      // assign the current edge value in the current spot
      current_network(sender,recipient) = edge_value;
    }

    // calculate theta times the current statistics.
    double to_return = wobj::CalculateNetworkStatistics(
      current_network, statistics_to_use, thetas, triples, pairs,
      alphas, together);

    return to_return;
  };


  // ***********************************************************************//
  //     Constructor for Parallel Token Topic Distribution Generator        //
  // ***********************************************************************//

  // create a RcppParallel::Worker struct that we can use to fill in the
  // entries in our token topic distribution in parallel
  struct Parallel_Integrand : public RcppParallel::Worker {

    // Instantiate all of our input variables which will then be initialized
    // via the Function initializer below. This is how the struct makes sure
    // that everything has the right type
    arma::mat current_network;
    arma::vec statistics_to_use;
    arma::vec thetas;
    arma::mat triples;
    arma::mat pairs;
    arma::vec alphas;
    int together;
    int sender;
    int recipient;
    arma::vec integration_interval;

    // We need to initialize the output vector to an RcppParallel::RVector
    // vector (which is compatible with an Rcpp::NumericVector but not an
    // arma::vec). This is what will be implicitly returned
    RcppParallel::RVector<double> return_dist;

    // Function initializer
    Parallel_Integrand(arma::mat current_network,
                       arma::vec statistics_to_use,
                       arma::vec thetas,
                       arma::mat triples,
                       arma::mat pairs,
                       arma::vec alphas,
                       int together,
                       int sender,
                       int recipient,
                       arma::vec integration_interval,
                       Rcpp::NumericVector return_dist)
      : current_network(current_network),
        statistics_to_use(statistics_to_use),
        thetas(thetas),
        triples(triples),
        pairs(pairs),
        alphas(alphas),
        together(together),
        sender(sender),
        recipient(recipient),
        integration_interval(integration_interval),
        return_dist(return_dist) {}

    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
        double integ = wobj::integrand(current_network,
                                     statistics_to_use,
                                     thetas,
                                     triples,
                                     pairs,
                                     alphas,
                                     together,
                                     sender,
                                     recipient,
                                     integration_interval[i]);
        return_dist[i] = integ;
      }
    }
  };

  // ***********************************************************************//
  //             Calculate token topic probabilities in parallel            //
  // ***********************************************************************//
  arma::vec parallel_integration(
      arma::mat current_network,
      arma::vec statistics_to_use,
      arma::vec thetas,
      arma::mat triples,
      arma::mat pairs,
      arma::vec alphas,
      int together,
      int sender,
      int recipient,
      arma::vec integration_interval) {

    int innterval_size = integration_interval.n_elem;
    // the vector that will be operated on by RcppParallel::parallelFor.
    // This vector must be an Rcpp::NumericVector or the parallelization
    // will not work. THis makes everything slower, but I have not found an
    // alternative
    Rcpp::NumericVector output_vec(innterval_size);
    // Once we have the token topic distribution back from the parallel for
    // function, we will put it in to an arma::vec of the same length.
    arma::vec return_vec = arma::zeros(innterval_size);

    // create the RcppParallel::Worker function which will iterate over the
    // distribution
    Parallel_Integrand Parallel_Integrand(
        current_network,
        statistics_to_use,
        thetas,
        triples,
        pairs,
        alphas,
        together,
        sender,
        recipient,
        integration_interval,
        output_vec);

    // Call our Parallel_Token_Topic_Distribution with parallelFor
    RcppParallel::parallelFor(0,
                              innterval_size,
                              Parallel_Integrand);

    // Take the NumericVector output and put it in the arma::vec to actually
    // return
    for (int i = 0; i < innterval_size; i++) {
      // write to output vector
      return_vec[i] = output_vec[i];
    }

    return return_vec;
  }

  double log_sum_exp_integrator (arma::mat current_network,
                                 arma::vec statistics_to_use,
                                 arma::vec thetas,
                                 arma::mat triples,
                                 arma::mat pairs,
                                 arma::vec alphas,
                                 int together,
                                 int sender,
                                 int recipient,
                                 arma::vec integration_interval,
                                 bool parallel) {

    int num_evaluations = integration_interval.n_elem;
    arma::vec integral_evaluations = arma::zeros(num_evaluations);

    if (parallel) {
      integral_evaluations = parallel_integration(
        current_network,
        statistics_to_use,
        thetas,
        triples,
        pairs,
        alphas,
        together,
        sender,
        recipient,
        integration_interval);
    } else {
      for (int i = 0; i < num_evaluations; ++i) {
        integral_evaluations[i] = integrand(current_network,
                                            statistics_to_use,
                                            thetas,
                                            triples,
                                            pairs,
                                            alphas,
                                            together,
                                            sender,
                                            recipient,
                                            integration_interval[i]);
      }
    }

    // find the max on the interval
    double max_val = arma::max(integral_evaluations);

    // now calculate log sum exp
    double sum_term = 0;
    for (int i = 0; i < num_evaluations; ++i) {
      sum_term += exp(integral_evaluations[i] - max_val);
    }

    // we are taking the average
    sum_term = sum_term/double(num_evaluations);

    sum_term = max_val + log(sum_term);

    return sum_term;

  }





} //end of wobj namespace

using std::log;

// [[Rcpp::export]]
double weighted_mple_objective ( int number_of_nodes,
                          arma::vec statistics_to_use,
                          arma::mat current_network,
                          arma::vec thetas,
                          arma::mat triples,
                          arma::mat pairs,
                          arma::vec alphas,
                          int together,
                          arma::vec integration_interval,
                          bool parallel) {

  double objective = 0;
  // Calculate the
  for (int i = 0; i < number_of_nodes; ++i) {
    for (int j = 0; j < number_of_nodes; ++j) {
      double temp1 = wobj::integrand(current_network,
                               statistics_to_use,
                               thetas,
                               triples,
                               pairs,
                               alphas,
                               together,
                               i,
                               j,
                               -1);

      double temp2 = wobj::log_sum_exp_integrator(current_network,
                               statistics_to_use,
                               thetas,
                               triples,
                               pairs,
                               alphas,
                               together,
                               i,
                               j,
                               integration_interval,
                               parallel);

      objective += temp1 - temp2;
    }
  }

  return objective;
}
