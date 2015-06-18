
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <math.h>
#include <cmath>


using namespace Rcpp;

namespace mjd {
    
    // Returns the erf() of a value (not super precice, but ok)
    double erf(double x)
    {  
     double y = 1.0 / ( 1.0 + 0.3275911 * x);   
     return 1 - (((((
            + 1.061405429  * y
            - 1.453152027) * y
            + 1.421413741) * y
            - 0.284496736) * y 
            + 0.254829592) * y) 
            * exp (-x * x);      
    }
    
    // Returns the probability of x, given the distribution described by mu and sigma.
    double pdf(double x, double mu, double sigma)
    {
      //Constants
      static const double pi = 3.14159265; 
      return exp( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
    }
    
    // Returns the probability of [-inf,x] of a gaussian distribution
    double cdf(double x, double mu, double sigma)
    {
        return 0.5 * (1 + mjd::erf((x - mu) / (sigma * sqrt(2.))));
    }
    
  
  // Function to calculate the number of out 2-stars
  double Out2Star(arma::mat net, 
                  arma::mat triples, 
                  double alpha,
                  int together) {
                      
    int number_of_triples = triples.n_rows;
    double st1 = 0;
    double st2 = 0;
    double st3 = 0;
    
    for (int i = 0; i < number_of_triples; ++i) {
      st1 += net(triples(i, 0), triples(i, 1)) * net(triples(i, 0), 
          triples(i, 2));
      st2 += net(triples(i, 1), triples(i, 0)) * net(triples(i, 1), 
          triples(i, 2));
      st3 += net(triples(i, 2), triples(i, 0)) * net(triples(i, 2), 
          triples(i, 1));
    }
    
    double to_return = 0;
    if (together == 1) {
      to_return = pow((st1 + st2 + st3), alpha);
    } else {
      to_return = pow(st1, alpha) + pow(st2, alpha) + pow(st3, alpha);
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
    
    for (int i = 0; i < number_of_triples; ++i) {
      st1 += net(triples(i, 2), triples(i, 0)) * net(triples(i, 1), 
          triples(i, 0));
      st2 += net(triples(i, 2), triples(i, 1)) * net(triples(i, 0), 
          triples(i, 1));
      st3 += net(triples(i, 0), triples(i, 2)) * net(triples(i, 1), 
          triples(i, 2));
    }
    
    double to_return = 0;
    if (together == 1) {
      to_return = pow((st1 + st2 + st3), alpha);
    } else {
      to_return = pow(st1, alpha) + pow(st2, alpha) + pow(st3, alpha);
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
    
    double to_return = 0;
    if (together == 1) {
      to_return = pow((st1 + st2 + st3 + st4 + st5 + st6), alpha);
    } else {
      to_return = pow(st1, alpha) + pow(st2, alpha) + pow(st3, alpha) 
          + pow(st4, alpha) + pow(st5, alpha) + pow(st6, alpha);
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
    
    for (int i = 0; i < number_of_triples; ++i) {
      st1 += net(triples(i, 0), triples(i, 1)) * net(triples(i, 1), 
          triples(i, 2)) * net(triples(i, 2), triples(i, 0));
      st2 += net(triples(i, 1), triples(i, 0)) * net(triples(i, 2), 
          triples(i,1)) * net(triples(i, 0), triples(i, 2));
    }
    
    double to_return = 0;
    if (together == 1) {
      to_return = pow((st1 + st2), alpha);
    } else {
      to_return = pow(st1, alpha) + pow(st2, alpha);
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
    for (int i = 0; i < number_of_pairs; ++i) {
        st1 += net(pairs(i, 0), pairs(i, 1)) * net(pairs(i, 1), pairs(i, 0));
    }
    
    double to_return = pow(st1, alpha);
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
  
  // Function that will calculate and save all of the h statistics for a network
  arma::vec save_network_statistics(arma::mat current_network,
                                    arma::mat triples,
                                    arma::mat pairs,
                                    arma::vec alphas,
                                    int together) {
        
    arma::vec to_return = arma::zeros(6);
    to_return[0] = Out2Star(current_network, triples, alphas[0], together);
    to_return[1] = In2Star(current_network, triples, alphas[1], together);
    to_return[2] = CTriads(current_network, triples, alphas[2], together);
    to_return[3] = Recip(current_network, pairs, alphas[3], together);
    to_return[4] = TTriads(current_network, triples, alphas[4], together);
    to_return[5] = EdgeDensity(current_network, pairs, alphas[5], together);
    return to_return;
  };
  

  
} //end of mjd namespace


// [[Rcpp::export]]
List MH (int number_of_iterations, 
          double shape_parameter,
          int number_of_nodes,
          arma::vec statistics_to_use,
          arma::mat initial_network,
          int take_sample_every,
          arma::vec thetas,
          arma::mat triples,
          arma::mat pairs,
          arma::vec alphas,
          int together,
          int seed) {
        
  // Allocate variables and data structures
  double variance = shape_parameter;
  int list_length = 4;
  List to_return(list_length);
  int number_of_samples_to_store = ceil (number_of_iterations / 
      take_sample_every);
  int number_of_thetas = 6;
  int MH_Counter = 0; 
  int Storage_Counter = 0; 
  arma::vec Accept_or_Reject = arma::zeros (number_of_iterations);
  arma::cube Network_Samples = arma::zeros (number_of_nodes, number_of_nodes,
      number_of_samples_to_store);
  arma::vec Mean_Edge_Weights = arma::zeros (number_of_samples_to_store);
  arma::mat Save_H_Statistics = arma::zeros (number_of_samples_to_store, 
      number_of_thetas);
  arma::mat current_edge_weights = initial_network;

  // Set RNG and define uniform distribution 
  boost::mt19937_64 generator(seed);
  boost::random::uniform_real_distribution< >  uniform_distribution(0.0,1.0);
  Function report("Report");

  // Outer loop over the number of samples
  for (int n = 0; n < number_of_iterations; ++n) {
     //report(variance);
     //std::normal_distribution<double> test(0,variance);
     //double testval = test(generator);
     //report(testval);
    double log_prob_accept = 0;
    arma::mat proposed_edge_weights = current_edge_weights;
    // Run loop to sample new edge weights
    for (int i = 0; i < number_of_nodes; ++i) {
      for (int j = 0; j < number_of_nodes; ++j) {
        if (i != j) {
            
            double log_probability_of_current_under_new = 0;
            double log_probability_of_new_under_current = 0;
            //draw a new edge value centered at the old edge value
            double current_edge_value = current_edge_weights(i,j);
            //draw from a truncated normal
            boost::normal_distribution<double> proposal(current_edge_value,variance);
            int in_zero_one = 0;
            double new_edge_value = 0.5;
            while(in_zero_one == 0){
                new_edge_value = proposal(generator);
                if(new_edge_value > 0 & new_edge_value < 1){
                    in_zero_one = 1;
                }
            }
            if (new_edge_value > 0.999) {
                new_edge_value = 0.999;
            }    
            if (new_edge_value < 0.001) {
                new_edge_value = 0.001;
            }
            //report(new_edge_value);
            
            // calculate the probability of the new edge under current beta dist
            double lower_bound = mjd::cdf(0,current_edge_value,variance);
            double upper_bound = mjd::cdf(1,current_edge_value,variance);
            double raw_prob = mjd::pdf(new_edge_value,current_edge_value,variance);
            double prob_new_edge_under_old = (raw_prob/(upper_bound - lower_bound));
            
            // calculate the probability of the current edge under new beta dist
            lower_bound = mjd::cdf(0,new_edge_value,variance);
            upper_bound = mjd::cdf(1,new_edge_value,variance);
            raw_prob = mjd::pdf(current_edge_value,new_edge_value,variance);
            double prob_old_edge_under_new = (raw_prob/(upper_bound - lower_bound));
            
            //save everything
            proposed_edge_weights(i,j) = new_edge_value;  
            log_probability_of_new_under_current = log(prob_new_edge_under_old);
            log_probability_of_current_under_new = log(prob_old_edge_under_new);
            
            // Calculate acceptance probability
            log_prob_accept += (log_probability_of_current_under_new
              - log_probability_of_new_under_current);
          
        } 
      }
    }
    
    
    double proposed_addition = mjd::CalculateNetworkStatistics(
              proposed_edge_weights, statistics_to_use, thetas, triples, pairs,
              alphas, together);
    double current_addition = mjd::CalculateNetworkStatistics(
              current_edge_weights, statistics_to_use, thetas, triples, pairs,
              alphas, together);
    //report("H stats, Old and New");
    //report(n);
    //report(log_prob_accept);
    //report(current_addition);
    //report(proposed_addition);
    log_prob_accept += (proposed_addition - current_addition);
    //report("With H statistics");
    //report(log_prob_accept);
    double rand_num = uniform_distribution(generator);
    double lud = 0;
    lud = log(rand_num);
    
    double accept_proportion = 0;
    // Accept or reject the new proposed positions
    if (log_prob_accept < lud) {
      accept_proportion +=0;
    } else {
      accept_proportion +=1;
      for (int i = 0; i < number_of_nodes; ++i) {
          for (int j = 0; j < number_of_nodes; ++j) {
              if (i != j) {
                  double temp = proposed_edge_weights(i, j);
                  current_edge_weights(i, j) = temp;
              }
          }
      }
    }
    
    Accept_or_Reject[n] = accept_proportion;
    Storage_Counter += 1;
    
    // Save network statistics
    if (Storage_Counter == take_sample_every) {
        
      arma::vec save_stats = mjd::save_network_statistics(current_edge_weights,
          triples, pairs, alphas, together);
          
      for (int m = 0; m < 6; ++m) {
        Save_H_Statistics(MH_Counter, m) = save_stats[m];
      }

      double mew = 0;
     
      for (int i = 0; i < number_of_nodes; ++i) {
        for (int j = 0; j < number_of_nodes; ++j) {
          if (i != j) {
              
            //we use this trick to break the referencing
            double temp = current_edge_weights(i, j);
            Network_Samples(i, j, MH_Counter) = temp;
            mew += temp;
          }
        }
      }
        
      mew = mew / double(number_of_nodes * (number_of_nodes - 1));
      Mean_Edge_Weights[MH_Counter] = mew;
      Storage_Counter = 0;
      MH_Counter += 1;
    }
  }
  
  // Save the data and then return
  to_return[0] = Accept_or_Reject;
  to_return[1] = Network_Samples;
  to_return[2] = Save_H_Statistics;
  to_return[3] = Mean_Edge_Weights;
  return to_return;    
}
