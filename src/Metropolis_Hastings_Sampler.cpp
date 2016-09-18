// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]

#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <boost/config/no_tr1/cmath.hpp>
#include <istream>
#include <iosfwd>
#include <boost/assert.hpp>
#include <boost/limits.hpp>
#include <boost/static_assert.hpp>
#include <boost/random/detail/config.hpp>
#include <boost/random/detail/operators.hpp>
#include <boost/random/uniform_01.hpp>


using namespace Rcpp;

namespace mjd {

  namespace random {
  // Distributed under the Boost Software License, Version 1.0.
  //    (See http://www.boost.org/LICENSE_1_0.txt)
  // Coppied from the Boost libraries because they use an assert statement
  // which will not pass r cmd check.
  // http://www.boost.org/doc/libs/1_53_0/boost/random/normal_distribution.hpp
  // deterministic Box-Muller method, uses trigonometric functions

  /**
  * Instantiations of class template normal_distribution model a
  * \random_distribution. Such a distribution produces random numbers
  * @c x distributed with probability density function
  * \f$\displaystyle p(x) =
  *   \frac{1}{\sqrt{2\pi\sigma}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}
  * \f$,
  * where mean and sigma are the parameters of the distribution.
  */
  template<class RealType = double>
  class normal_distribution
  {
  public:
    typedef RealType input_type;
    typedef RealType result_type;

    class param_type {
    public:
      typedef normal_distribution distribution_type;

      /**
      * Constructs a @c param_type with a given mean and
      * standard deviation.
      *
      * Requires: sigma >= 0
      */
      explicit param_type(RealType mean_arg = RealType(0.0),
                          RealType sigma_arg = RealType(1.0))
        : _mean(mean_arg),
          _sigma(sigma_arg)
      {}

      /** Returns the mean of the distribution. */
      RealType mean() const { return _mean; }

      /** Returns the standand deviation of the distribution. */
      RealType sigma() const { return _sigma; }

      /** Writes a @c param_type to a @c std::ostream. */
      BOOST_RANDOM_DETAIL_OSTREAM_OPERATOR(os, param_type, parm)
      { os << parm._mean << " " << parm._sigma ; return os; }

      /** Reads a @c param_type from a @c std::istream. */
      BOOST_RANDOM_DETAIL_ISTREAM_OPERATOR(is, param_type, parm)
      { is >> parm._mean >> std::ws >> parm._sigma; return is; }

      /** Returns true if the two sets of parameters are the same. */
      BOOST_RANDOM_DETAIL_EQUALITY_OPERATOR(param_type, lhs, rhs)
      { return lhs._mean == rhs._mean && lhs._sigma == rhs._sigma; }

      /** Returns true if the two sets of parameters are the different. */
      BOOST_RANDOM_DETAIL_INEQUALITY_OPERATOR(param_type)

    private:
      RealType _mean;
      RealType _sigma;
    };

    /**
    * Constructs a @c normal_distribution object. @c mean and @c sigma are
    * the parameters for the distribution.
    *
    * Requires: sigma >= 0
    */
    explicit normal_distribution(const RealType& mean_arg = RealType(0.0),
                                 const RealType& sigma_arg = RealType(1.0))
      : _mean(mean_arg), _sigma(sigma_arg),
        _r1(0), _r2(0), _cached_rho(0), _valid(false)
    {}

    /**
    * Constructs a @c normal_distribution object from its parameters.
    */
    explicit normal_distribution(const param_type& parm)
      : _mean(parm.mean()), _sigma(parm.sigma()),
        _r1(0), _r2(0), _cached_rho(0), _valid(false)
    {}

    /**  Returns the mean of the distribution. */
    RealType mean() const { return _mean; }
    /** Returns the standard deviation of the distribution. */
    RealType sigma() const { return _sigma; }

    /** Returns the smallest value that the distribution can produce. */
    RealType min BOOST_PREVENT_MACRO_SUBSTITUTION () const
    { return -std::numeric_limits<RealType>::infinity(); }
    /** Returns the largest value that the distribution can produce. */
    RealType max BOOST_PREVENT_MACRO_SUBSTITUTION () const
    { return std::numeric_limits<RealType>::infinity(); }

    /** Returns the parameters of the distribution. */
    param_type param() const { return param_type(_mean, _sigma); }
    /** Sets the parameters of the distribution. */
    void param(const param_type& parm)
    {
      _mean = parm.mean();
      _sigma = parm.sigma();
      _valid = false;
    }

    /**
    * Effects: Subsequent uses of the distribution do not depend
    * on values produced by any engine prior to invoking reset.
    */
    void reset() { _valid = false; }

    /**  Returns a normal variate. */
    template<class Engine>
    result_type operator()(Engine& eng)
    {
      using std::sqrt;
      using std::log;
      using std::sin;
      using std::cos;

      if(!_valid) {
        _r1 = boost::uniform_01<RealType>()(eng);
        _r2 = boost::uniform_01<RealType>()(eng);
        _cached_rho = sqrt(-result_type(2) * log(result_type(1)-_r2));
        _valid = true;
      } else {
        _valid = false;
      }
      // Can we have a boost::mathconst please?
      const result_type pi = result_type(3.14159265358979323846);

      return _cached_rho * (_valid ?
                            cos(result_type(2)*pi*_r1) :
                              sin(result_type(2)*pi*_r1))
        * _sigma + _mean;
    }

    /** Returns a normal variate with parameters specified by @c param. */
    template<class URNG>
    result_type operator()(URNG& urng, const param_type& parm)
    {
      return normal_distribution(parm)(urng);
    }

    /** Writes a @c normal_distribution to a @c std::ostream. */
    BOOST_RANDOM_DETAIL_OSTREAM_OPERATOR(os, normal_distribution, nd)
    {
      os << nd._mean << " " << nd._sigma << " "
         << nd._valid << " " << nd._cached_rho << " " << nd._r1;
      return os;
    }

    /** Reads a @c normal_distribution from a @c std::istream. */
    BOOST_RANDOM_DETAIL_ISTREAM_OPERATOR(is, normal_distribution, nd)
    {
      is >> std::ws >> nd._mean >> std::ws >> nd._sigma
         >> std::ws >> nd._valid >> std::ws >> nd._cached_rho
         >> std::ws >> nd._r1;
         return is;
    }

    /**
    * Returns true if the two instances of @c normal_distribution will
    * return identical sequences of values given equal generators.
    */
    BOOST_RANDOM_DETAIL_EQUALITY_OPERATOR(normal_distribution, lhs, rhs)
    {
      return lhs._mean == rhs._mean && lhs._sigma == rhs._sigma
      && lhs._valid == rhs._valid
      && (!lhs._valid || (lhs._r1 == rhs._r1 && lhs._r2 == rhs._r2));
    }

    /**
    * Returns true if the two instances of @c normal_distribution will
    * return different sequences of values given equal generators.
    */
    BOOST_RANDOM_DETAIL_INEQUALITY_OPERATOR(normal_distribution)

  private:
    RealType _mean, _sigma;
    RealType _r1, _r2, _cached_rho;
    bool _valid;

  };

} // namespace random

using random::normal_distribution;
using std::pow;
using std::exp;
using std::sqrt;

// add in the functions I wrote for correlation networks
arma::mat partials_to_correlations(arma::mat partial_correlations){
  int nrow = partial_correlations.n_rows;
  arma::mat correlations = arma::ones(nrow,nrow);
  correlations.diag(1) = partial_correlations.diag(1);
  correlations.diag(-1) = partial_correlations.diag(-1);

  for (int k = 2; k < (nrow); ++k) {
    for (int i = 1; i < (nrow - k + 1); ++i) {
      arma::mat R2 = correlations.submat(i,i,(i + k -2),(i + k -2));
      arma::mat R2_ones = arma::ones(R2.n_cols,R2.n_cols);
      arma::rowvec r1 = correlations(i-1,arma::span(i,(i + k -2)));
      arma::rowvec r3 =  correlations((i + k -1),arma::span(i,(i + k- 2)));
      arma::mat R2_solved = arma::inv(R2);

      arma::mat temp1 = 1 - r1 * R2_solved * r1.t();
      arma::mat temp2 =  1 - r3 * R2_solved * r3.t();
      arma::mat D =  sqrt(temp1 * temp2);
      arma::mat temp3 = r1 * R2_solved * r3.t() + partial_correlations(i -1, i + k-1)*D;
      correlations((i-1), (i + k-1)) = temp3(0,0);
      correlations((i + k-1), (i-1)) = correlations(i-1, i + k-1);
    }
  }
  return correlations;
}


arma::mat bounded_to_correlations(arma::mat bounded_network){
  //transform back to the partial space
  arma::mat partials = 2 * bounded_network -1;
  int temp = partials.n_rows;
  partials.diag() = arma::ones(temp);
  //transform to correlation space
  arma::mat correlations =  mjd::partials_to_correlations(partials);
  return correlations;
}


double jacobian(arma::mat partial_correlations){
  int nrow = partial_correlations.n_rows;
  arma::vec corrs_1 = partial_correlations.diag(1);
  arma::vec temp = pow(corrs_1,2);
  arma::vec temp2 = pow((1 - temp),(nrow -2));
  double prod_1 =  arma::prod(temp2);
  double prod_2 = 1;

  for (int k = 2; k < (nrow -1); ++k) {
    for (int i = 0; i < (nrow - k); ++i) {
      double temp3 = pow(partial_correlations(i,(i + k)),2);
      prod_2 = prod_2 * pow((1 - temp3),(nrow - 1 - k));
    }
  }

  double temp4 = pow(prod_1,(nrow - 2));
  double result = pow(temp4*prod_2,0.5);
  return result;
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


struct Parallel_CalculateNetworkStatistics : public RcppParallel::Worker {
  arma::mat current_network;
  arma::vec statistics_to_use;
  arma::vec thetas;
  arma::mat triples;
  arma::mat pairs;
  arma::vec alphas;
  int together;
  RcppParallel::RVector<double> return_dist;

  Parallel_CalculateNetworkStatistics(arma::mat current_network,
                                      arma::vec statistics_to_use,
                                      arma::vec thetas,
                                      arma::mat triples,
                                      arma::mat pairs,
                                      arma::vec alphas,
                                      int together,
                                      Rcpp::NumericVector return_dist)
    : current_network(current_network),
      statistics_to_use(statistics_to_use),
      thetas(thetas),
      triples(triples),
      pairs(pairs),
      alphas(alphas),
      together(together),
      return_dist(return_dist) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      if (i == 0 & statistics_to_use[0] == 1) {
        return_dist[i] = thetas[0] * Out2Star(current_network, triples, alphas[0],
                                              together);
      }
      if (i == 1 & statistics_to_use[1] == 1) {
        return_dist[i] = thetas[1] * In2Star(current_network, triples, alphas[1],
                                             together);
      }
      if (i == 2 & statistics_to_use[2] == 1) {
        return_dist[i] = thetas[2] * CTriads(current_network, triples, alphas[2],
                                             together);
      }
      if (i == 3 & statistics_to_use[3] == 1) {
        return_dist[i] = thetas[3] * Recip(current_network, pairs, alphas[3],
                                           together);
      }
      if (i == 4 & statistics_to_use[4] == 1) {
        return_dist[i] = thetas[4] * TTriads(current_network, triples, alphas[4],
                                             together);
      }
      if (i == 5 & statistics_to_use[5] == 1) {
        return_dist[i] = thetas[5] * EdgeDensity(current_network, pairs,
                                                 alphas[5], together);
      }
    }
  }
};

double parallel_CalculateNetworkStatistics(
    arma::mat current_network,
    arma::vec statistics_to_use,
    arma::vec thetas,
    arma::mat triples,
    arma::mat pairs,
    arma::vec alphas,
    int together) {

  int number_of_stats = statistics_to_use.n_elem;
  Rcpp::NumericVector output_vec(number_of_stats);

  Parallel_CalculateNetworkStatistics Parallel_CalculateNetworkStatistics(
      current_network,
      statistics_to_use,
      thetas,
      triples,
      pairs,
      alphas,
      together,
      output_vec);

  RcppParallel::parallelFor(0,
                            number_of_stats,
                            Parallel_CalculateNetworkStatistics);

  double return_val = 0;
  for (int i = 0; i < number_of_stats; i++) {
    return_val += output_vec[i];
  }

  return return_val;
}


// Function that will calculate h statistics
double CalculateNetworkStatistics(arma::mat current_network,
                                  arma::vec statistics_to_use,
                                  arma::vec thetas,
                                  arma::mat triples,
                                  arma::mat pairs,
                                  arma::vec alphas,
                                  int together,
                                  bool parallel) {

  double to_return = 0;
  if (parallel) {
    to_return = mjd::parallel_CalculateNetworkStatistics(
      current_network,
      statistics_to_use,
      thetas,
      triples,
      pairs,
      alphas,
      together);
  } else {
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


// Function that will calculate h statistics
arma::vec h_function_and_statistics(arma::mat current_network,
                                    arma::vec statistics_to_use,
                                    arma::vec thetas,
                                    arma::mat triples,
                                    arma::mat pairs,
                                    arma::vec alphas,
                                    int together) {

  arma::vec to_return = arma::zeros(7);
  to_return[0] = Out2Star(current_network, triples, alphas[0], together);
  to_return[1] = In2Star(current_network, triples, alphas[1], together);
  to_return[2] = CTriads(current_network, triples, alphas[2], together);
  to_return[3] = Recip(current_network, pairs, alphas[3], together);
  to_return[4] = TTriads(current_network, triples, alphas[4], together);
  to_return[5] = EdgeDensity(current_network, pairs, alphas[5], together);

  double h_value = 0;
  if (statistics_to_use[0] == 1) {
    h_value += thetas[0] * to_return[0];
  }
  if (statistics_to_use[1] == 1) {
    h_value += thetas[1] * to_return[1];
  }
  if (statistics_to_use[2] == 1) {
    h_value += thetas[2] * to_return[2];
  }
  if (statistics_to_use[3] == 1) {
    h_value += thetas[3] * to_return[3];
  }
  if (statistics_to_use[4] == 1) {
    h_value += thetas[4] * to_return[4];
  }
  if (statistics_to_use[5] == 1) {
    h_value += thetas[5] * to_return[5];
  }
  // put in the h value
  to_return[6] = h_value;
  return to_return;
};

} //end of mjd namespace

using std::log;

// [[Rcpp::export]]
List Metropolis_Hastings_Sampler (int number_of_iterations,
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
          int seed,
          int number_of_samples_to_store,
          int using_correlation_network,
          int undirect_network,
          bool parallel) {

  // Allocate variables and data structures
  double variance = shape_parameter;
  int list_length = 9;
  List to_return(list_length);
  int number_of_thetas = statistics_to_use.n_elem;
  int MH_Counter = 0;
  int Storage_Counter = 0;
  bool network_did_not_change = false;
  double previous_h_function_value = 0;
  arma::vec Accept_or_Reject = arma::zeros (number_of_iterations);
  arma::vec Log_Prob_Accept = arma::zeros (number_of_iterations);
  arma::vec P_Ratios = arma::zeros (number_of_iterations);
  arma::vec Q_Ratios = arma::zeros (number_of_iterations);
  arma::vec Proposed_Density = arma::zeros (number_of_iterations);
  arma::vec Current_Density = arma::zeros (number_of_iterations);
  arma::cube Network_Samples = arma::zeros (number_of_nodes, number_of_nodes,
      number_of_samples_to_store);
  arma::vec Mean_Edge_Weights = arma::zeros (number_of_samples_to_store);
  arma::mat Save_H_Statistics = arma::zeros (number_of_samples_to_store,
      number_of_thetas);
  arma::mat current_edge_weights = initial_network;
  arma::mat corr_current_edge_weights = arma::zeros (number_of_nodes, number_of_nodes);

  // deal with the case where we have a correlation network.
  if(using_correlation_network == 1){
    current_edge_weights.diag() = arma::ones(number_of_nodes);
    undirect_network = 1;
  }


  // Set RNG and define uniform distribution
  boost::mt19937 generator(seed);
  boost::uniform_01<double> uniform_distribution;
  // Outer loop over the number of samples
  for (int n = 0; n < number_of_iterations; ++n) {
    double log_prob_accept = 0;
    arma::mat proposed_edge_weights = current_edge_weights;

    // deal with the case where we have an undirected network.
    if(undirect_network == 1){
      // Run loop to sample new edge weights
      for (int i = 0; i < number_of_nodes; ++i) {
        for (int j = 0; j < i; ++j) {
          if (i != j) {

            double log_probability_of_current_under_new = 0;
            double log_probability_of_new_under_current = 0;
            //draw a new edge value centered at the old edge value
            double current_edge_value = current_edge_weights(i,j);
            //draw from a truncated normal
            mjd::normal_distribution<double> proposal(current_edge_value,variance);
            int in_zero_one = 0;
            //NumericVector new_edge_value = 0.5;
            double new_edge_value = 0.5;
            while(in_zero_one == 0){
              new_edge_value = proposal(generator);
              if(new_edge_value > 0 & new_edge_value < 1){
                in_zero_one = 1;
              }
            }
            // if (new_edge_value > 0.999) {
            //   new_edge_value = 0.999;
            // }
            // if (new_edge_value < 0.001) {
            //   new_edge_value= 0.001;
            // }
            // calculate the probability of the new edge under current beta dist
            double lower_bound = R::pnorm(0,current_edge_value,variance, 1, 0);
            double upper_bound = R::pnorm(1,current_edge_value,variance, 1, 0);
            double raw_prob = R::dnorm(new_edge_value,current_edge_value,variance,0);
            double prob_new_edge_under_old = (raw_prob/(upper_bound - lower_bound));
            // calculate the probability of the current edge under new beta dist
            lower_bound = R::pnorm(0,new_edge_value,variance, 1, 0);
            upper_bound = R::pnorm(1,new_edge_value,variance, 1, 0);
            raw_prob = R::dnorm(current_edge_value,new_edge_value,variance,0);
            double prob_old_edge_under_new = (raw_prob/(upper_bound - lower_bound));
            //save everything
            proposed_edge_weights(i,j) = new_edge_value;
            proposed_edge_weights(j,i) = new_edge_value;
            log_probability_of_new_under_current = log(prob_new_edge_under_old);
            log_probability_of_current_under_new = log(prob_old_edge_under_new);

            // Calculate acceptance probability
            log_prob_accept += (log_probability_of_current_under_new
                                  - log_probability_of_new_under_current);

          }
        }
      }

    }else{
      // int counter =  0;
      // Run loop to sample new edge weights
      for (int i = 0; i < number_of_nodes; ++i) {
        for (int j = 0; j < number_of_nodes; ++j) {
          if (i != j) {
            double log_probability_of_current_under_new = 0;
            double log_probability_of_new_under_current = 0;
            //draw a new edge value centered at the old edge value
            double current_edge_value = current_edge_weights(i,j);
            //draw from a truncated normal
            mjd::normal_distribution<double> proposal(current_edge_value,variance);
            int in_zero_one = 0;
            //NumericVector new_edge_value = 0.5;
            double new_edge_value = 0.5;
            while(in_zero_one == 0){
              new_edge_value = proposal(generator);
              if(new_edge_value > 0 & new_edge_value < 1){
                in_zero_one = 1;
              }
            }

            // calculate the probability of the new edge under current beta dist
            double lower_bound = R::pnorm(0,current_edge_value,variance, 1, 0);
            double upper_bound = R::pnorm(1,current_edge_value,variance, 1, 0);
            double raw_prob = R::dnorm(new_edge_value,current_edge_value,variance,0);
            double prob_new_edge_under_old = (raw_prob/(upper_bound - lower_bound));

            // calculate the probability of the current edge under new beta dist
            lower_bound = R::pnorm(0,new_edge_value,variance, 1, 0);
            upper_bound = R::pnorm(1,new_edge_value,variance, 1, 0);
            raw_prob = R::dnorm(current_edge_value,new_edge_value,variance,0);
            double prob_old_edge_under_new = (raw_prob/(upper_bound - lower_bound));

            //if(counter < 100){
            //  Rcpp::Rcout << "lower_bound " << lower_bound <<
            //    " upper_bound  " << upper_bound <<
            //      " raw_prob  " << raw_prob << std::endl;
            //}

            //save everything
            proposed_edge_weights(i,j) = new_edge_value;
            log_probability_of_new_under_current = log(prob_new_edge_under_old);
            log_probability_of_current_under_new = log(prob_old_edge_under_new);

            // Calculate acceptance probability
            log_prob_accept += (log_probability_of_current_under_new
                                  - log_probability_of_new_under_current);

            //if(counter < 100){
            //  Rcpp::Rcout << "prob_old_edge_under_new " << prob_old_edge_under_new <<
            //    "prob_new_edge_under_old  " << prob_new_edge_under_old << std::endl;
            //  Rcpp::Rcout << "Isfinite: " << std::isfinite(log_prob_accept) <<
            //    "Isnan: " << std::isnan(log_prob_accept)<<
            //      "value: " << log_prob_accept << std::endl;
            //}
            //counter += 1
          }
        }
      }
    } //end of condition for whether we are using a correlation network

    double proposed_addition = 0;
    double current_addition = 0;

    if(using_correlation_network == 1){
      arma::mat corr_proposed_edge_weights = mjd::bounded_to_correlations(proposed_edge_weights);
      arma::mat corr_current_edge_weights = mjd::bounded_to_correlations(current_edge_weights);
      proposed_addition = mjd::CalculateNetworkStatistics(
        corr_proposed_edge_weights, statistics_to_use, thetas, triples, pairs,
        alphas, together, parallel);
      current_addition = mjd::CalculateNetworkStatistics(
        corr_current_edge_weights, statistics_to_use, thetas, triples, pairs,
        alphas, together, parallel);
    }else{
      proposed_addition = mjd::CalculateNetworkStatistics(
        proposed_edge_weights, statistics_to_use, thetas, triples, pairs,
        alphas, together, parallel);
      // only calculate the h function if we updated the network last round
      // otherwise use the cached value.
      if (network_did_not_change) {
        current_addition = previous_h_function_value;
      } else {
        current_addition = mjd::CalculateNetworkStatistics(
          current_edge_weights, statistics_to_use, thetas, triples, pairs,
          alphas, together, parallel);
        previous_h_function_value = current_addition ;
      }
    }

    // store some additional diagnostics h value is the last entry
    P_Ratios[n] = (proposed_addition - current_addition);
    Q_Ratios[n] = log_prob_accept;

    double total_edges = double(number_of_nodes * (number_of_nodes - 1));
    double temp1 = arma::accu(proposed_edge_weights);
    Proposed_Density[n] = temp1/total_edges;
    double temp2 = arma::accu(current_edge_weights);
    Current_Density[n] = temp2/total_edges;

    log_prob_accept += (proposed_addition- current_addition);

    if(using_correlation_network == 1){
      // now add in the bit about Jacobians
      double numerator = log(mjd::jacobian(2*proposed_edge_weights-1));
      double denominator = log(mjd::jacobian(2*current_edge_weights-1));
      log_prob_accept += numerator - denominator;
    }

    double rand_num = uniform_distribution(generator);
    double lud = 0;
    lud = log(rand_num);

    double accept_proportion = 0;
    // Accept or reject the new proposed positions
    if (log_prob_accept < lud) {
      accept_proportion +=0;
      network_did_not_change = true;
    } else {
      accept_proportion +=1;
      network_did_not_change = false;
      for (int i = 0; i < number_of_nodes; ++i) {
          for (int j = 0; j < number_of_nodes; ++j) {
              if (i != j) {
                  double temp = proposed_edge_weights(i, j);
                  current_edge_weights(i, j) = temp;
              }
          }
      }
    }

    Log_Prob_Accept[n] = log_prob_accept;
    Accept_or_Reject[n] = accept_proportion;
    Storage_Counter += 1;

    // Save network statistics
    if (Storage_Counter == take_sample_every) {

      if(using_correlation_network == 1){
        arma::mat corr_current_edge_weights = mjd::bounded_to_correlations(current_edge_weights);
        arma::vec save_stats = mjd::save_network_statistics(corr_current_edge_weights,
                                                            triples, pairs, alphas, together);
        for (int m = 0; m < 6; ++m) {
          Save_H_Statistics(MH_Counter, m) = save_stats[m];
        }
      }else{
        arma::vec save_stats = mjd::save_network_statistics(current_edge_weights,
                                                            triples, pairs, alphas, together);

        for (int m = 0; m < 6; ++m) {
          Save_H_Statistics(MH_Counter, m) = save_stats[m];
        }
      }


      double mew = 0;
      if(using_correlation_network == 1){
        corr_current_edge_weights = mjd::bounded_to_correlations(current_edge_weights);
      }

      for (int i = 0; i < number_of_nodes; ++i) {
        for (int j = 0; j < number_of_nodes; ++j) {
          if (i != j) {
            if(using_correlation_network == 1){
              //we use this trick to break the referencing
              double temp = corr_current_edge_weights(i, j);
              Network_Samples(i, j, MH_Counter) = temp;
              mew += temp;
            }else{
              //we use this trick to break the referencing
              double temp = current_edge_weights(i, j);
              Network_Samples(i, j, MH_Counter) = temp;
              mew += temp;
            }
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
  to_return[4] = Log_Prob_Accept;
  to_return[5] = P_Ratios;
  to_return[6] = Q_Ratios;
  to_return[7] = Proposed_Density;
  to_return[8] = Current_Density;
  return to_return;
}
