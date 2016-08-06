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
using namespace arma;

namespace gergm {

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




// Function to calculate the number of transitive triads
arma::vec indiviual_triad_values(arma::mat net,
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

  arma::vec to_return = arma::zeros(number_of_triples);
  if (together == 1) {
    for (int i = 0; i < number_of_triples; ++i) {
      st1 = net(triples(i, 0), triples(i, 1)) * net(triples(i, 1),
                 triples(i, 2)) * net(triples(i, 0), triples(i, 2));
      st2 = net(triples(i, 0), triples(i, 1)) * net(triples(i, 2),
                 triples(i, 1)) * net(triples(i, 2), triples(i, 0));
      st3 = net(triples(i, 0), triples(i, 1)) * net(triples(i, 2),
                 triples(i, 1)) * net(triples(i, 0), triples(i, 2));
      st4 = net(triples(i, 1), triples(i, 0)) * net(triples(i, 1),
                 triples(i, 2)) * net(triples(i, 2), triples(i, 0));
      st5 = net(triples(i, 1), triples(i, 0)) * net(triples(i, 1),
                 triples(i, 2)) * net(triples(i, 0), triples(i, 2));
      st6 = net(triples(i, 1), triples(i, 0)) * net(triples(i, 2),
                 triples(i, 1)) * net(triples(i, 2), triples(i, 0));
      to_return[i] = pow((st1 + st2 + st3 + st4 + st5 + st6), alpha);
    }

  } else {
    for (int i = 0; i < number_of_triples; ++i) {
      st1 = pow(net(triples(i, 0), triples(i, 1)), alpha) *
        pow(net(triples(i, 1), triples(i, 2)), alpha) *
        pow(net(triples(i, 0), triples(i, 2)), alpha);
      st2 = pow(net(triples(i, 0), triples(i, 1)), alpha) *
        pow(net(triples(i, 2), triples(i, 1)), alpha) *
        pow(net(triples(i, 2), triples(i, 0)), alpha);
      st3 = pow(net(triples(i, 0), triples(i, 1)), alpha) *
        pow(net(triples(i, 2),triples(i, 1)), alpha) *
        pow(net(triples(i, 0), triples(i, 2)), alpha);
      st4 = pow(net(triples(i, 1), triples(i, 0)), alpha) *
        pow(net(triples(i, 1), triples(i, 2)), alpha) *
        pow(net(triples(i, 2), triples(i, 0)), alpha);
      st5 = pow(net(triples(i, 1), triples(i, 0)), alpha) *
        pow(net(triples(i, 1), triples(i, 2)), alpha) *
        pow(net(triples(i, 0), triples(i, 2)), alpha);
      st6 = pow(net(triples(i, 1), triples(i, 0)), alpha) *
        pow(net(triples(i, 2), triples(i, 1)), alpha) *
        pow(net(triples(i, 2), triples(i, 0)), alpha);
      to_return[i] = st1 + st2 + st3 + st4 + st5 + st6;
    }
  }
  return to_return;
};

// get triad weights for resampling triads for approximate MH
arma::vec triad_weights (
    arma::mat net,
    arma::Mat<double> triples,
    double alpha,
    int together,
    double smoothing_parameter) {

  int number_of_triples = triples.n_rows;

  arma::vec triad_values = gergm::indiviual_triad_values(net,
                                                         triples,
                                                         alpha,
                                                         together);

  double mean_val = arma::accu(triad_values)/double(number_of_triples);

  // loop over all values and calculate the absolute difference
  for (int i = 0; i < number_of_triples; ++i) {
    triad_values[i] = mean_val - triad_values[i];
  }
  // make differences absolute
  triad_values = arma::abs(triad_values);

  //get the maximum difference and multiply by smoothing_parameter > 1
  double max_diff = arma::max(triad_values) * smoothing_parameter;

  // get the un-normalized weights
  for (int i = 0; i < number_of_triples; ++i) {
    triad_values[i] = max_diff - triad_values[i];
  }

  // normalize
  double sum = arma::accu(triad_values);
  for (int i = 0; i < number_of_triples; ++i) {
    triad_values[i] = triad_values[i]/sum;
  }

  return triad_values;
}


// Function that will calculate h statistics
double calculate_individual_statistic(arma::mat current_network,
                                  int base_statistic_index,
                                  arma::Mat<double> triples,
                                  arma::Mat<double> pairs,
                                  double alpha,
                                  int together,
                                  arma::uvec use_selected_rows,
                                  bool use_all_rows) {

  double to_return = 0;

  if (base_statistic_index == 0) {
    if (!use_all_rows) {
      arma::mat pass_in = triples.rows(use_selected_rows);
      to_return = Out2Star(current_network,
                           pass_in,
                           alpha,
                           together);
    } else {
      // the base case for normal statistics
      to_return = Out2Star(current_network,
                           triples,
                           alpha,
                           together);
    }
  } // end of Out 2 stars conditional
  if (base_statistic_index == 1) {
    if (!use_all_rows) {
      arma::mat pass_in = triples.rows(use_selected_rows);
      to_return = In2Star(current_network,
                           pass_in,
                           alpha,
                           together);
    } else {
      // the base case for normal statistics
      to_return = In2Star(current_network,
                           triples,
                           alpha,
                           together);
    }
  }
  if (base_statistic_index == 2) {
    if (!use_all_rows) {
      arma::mat pass_in = triples.rows(use_selected_rows);
      to_return = CTriads(current_network,
                          pass_in,
                          alpha,
                          together);
    } else {
      // the base case for normal statistics
      to_return = CTriads(current_network,
                          triples,
                          alpha,
                          together);
    }
  }
  if (base_statistic_index == 3) {
    if (!use_all_rows) {
      arma::mat pass_in = pairs.rows(use_selected_rows);
      to_return = Recip(current_network,
                          pass_in,
                          alpha,
                          together);
    } else {
      // the base case for normal statistics
      to_return = Recip(current_network,
                          pairs,
                          alpha,
                          together);
    }
  }
  if (base_statistic_index == 4) {
    if (!use_all_rows) {
      arma::mat pass_in = triples.rows(use_selected_rows);
      to_return = TTriads(current_network,
                          pass_in,
                          alpha,
                          together);
    } else {
      // the base case for normal statistics
      to_return = TTriads(current_network,
                          triples,
                          alpha,
                          together);
    }
  }
  if (base_statistic_index == 5) {
    if (!use_all_rows) {
      arma::mat pass_in = pairs.rows(use_selected_rows);
      to_return = EdgeDensity(current_network,
                          pass_in,
                          alpha,
                          together);
    } else {
      // the base case for normal statistics
      to_return = EdgeDensity(current_network,
                              pairs,
                              alpha,
                              together);
    }
  }

  return to_return;
}


// Function to set up h statistic calculation
double get_individual_statistic_value(arma::mat current_network,
                                      arma::vec statistics_to_use,
                                      int index,
                                      arma::Mat<double> triples,
                                      arma::Mat<double> pairs,
                                      double alpha,
                                      int together,
                                      arma::umat selected_rows_matrix,
                                      arma::vec rows_to_use,
                                      arma::vec non_base_statistic_indicator,
                                      arma::Mat<double> random_triad_samples,
                                      arma::Mat<double> random_dyad_samples,
                                      bool use_triad_sampling) {

    // some notes on particular arguments:
    //
    // use_selected_rows -- for statistics we wish to calculate on some subnetwork,
    // this matrix will provide the row indices. The matrix will have one column
    // for each statistic to be included in the model and will have a number of
    // rows equal to the maximum number of indicies indicated for any particular
    // statistic.
    //
    // rows_to_use -- tells us how many rows in each column shoul actually be used.
    // this is meant to avoid having to pass in a list which would be horribly slow.
    //
    // statistics_to_use -- a vector for each statistic included in the model
    // indicating which base statistic it corresponds to.
    //
    // index  -- the numeric index of which entry in statistics_to_use we are
    // operating on

    // get the current statistic index
    int base_statistic_index = statistics_to_use[index];

    // determine whether we are going to use all rows or not
    bool use_all_rows = true;

    // set the use_selected_rows vector even if we are not going to use them.
    arma::uvec use_selected_rows = selected_rows_matrix.col(index);
    int from = 0;
    int to = rows_to_use[index];
    //Rcpp::Rcout << "from " << from << std::endl;
    //Rcpp::Rcout << "to " << to << std::endl;
    use_selected_rows = use_selected_rows.subvec(from, to);
    //Rcpp::Rcout << "use_selected_rows " << use_selected_rows << std::endl;

    double to_return = 0;
    // if we are not using a base statistic, set use_all_rows to false.
    if (non_base_statistic_indicator[index] == 1) {
      use_all_rows = false;
      to_return = calculate_individual_statistic(current_network,
                                                 base_statistic_index,
                                                 triples,
                                                 pairs,
                                                 alpha,
                                                 together,
                                                 use_selected_rows,
                                                 use_all_rows);
    } else {
      // if we are using a base statistic and are using triad down-sampling then
      // set use selected rows to this.
      if (use_triad_sampling) {
        // Function that will calculate h statistics
        to_return = calculate_individual_statistic(current_network,
                                                   base_statistic_index,
                                                   random_triad_samples,
                                                   random_dyad_samples,
                                                   alpha,
                                                   together,
                                                   use_selected_rows,
                                                   use_all_rows);
      } else {
        // Function that will calculate h statistics
        to_return = calculate_individual_statistic(current_network,
                                                   base_statistic_index,
                                                   triples,
                                                   pairs,
                                                   alpha,
                                                   together,
                                                   use_selected_rows,
                                                   use_all_rows);
      }
    }

    return to_return;
  }


struct Parallel_CalculateNetworkStatistics : public RcppParallel::Worker {
  arma::mat current_network;
  arma::vec statistics_to_use;
  arma::Mat<double> triples;
  arma::Mat<double> pairs;
  arma::vec alphas;
  int together;
  arma::umat selected_rows_matrix;
  arma::vec rows_to_use;
  arma::vec non_base_statistic_indicator;
  arma::Mat<double> random_triad_samples;
  arma::Mat<double> random_dyad_samples;
  bool use_triad_sampling;
  arma::vec thetas;
  RcppParallel::RVector<double> return_dist;

  Parallel_CalculateNetworkStatistics(arma::mat current_network,
                                      arma::vec statistics_to_use,
                                      arma::Mat<double> triples,
                                      arma::Mat<double> pairs,
                                      arma::vec alphas,
                                      int together,
                                      arma::umat selected_rows_matrix,
                                      arma::vec rows_to_use,
                                      arma::vec non_base_statistic_indicator,
                                      arma::Mat<double> random_triad_samples,
                                      arma::Mat<double> random_dyad_samples,
                                      bool use_triad_sampling,
                                      arma::vec thetas,
                                      Rcpp::NumericVector return_dist)
    : current_network(current_network),
      statistics_to_use(statistics_to_use),
      triples(triples),
      pairs(pairs),
      alphas(alphas),
      together(together),
      selected_rows_matrix(selected_rows_matrix),
      rows_to_use(rows_to_use),
      non_base_statistic_indicator(non_base_statistic_indicator),
      random_triad_samples(random_triad_samples),
      random_dyad_samples(random_dyad_samples),
      use_triad_sampling(use_triad_sampling),
      thetas(thetas),
      return_dist(return_dist) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      return_dist[i] = thetas[i] * get_individual_statistic_value(
        current_network,
        statistics_to_use,
        i,
        triples,
        pairs,
        alphas[i],
        together,
        selected_rows_matrix,
        rows_to_use,
        non_base_statistic_indicator,
        random_triad_samples,
        random_dyad_samples,
        use_triad_sampling);
    }
  }
};

double parallel_CalculateNetworkStatistics(
    arma::mat current_network,
    arma::vec statistics_to_use,
    arma::Mat<double> triples,
    arma::Mat<double> pairs,
    arma::vec alphas,
    int together,
    arma::umat selected_rows_matrix,
    arma::vec rows_to_use,
    arma::vec non_base_statistic_indicator,
    arma::Mat<double> random_triad_samples,
    arma::Mat<double> random_dyad_samples,
    bool use_triad_sampling,
    arma::vec thetas) {

  int number_of_stats = statistics_to_use.n_elem;
  Rcpp::NumericVector output_vec(number_of_stats);

  Parallel_CalculateNetworkStatistics Parallel_CalculateNetworkStatistics(
      current_network,
      statistics_to_use,
      triples,
      pairs,
      alphas,
      together,
      selected_rows_matrix,
      rows_to_use,
      non_base_statistic_indicator,
      random_triad_samples,
      random_dyad_samples,
      use_triad_sampling,
      thetas,
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
                                  arma::Mat<double> triples,
                                  arma::Mat<double> pairs,
                                  arma::vec alphas,
                                  int together,
                                  bool parallel,
                                  arma::umat selected_rows_matrix,
                                  arma::vec rows_to_use,
                                  arma::vec non_base_statistic_indicator,
                                  arma::Mat<double> random_triad_samples,
                                  arma::Mat<double> random_dyad_samples,
                                  bool use_triad_sampling) {

  // this is the number of statistics we will be operating on in sampling
  int number_of_thetas = statistics_to_use.n_elem;

  double to_return = 0;
  if (parallel) {
    to_return = gergm::parallel_CalculateNetworkStatistics(
      current_network,
      statistics_to_use,
      triples,
      pairs,
      alphas,
      together,
      selected_rows_matrix,
      rows_to_use,
      non_base_statistic_indicator,
      random_triad_samples,
      random_dyad_samples,
      use_triad_sampling,
      thetas);
  } else {
    for (int i = 0; i < number_of_thetas; ++i) {
      to_return += thetas[i] * get_individual_statistic_value(
          current_network,
          statistics_to_use,
          i,
          triples,
          pairs,
          alphas[i],
          together,
          selected_rows_matrix,
          rows_to_use,
          non_base_statistic_indicator,
          random_triad_samples,
          random_dyad_samples,
          use_triad_sampling);
    }
  }
  return to_return;
};

  // Function that will calculate and save all of the h statistics for a network
  arma::vec save_network_statistics(arma::mat current_network,
                                    arma::vec statistics_to_use,
                                    arma::vec base_statistics_to_save,
                                    arma::vec base_statistic_alphas,
                                    arma::Mat<double> triples,
                                    arma::Mat<double> pairs,
                                    arma::vec alphas,
                                    int together,
                                    arma::umat save_statistics_selected_rows_matrix,
                                    arma::vec rows_to_use,
                                    int num_non_base_statistics,
                                    arma::vec non_base_statistic_indicator) {

    // determine which statistics are non-base statistics, then append these on to
    // the end of the vector of all base statistics that we are going to calculate
    // and return for this
    int num_base_statistics_to_save = base_statistics_to_save.n_elem;
    int statistics_to_save = num_non_base_statistics + num_base_statistics_to_save;
    //combine everything appropriately so we can calculate on the extended set of
    //statistics.
    arma::vec statistic_values = arma::zeros(statistics_to_save);

    // take all of the base_statistics_to_save values, then add on the non-base
    // statistics at the end from the statistics_to_use vector we use in the
    // normal MH updates.
    arma::vec combined_statistics_to_use = arma::zeros(statistics_to_save);
    arma::vec combined_alphas = arma::zeros(statistics_to_save);
    arma::vec combined_rows_to_use = arma::zeros(statistics_to_save);
    arma::vec combined_non_base_statistic_indicator = arma::zeros(statistics_to_save);

    // loop through and populate with base statistic values
    for (int i = 0; i < num_base_statistics_to_save; ++i) {
      combined_statistics_to_use[i] = base_statistics_to_save[i];
      combined_alphas[i] = base_statistic_alphas[i];
    }
    int counter = num_base_statistics_to_save;
    for (int i = 0; i < non_base_statistic_indicator.n_elem; ++i) {
      if (non_base_statistic_indicator[i] == 1) {
        combined_statistics_to_use[counter] = statistics_to_use[i];
        combined_alphas[counter] = alphas[i];
        combined_rows_to_use[counter] = rows_to_use[i];
        combined_non_base_statistic_indicator[counter] = 1;
        counter += 1;
      }
    }

    // figure out the number of rows to use
    arma::Mat<double> proxy_random_triad_samples(2,statistics_to_save);
    arma::Mat<double> proxy_random_dyad_samples(2,statistics_to_save);

    // loop through and calculate the statistics we are going to save and store
    // them in a vector.
    for (int i = 0; i < statistics_to_save; ++i) {
      statistic_values[i] = get_individual_statistic_value(
        current_network,
        combined_statistics_to_use,
        i,
        triples,
        pairs,
        combined_alphas[i],
        together,
        save_statistics_selected_rows_matrix,
        combined_rows_to_use,
        combined_non_base_statistic_indicator,
        proxy_random_triad_samples,
        proxy_random_dyad_samples,
        false);
    }


    return statistic_values;
  };

  // Function that will calculate h statistics
  double integrand(arma::mat current_network,
                   arma::vec statistics_to_use,
                   arma::vec thetas,
                   arma::Mat<double> triples,
                   arma::Mat<double> pairs,
                   arma::umat save_statistics_selected_rows_matrix,
                   arma::vec rows_to_use,
                   arma::vec base_statistics_to_save,
                   arma::vec base_statistic_alphas,
                   int num_non_base_statistics,
                   arma::vec non_base_statistic_indicator,
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
    arma::vec save_stats = gergm::save_network_statistics(
      current_network,
      statistics_to_use,
      base_statistics_to_save,
      base_statistic_alphas,
      triples,
      pairs,
      alphas,
      together,
      save_statistics_selected_rows_matrix,
      rows_to_use,
      num_non_base_statistics,
      non_base_statistic_indicator);

    // now multiply thetas by the statistics
    int num_theta = thetas.n_elem;
    double to_return = 0;
    for (int i = 0; i < num_theta; ++i) {
      to_return += thetas[i] * save_stats[i];
    }
    return to_return;
  };


  // create a RcppParallel::Worker struct that we can use to fill in the
  // entries in our token topic distribution in parallel
  struct Parallel_Integrand : public RcppParallel::Worker {

    // Instantiate all of our input variables which will then be initialized
    // via the Function initializer below. This is how the struct makes sure
    // that everything has the right type
    arma::mat current_network;
    arma::vec statistics_to_use;
    arma::vec thetas;
    arma::Mat<double> triples;
    arma::Mat<double> pairs;
    arma::umat save_statistics_selected_rows_matrix;
    arma::vec rows_to_use;
    arma::vec base_statistics_to_save;
    arma::vec base_statistic_alphas;
    int num_non_base_statistics;
    arma::vec non_base_statistic_indicator;
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
                       arma::Mat<double> triples,
                       arma::Mat<double> pairs,
                       arma::umat save_statistics_selected_rows_matrix,
                       arma::vec rows_to_use,
                       arma::vec base_statistics_to_save,
                       arma::vec base_statistic_alphas,
                       int num_non_base_statistics,
                       arma::vec non_base_statistic_indicator,
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
        save_statistics_selected_rows_matrix(save_statistics_selected_rows_matrix),
        rows_to_use(rows_to_use),
        base_statistics_to_save(base_statistics_to_save),
        base_statistic_alphas(base_statistic_alphas),
        num_non_base_statistics(num_non_base_statistics),
        non_base_statistic_indicator(non_base_statistic_indicator),
        alphas(alphas),
        together(together),
        sender(sender),
        recipient(recipient),
        integration_interval(integration_interval),
        return_dist(return_dist) {}

    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
        double integ = gergm::integrand(current_network,
                                        statistics_to_use,
                                        thetas,
                                        triples,
                                        pairs,
                                        save_statistics_selected_rows_matrix,
                                        rows_to_use,
                                        base_statistics_to_save,
                                        base_statistic_alphas,
                                        num_non_base_statistics,
                                        non_base_statistic_indicator,
                                        alphas,
                                        together,
                                        sender,
                                        recipient,
                                        integration_interval[i]);
        return_dist[i] = integ;
      }
    }
  };


  arma::vec parallel_integration(
      arma::mat current_network,
      arma::vec statistics_to_use,
      arma::vec thetas,
      arma::Mat<double> triples,
      arma::Mat<double> pairs,
      arma::umat save_statistics_selected_rows_matrix,
      arma::vec rows_to_use,
      arma::vec base_statistics_to_save,
      arma::vec base_statistic_alphas,
      int num_non_base_statistics,
      arma::vec non_base_statistic_indicator,
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
        save_statistics_selected_rows_matrix,
        rows_to_use,
        base_statistics_to_save,
        base_statistic_alphas,
        num_non_base_statistics,
        non_base_statistic_indicator,
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
  };

  double log_sum_exp_integrator (arma::mat current_network,
                                 arma::vec statistics_to_use,
                                 arma::vec thetas,
                                 arma::Mat<double> triples,
                                 arma::Mat<double> pairs,
                                 arma::umat save_statistics_selected_rows_matrix,
                                 arma::vec rows_to_use,
                                 arma::vec base_statistics_to_save,
                                 arma::vec base_statistic_alphas,
                                 int num_non_base_statistics,
                                 arma::vec non_base_statistic_indicator,
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
        save_statistics_selected_rows_matrix,
        rows_to_use,
        base_statistics_to_save,
        base_statistic_alphas,
        num_non_base_statistics,
        non_base_statistic_indicator,
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
                                            save_statistics_selected_rows_matrix,
                                            rows_to_use,
                                            base_statistics_to_save,
                                            base_statistic_alphas,
                                            num_non_base_statistics,
                                            non_base_statistic_indicator,
                                            alphas,
                                            together,
                                            sender,
                                            recipient,
                                            integration_interval[i]);
      }
    }

    //Rcpp::Rcout << "integral_evaluations " << integral_evaluations << std::endl;

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

  };


} //end of gergm namespace

using std::log;

// [[Rcpp::export]]
List Extended_Metropolis_Hastings_Sampler (int number_of_iterations,
                                  double shape_parameter,
                                  int number_of_nodes,
                                  arma::vec statistics_to_use,
                                  arma::mat initial_network,
                                  int take_sample_every,
                                  arma::vec thetas,
                                  arma::Mat<double> triples,
                                  arma::Mat<double> pairs,
                                  arma::vec alphas,
                                  int together,
                                  int seed,
                                  int number_of_samples_to_store,
                                  int undirect_network,
                                  bool parallel,
                                  arma::umat use_selected_rows,
                                  arma::umat save_statistics_selected_rows_matrix,
                                  arma::vec rows_to_use,
                                  arma::vec base_statistics_to_save,
                                  arma::vec base_statistic_alphas,
                                  int num_non_base_statistics,
                                  arma::vec non_base_statistic_indicator,
                                  double p_ratio_multaplicative_factor,
                                  Rcpp::List random_triad_sample_list,
                                  Rcpp::List random_dyad_sample_list,
                                  bool use_triad_sampling,
                                  int num_unique_random_triad_samples) {

  // Allocate variables and data structures
  double variance = shape_parameter;
  // the list we will put stuff in to return it to R
  int list_length = 9;
  List to_return(list_length);
  // this is the number of statistics we will be saving (all selected base + non base)
  int statistics_to_save = num_non_base_statistics +
    base_statistics_to_save.n_elem;

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
                                             statistics_to_save);
  arma::mat current_edge_weights = initial_network;
  arma::mat corr_current_edge_weights = arma::zeros (number_of_nodes, number_of_nodes);

  // values for stochastic MH
  arma::Mat<double> random_triad_samples(2,2);
  arma::Mat<double> random_dyad_samples(2,2);
  int update_triad_samples_every = 10;
  int triad_sample_update_counter = 0;
  int random_triad_sample_counter = 0;
  if (use_triad_sampling) {
    arma::Mat<double> temp = random_triad_sample_list[random_triad_sample_counter];
    random_triad_samples = temp;
    arma::Mat<double> temp2 = random_dyad_sample_list[random_triad_sample_counter];
    random_dyad_samples = temp2;
    random_triad_sample_counter += 1;
  }


  // Set RNG and define uniform distribution
  boost::mt19937 generator(seed);
  boost::uniform_01<double> uniform_distribution;
  // Outer loop over the number of samples
  for (int n = 0; n < number_of_iterations; ++n) {
    //Rcpp::Rcout << "Iteration: " << n << std::endl;
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
            gergm::normal_distribution<double> proposal(current_edge_value,variance);
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
            gergm::normal_distribution<double> proposal(current_edge_value,variance);
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
    } //end of condition for whether we are using a correlation network

    double proposed_addition = 0;
    double current_addition = 0;

    // if we are using random triad sampling, then extract the appropriate
    // triples and pairs matrices
    if (update_triad_samples_every == triad_sample_update_counter) {
      if (use_triad_sampling) {
        arma::Mat<double> temp = random_triad_sample_list[random_triad_sample_counter];
        random_triad_samples = temp;
        arma::Mat<double> temp2 = random_dyad_sample_list[random_triad_sample_counter];
        random_dyad_samples = temp2;
      }
      triad_sample_update_counter = 0;
      // increment the random sample counter
      random_triad_sample_counter += 1;
      // if it is equal to the number of different samples we are working with,
      // then reset it
      if (num_unique_random_triad_samples == random_triad_sample_counter) {
        random_triad_sample_counter = 0;
      }
    }
    triad_sample_update_counter += 1;

      proposed_addition = gergm::CalculateNetworkStatistics(
        proposed_edge_weights,
        statistics_to_use,
        thetas,
        triples,
        pairs,
        alphas,
        together,
        parallel,
        use_selected_rows,
        rows_to_use,
        non_base_statistic_indicator,
        random_triad_samples,
        random_dyad_samples,
        use_triad_sampling);
      // only calculate the h function if we updated the network last round
      // otherwise use the cached value.
      if (network_did_not_change) {
        current_addition = previous_h_function_value;
      } else {
        current_addition = gergm::CalculateNetworkStatistics(
          current_edge_weights,
          statistics_to_use,
          thetas,
          triples,
          pairs,
          alphas,
          together,
          parallel,
          use_selected_rows,
          rows_to_use,
          non_base_statistic_indicator,
          random_triad_samples,
          random_dyad_samples,
          use_triad_sampling);
        previous_h_function_value = current_addition ;
      }
    

    // store some additional diagnostics h value is the last entry
    P_Ratios[n] = p_ratio_multaplicative_factor * (proposed_addition -
      current_addition);
    Q_Ratios[n] = log_prob_accept;

    double total_edges = double(number_of_nodes * (number_of_nodes - 1));
    double temp1 = arma::accu(proposed_edge_weights);
    Proposed_Density[n] = temp1/total_edges;
    double temp2 = arma::accu(current_edge_weights);
    Current_Density[n] = temp2/total_edges;

    // now we add in a p-ratio multaplicative factor incase we are randomly
    // downsampling
    log_prob_accept += p_ratio_multaplicative_factor * (proposed_addition -
      current_addition);

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
      //Rcpp::Rcout << "Iteration: " << n << std::endl;
        arma::vec save_stats = gergm::save_network_statistics(
          current_edge_weights,
          statistics_to_use,
          base_statistics_to_save,
          base_statistic_alphas,
          triples,
          pairs,
          alphas,
          together,
          save_statistics_selected_rows_matrix,
          rows_to_use,
          num_non_base_statistics,
          non_base_statistic_indicator);
        for (int m = 0; m < statistics_to_save; ++m) {
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
  to_return[4] = Log_Prob_Accept;
  to_return[5] = P_Ratios;
  to_return[6] = Q_Ratios;
  to_return[7] = Proposed_Density;
  to_return[8] = Current_Density;
  return to_return;
}


// [[Rcpp::export]]
arma::vec h_statistics (arma::vec statistics_to_use,
                       arma::mat current_edge_weights,
                       arma::Mat<double> triples,
                       arma::Mat<double> pairs,
                       arma::vec alphas,
                       int together,
                       arma::umat save_statistics_selected_rows_matrix,
                       arma::vec rows_to_use,
                       arma::vec base_statistics_to_save,
                       arma::vec base_statistic_alphas,
                       int num_non_base_statistics,
                       arma::vec non_base_statistic_indicator) {

    arma::vec save_stats = gergm::save_network_statistics(
      current_edge_weights,
      statistics_to_use,
      base_statistics_to_save,
      base_statistic_alphas,
      triples,
      pairs,
      alphas,
      together,
      save_statistics_selected_rows_matrix,
      rows_to_use,
      num_non_base_statistics,
      non_base_statistic_indicator);

  return save_stats;
}


// [[Rcpp::export]]
double extended_weighted_mple_objective ( int number_of_nodes,
                                 arma::vec statistics_to_use,
                                 arma::mat current_network,
                                 arma::Mat<double> triples,
                                 arma::Mat<double> pairs,
                                 arma::umat save_statistics_selected_rows_matrix,
                                 arma::vec rows_to_use,
                                 arma::vec base_statistics_to_save,
                                 arma::vec base_statistic_alphas,
                                 int num_non_base_statistics,
                                 arma::vec non_base_statistic_indicator,
                                 arma::vec thetas,
                                 arma::vec alphas,
                                 int together,
                                 arma::vec integration_interval,
                                 bool parallel) {

  double objective = 0;
  // Calculate the
  for (int i = 0; i < number_of_nodes; ++i) {
    for (int j = 0; j < number_of_nodes; ++j) {
      double temp1 = gergm::integrand(current_network,
                                     statistics_to_use,
                                     thetas,
                                     triples,
                                     pairs,
                                     save_statistics_selected_rows_matrix,
                                     rows_to_use,
                                     base_statistics_to_save,
                                     base_statistic_alphas,
                                     num_non_base_statistics,
                                     non_base_statistic_indicator,
                                     alphas,
                                     together,
                                     i,
                                     j,
                                     -1);

      double temp2 = gergm::log_sum_exp_integrator(current_network,
                                                  statistics_to_use,
                                                  thetas,
                                                  triples,
                                                  pairs,
                                                  save_statistics_selected_rows_matrix,
                                                  rows_to_use,
                                                  base_statistics_to_save,
                                                  base_statistic_alphas,
                                                  num_non_base_statistics,
                                                  non_base_statistic_indicator,
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


// [[Rcpp::export]]
arma::vec get_indiviual_triad_values (
                        arma::mat net,
                        arma::Mat<double> triples,
                        double alpha,
                        int together) {

arma::vec triad_values = gergm::indiviual_triad_values(net,
                                 triples,
                                 alpha,
                                 together);

  return triad_values;
}

// [[Rcpp::export]]
arma::vec get_triad_weights (
    arma::mat net,
    arma::Mat<double> triples,
    double alpha,
    int together,
    double smoothing_parameter) {

  arma::vec triad_weights = gergm::triad_weights (
    net,
    triples,
    alpha,
    together,
    smoothing_parameter);

  return triad_weights;
}


// [[Rcpp::export]]
List Individual_Edge_Conditional_Prediction (
    int number_of_iterations,
    double shape_parameter,
    int number_of_nodes,
    arma::vec statistics_to_use,
    arma::mat initial_network,
    int take_sample_every,
    arma::vec thetas,
    arma::Mat<double> triples,
    arma::Mat<double> pairs,
    arma::vec alphas,
    int together,
    int seed,
    int number_of_samples_to_store,
    int undirect_network,
    bool parallel,
    arma::umat use_selected_rows,
    arma::umat save_statistics_selected_rows_matrix,
    arma::vec rows_to_use,
    arma::vec base_statistics_to_save,
    arma::vec base_statistic_alphas,
    int num_non_base_statistics,
    arma::vec non_base_statistic_indicator,
    double p_ratio_multaplicative_factor,
    Rcpp::List random_triad_sample_list,
    Rcpp::List random_dyad_sample_list,
    bool use_triad_sampling,
    int num_unique_random_triad_samples,
    int i,
    int j) {

  // Allocate variables and data structures
  double variance = shape_parameter;
  // the list we will put stuff in to return it to R
  int list_length = 9;
  List to_return(list_length);
  // this is the number of statistics we will be saving (all selected base + non base)
  int statistics_to_save = num_non_base_statistics +
    base_statistics_to_save.n_elem;

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
                                             statistics_to_save);
  arma::mat current_edge_weights = initial_network;
  arma::mat corr_current_edge_weights = arma::zeros (number_of_nodes, number_of_nodes);

  // values for stochastic MH
  arma::Mat<double> random_triad_samples(2,2);
  arma::Mat<double> random_dyad_samples(2,2);
  int update_triad_samples_every = 10;
  int triad_sample_update_counter = 0;
  int random_triad_sample_counter = 0;
  if (use_triad_sampling) {
    arma::Mat<double> temp = random_triad_sample_list[random_triad_sample_counter];
    random_triad_samples = temp;
    arma::Mat<double> temp2 = random_dyad_sample_list[random_triad_sample_counter];
    random_dyad_samples = temp2;
    random_triad_sample_counter += 1;
  }

  // Set RNG and define uniform distribution
  boost::mt19937 generator(seed);
  boost::uniform_01<double> uniform_distribution;
  // Outer loop over the number of samples
  for (int n = 0; n < number_of_iterations; ++n) {
    //Rcpp::Rcout << "Iteration: " << n << std::endl;
    double log_prob_accept = 0;
    arma::mat proposed_edge_weights = current_edge_weights;

    // deal with the case where we have an undirected network.
    if(undirect_network == 1){
      // Run loop to sample new edge weights

      double log_probability_of_current_under_new = 0;
      double log_probability_of_new_under_current = 0;
      //draw a new edge value centered at the old edge value
      double current_edge_value = current_edge_weights(i,j);
      //draw from a truncated normal
      gergm::normal_distribution<double> proposal(current_edge_value,variance);
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
      //save everything
      proposed_edge_weights(i,j) = new_edge_value;
      proposed_edge_weights(j,i) = new_edge_value;
      log_probability_of_new_under_current = log(prob_new_edge_under_old);
      log_probability_of_current_under_new = log(prob_old_edge_under_new);

      // Calculate acceptance probability
      log_prob_accept += (log_probability_of_current_under_new
                            - log_probability_of_new_under_current);

    }else{
      // int counter =  0;
      // Run loop to sample new edge weights
      double log_probability_of_current_under_new = 0;
      double log_probability_of_new_under_current = 0;
      //draw a new edge value centered at the old edge value
      double current_edge_value = current_edge_weights(i,j);
      //draw from a truncated normal
      gergm::normal_distribution<double> proposal(current_edge_value,variance);
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

      //save everything
      proposed_edge_weights(i,j) = new_edge_value;
      log_probability_of_new_under_current = log(prob_new_edge_under_old);
      log_probability_of_current_under_new = log(prob_old_edge_under_new);

      // Calculate acceptance probability
      log_prob_accept += (log_probability_of_current_under_new
                            - log_probability_of_new_under_current);
    } //end of condition for whether we are using a correlation network

    double proposed_addition = 0;
    double current_addition = 0;

    // if we are using random triad sampling, then extract the appropriate
    // triples and pairs matrices
    if (update_triad_samples_every == triad_sample_update_counter) {
      if (use_triad_sampling) {
        arma::Mat<double> temp = random_triad_sample_list[random_triad_sample_counter];
        random_triad_samples = temp;
        arma::Mat<double> temp2 = random_dyad_sample_list[random_triad_sample_counter];
        random_dyad_samples = temp2;
      }
      triad_sample_update_counter = 0;
      // increment the random sample counter
      random_triad_sample_counter += 1;
      // if it is equal to the number of different samples we are working with,
      // then reset it
      if (num_unique_random_triad_samples == random_triad_sample_counter) {
        random_triad_sample_counter = 0;
      }
    }
    triad_sample_update_counter += 1;

      proposed_addition = gergm::CalculateNetworkStatistics(
        proposed_edge_weights,
        statistics_to_use,
        thetas,
        triples,
        pairs,
        alphas,
        together,
        parallel,
        use_selected_rows,
        rows_to_use,
        non_base_statistic_indicator,
        random_triad_samples,
        random_dyad_samples,
        use_triad_sampling);
      // only calculate the h function if we updated the network last round
      // otherwise use the cached value.
      if (network_did_not_change) {
        current_addition = previous_h_function_value;
      } else {
        current_addition = gergm::CalculateNetworkStatistics(
          current_edge_weights,
          statistics_to_use,
          thetas,
          triples,
          pairs,
          alphas,
          together,
          parallel,
          use_selected_rows,
          rows_to_use,
          non_base_statistic_indicator,
          random_triad_samples,
          random_dyad_samples,
          use_triad_sampling);
        previous_h_function_value = current_addition ;
      }
    

    // store some additional diagnostics h value is the last entry
    P_Ratios[n] = p_ratio_multaplicative_factor * (proposed_addition -
      current_addition);
    Q_Ratios[n] = log_prob_accept;

    double total_edges = double(number_of_nodes * (number_of_nodes - 1));
    double temp1 = arma::accu(proposed_edge_weights);
    Proposed_Density[n] = temp1/total_edges;
    double temp2 = arma::accu(current_edge_weights);
    Current_Density[n] = temp2/total_edges;

    // now we add in a p-ratio multaplicative factor incase we are randomly
    // downsampling
    log_prob_accept += p_ratio_multaplicative_factor * (proposed_addition -
      current_addition);

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
      //Rcpp::Rcout << "Iteration: " << n << std::endl;
        arma::vec save_stats = gergm::save_network_statistics(
          current_edge_weights,
          statistics_to_use,
          base_statistics_to_save,
          base_statistic_alphas,
          triples,
          pairs,
          alphas,
          together,
          save_statistics_selected_rows_matrix,
          rows_to_use,
          num_non_base_statistics,
          non_base_statistic_indicator);
        for (int m = 0; m < statistics_to_save; ++m) {
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
  to_return[4] = Log_Prob_Accept;
  to_return[5] = P_Ratios;
  to_return[6] = Q_Ratios;
  to_return[7] = Proposed_Density;
  to_return[8] = Current_Density;
  return to_return;
}
