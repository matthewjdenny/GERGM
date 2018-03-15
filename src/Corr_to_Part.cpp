#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat Corr_to_Part (int d,
                        arma::mat correlations,
                        arma::mat partials) {

  for (int k = 2; k < (d + 1); ++k) {
    for (int i = 1; i < (d-k + 1); ++i) {

      arma::mat R2 = correlations.submat((i ),  (i ), (i + k - 2), (i + k - 2));
      // R2 <- correlations[(i + 1):(i + k - 1), (i + 1):(i + k - 1)]
      arma::mat r1 = correlations.submat(i-1,  (i), i-1, (i + k -2));
      // r1 <- correlations[i, (i + 1):(i + k - 1)]
      arma::mat r3 = correlations.submat((i + k-1),  (i), (i + k -1), (i + k -2));
      // r3 <- correlations[i + k, (i + 1):(i + k - 1)]

      arma::mat solved = arma::inv_sympd(R2);
      arma::mat temp = r1 * solved * r1.t();
      arma::mat temp5 = solved * r1.t();
      arma::mat temp2 = r3 * solved * r3.t();

      double D = (1 - temp(0,0)) * (1 - temp2(0,0));
      D = std::sqrt(D);

      arma::mat temp3 = r1 * solved * r3.t();
      double temp4 = (correlations(i -1, (i + k -1)) -  temp3(0,0)) / D;
      partials(i -1, (i + k - 1)) = temp4;
      partials((i + k - 1), i -1) = temp4;
    }
  }

  return partials;
}
