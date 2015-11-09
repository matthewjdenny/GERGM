// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
using std::pow;
using std::sqrt;
using namespace Rcpp;

// [[Rcpp::export]]
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
      arma::mat R2_solved = R2.i();

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

// namespace mjd {
//
// }

// // [[Rcpp::export]]
// double speed_test2(arma::mat partial_correlations){
//   double temp = 0;
//   for (int k = 0; k < 1000000; ++k) {
//     temp = mjd::jacobian(partial_correlations);
//   }
//   return temp;
// }
