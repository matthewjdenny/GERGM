// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
using std::pow;
using namespace Rcpp;

namespace mjd {
  double jacobian(arma::mat partial_correlations){
    int nrow = partial_correlations.n_rows;
    arma::vec corrs_1 = partial_correlations.diag(1);
    arma::vec temp = pow(corrs_1,2);
    arma::vec temp2 = pow((1 - temp),(nrow -2));
    double prod_1 =  arma::prod(temp2);
    double prod_2 = 1;

    for (int k = 2; k < (nrow -1); ++k) {
      for (int i = 0; i < (nrow - k); ++i) {
        //Rcpp::Rcout << i+1 << ","<< i+k+1 << std::endl;
        double temp3 = pow(partial_correlations(i,(i + k)),2);
        prod_2 = prod_2 * pow((1 - temp3),(nrow - 1 - k));
      }
    }

    double temp4 = pow(prod_1,(nrow - 2));
    double result = pow(temp4*prod_2,0.5);
    return result;
  }
}

// // [[Rcpp::export]]
// double speed_test(arma::mat partial_correlations){
//   double temp = 0;
//   for (int k = 0; k < 1000000; ++k) {
//     temp = mjd::jacobian(partial_correlations);
//   }
//   return temp;
// }
