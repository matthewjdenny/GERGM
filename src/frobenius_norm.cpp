#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

using std::abs;
using std::sqrt;
using std::pow;

// [[Rcpp::export]]
double frobenius_norm(
        arma::mat mat1,
        arma::mat mat2
){

    int nrow = mat1.n_rows;
    int ncol = mat1.n_cols;

    double distance = 0;

    // calculate df/idf
    for(int i = 0; i < nrow; ++i){
        for(int j = 0; j < ncol; ++j){
            distance += pow(abs(mat1(i,j) - mat2(i,j)),2);
        }
    }
    distance = sqrt(distance);
    return distance;
}