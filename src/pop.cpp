#include "ipsecr.h"

//==============================================================================

// https://stackoverflow.com/questions/24618370/using-rmultinom-with-rcpp
Rcpp::IntegerVector oneMultinomCall(Rcpp::NumericVector probs, int N) {
    int k = probs.size();
    Rcpp::IntegerVector ans(k);
    rmultinom(N, probs.begin(), k, ans.begin());
    return(ans);
}
//===============================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix popcpp (
        const Rcpp::NumericMatrix &mask,  // x-y coord
        Rcpp::NumericVector &prob,     // cell probability)
        double &maskspacing,
        int &N) {
    Rcpp::RNGScope scope;             // Rcpp initialise and finalise random seed 
    
    int M = mask.nrow();
    Rcpp::IntegerVector nm(M);
    Rcpp::NumericMatrix animals (N,2);
    // not to be called with N < 1
    if (N<1) Rcpp::stop ("zero population requested");
    nm = oneMultinomCall(prob, N);
    int n = 0;
    for (int m = 0; m < M; m++) {
        for (int i = 0; i < nm[m]; i++) {
            for (int j = 0; j<2; j++) {
                animals(n,j) = mask(m,j) + (R::runif(0,1)-0.5)*maskspacing;
            }
            n++;
        }
    }
    return (animals);
}
//==============================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix popevencpp (
        const Rcpp::NumericMatrix &bounds,  // x-y coord
        int &N) {
    Rcpp::RNGScope scope;             // Rcpp initialise and finalise random seed 
    // not to be called with N < 1
    if (N<1) Rcpp::stop ("zero population requested");
    Rcpp::NumericMatrix animals (N,2);
    std::fill(animals.begin(), animals.end(), Rcpp::NumericVector::get_na() ) ;

    double dx = bounds(1,0) - bounds(0,0);
    double dy = bounds(1,1) - bounds(0,1);
    double area = dx * dy;
    double spacing = std::sqrt(area/N);
    int nx = trunc(dx/spacing);
    int ny = trunc(dy/spacing);
    double marginx = (dx - (nx-1)*spacing) / 2;
    double marginy = (dy - (ny-1)*spacing) / 2;
    double x,y;
    marginx = marginx * R::runif(0,1) + bounds(0,0);
    marginy = marginy * R::runif(0,1) + bounds(0,1);
    for (int i = 0; i<nx; i++) {
        x = marginx + spacing * i;
        for (int j = 0; j<ny; j++) {
            y = marginy + spacing * j;
            animals(i*ny+j,0) = x;
            animals(i*ny+j,1) = y;
        }
    }
    return (animals);
}
//==============================================================================
