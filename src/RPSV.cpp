#include "ipsecr.h"

//==============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector rpsvcpp (
        const Rcpp::IntegerMatrix &sk,  
        const Rcpp::NumericMatrix &traps) {
    
    int S = sk.nrow();
    int K = sk.ncol();
    int n = 0;
    int i,s,k,cx;
    double sumx = 0;
    double sumy = 0;
    double sumx2 = 0;
    double sumy2 = 0;
    
    for (s=0; s<S; s++) {
        for (k=0; k<K; k++) {
            cx = sk(s,k);
            if (cx>0) {
                n += cx;
                for (i=0; i<cx; i++) {
                    sumx  += traps(k,0); 
                    sumy  += traps(k,1); 
                    sumx2 += traps(k,0) * traps(k,0); 
                    sumy2 += traps(k,1) * traps(k,1);
                }
            }
        }
    }
    Rcpp::NumericVector stats(3);
    stats[0] = n-1;
    stats[1] = sumx2 - sumx*sumx/n;
    stats[2] = sumy2 - sumy*sumy/n;
    
    return (stats);
}
//==============================================================================
