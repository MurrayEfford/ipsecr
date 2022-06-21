// include guard
#ifndef __ipsecr_h_INCLUDED__   // if ipsecr.h hasn't been included yet...
#define __ipsecr_h_INCLUDED__   // #define this so the compiler knows it has been included

//------------------------------------------------------------------------------
// BOOST used for statistical distributions from secr 4.4.6 October 2021
// return NAN for invalid inputs
// see https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/stat_tut/weg/error_eg.html
// and https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/pol_tutorial/changing_policy_defaults.html

#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
// must follow define domain error policy...
#include <boost/math/distributions.hpp>
//------------------------------------------------------------------------------

#include <R.h>       // random numbers

#include <Rcpp.h>

// constants
#define fuzz 1e-200
#define huge 1e10

//-------------------
// data structures   
//-------------------
struct rpoint {
    double x;
    double y;
};

//--------------------------------------------------------------------------

int i3 (int i, int j, int k, int ii, int jj);
int i4 (int i, int j, int k, int l, int ii, int jj, int kk);
int par3 (int fn);

//--------------------------------------------------------------------------

// probability of count with distribution specified by binomN 
double countp (int count, int binomN, double lambda);
//--------------------------------------------------------------------------

double expmin (double x);

double rcount (const int binomN, const double lambda, const double Tsk);
//---------------------------------------------------------------------

double randomtime (double p);
//---------------------------------------------------------------------

double gr (
    const int fn, 
    Rcpp::NumericVector gsb, 
    const rpoint xy, 
    const rpoint animal);
//---------------------------------------------------------------------

double hazard (double pp);

//---------------------------------------------------------------------
// 
double gpois (int count, double lambda);
double gbinom(int count, int size, double p);
double pski ( int binomN, int count, double Tski, double g, double pI);

//--------------------------------------------------------------------------

double d2cpp (
    const int k, 
    const int m, 
    const Rcpp::NumericMatrix &A1, 
    const Rcpp::NumericMatrix &A2);

// Functions to characterize detector type 

bool allcapped  (const Rcpp::IntegerVector detect);
bool allmulti (const Rcpp::IntegerVector detect);

bool anyvarying (const int nc, const int ss, const int nk, const int nmix,
                 const Rcpp::IntegerVector &PIA0);
bool anyb (
    const Rcpp::NumericMatrix &gsbval, 
    const Rcpp::NumericMatrix &gsb0val);

// miscellaneous functions

int bswitch (
    const int btype, 
    const int N, 
    const int i, 
    const int k, 
    const std::vector<int> &caughtbefore);

Rcpp::IntegerVector oneMultinomCall(Rcpp::NumericVector probs, int N);
    
#endif  // __ipsecr_h_INCLUDED__

