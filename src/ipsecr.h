// include guard
#ifndef __ipsecr_h_INCLUDED__   // if ipsecr.h hasn't been included yet...
#define __ipsecr_h_INCLUDED__   // #define this so the compiler knows it has been included

//------------------------------------------------------------------------------

// BOOST used for statistical distributions
// return NAN for invalid inputs
// see https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/stat_tut/weg/error_eg.html
// and https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/pol_tutorial/changing_policy_defaults.html

#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
// must follow define domain error policy...
#include <boost/math/distributions.hpp>
//------------------------------------------------------------------------------

// RcppArmadillo 2022-08-25

// RcppArmadillo.h pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

//------------------------------------------------------------------------------

// constants
#define huge 1e10

#endif  // __ipsecr_h_INCLUDED__

