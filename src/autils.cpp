#include "ipsecr.h"

using namespace std;
using namespace Rcpp;

double minimumexp = -100;

//--------------------------------------------------------------------------

double expmin (double x)
{
    if (x < minimumexp)
        return(0);
    else
        return(exp(x));
}
//--------------------------------------------------------------------------

// index to vector element corresponding to cell i,j,k in 3D array
// stored in column-major order 

int i3 (int i, int j, int k, int ii, int jj) {
    return(ii * (jj * k + j) + i);
}
//--------------------------------------------------------------------------

// index to vector element corresponding to cell i,j,k,l in 4D array
// stored in column-major order 

int i4 (int i, int j, int k, int l, int ii, int jj, int kk) {
    return (ii *(jj*(kk*l + k) + j) + i);
}
//--------------------------------------------------------------------------

double d2cpp (
        const int k,
        const int m,
        const NumericMatrix &A1,
        const NumericMatrix &A2)
    // return squared distance between two points given by row k in A1
    // and row m in A2, where A1 and A2 have respectively A1rows and A2rows
{
    return(
        (A1(k,0) - A2(m,0)) * (A1(k,0) - A2(m,0)) +
            (A1(k,1) - A2(m,1)) * (A1(k,1) - A2(m,1))
    );
}
//--------------------------------------------------------------------------

// customised dpois
double gpois (int count, double lambda)
{
    double x;
    if ((count < 0) || (count>0 && lambda <= 0)) {
        return(0);
    }
    else if (count == 0) {
        return (exp(-lambda));
    }
    else {
        boost::math::poisson_distribution<> pois(lambda);
        x = boost::math::pdf(pois, count);
        return (x);
    }
}
//--------------------------------------------------------------------------

// customised dbinom
double gbinom(int count, int size, double p)
{
    double x, q;
    int i;
    if ((count < 0) || (count > 0 && p <= 0)) {
        x = 0;
    }
    else if (count == 0) {
        q = 1 - p;
        x = q;
        for (i=1; i< size; i++) x *= q;
    }
    else {
        boost::math::binomial_distribution<> bin(size, p);
        x = boost::math::pdf(bin, count);
    }
    return (x);
}
//--------------------------------------------------------------------------

// probability of count with distribution specified by binomN
double countp (int count, int binomN, double lambda) {
    // Poisson
    if (binomN == 0) {
        if (count == 0)
            return (exp(-lambda));
        else {
            // return (R::dpois(count, lambda, 0));
            boost::math::poisson_distribution<> pois(lambda);
            return (boost::math::pdf(pois, count));
        }
    }

    // Bernoulli
    else if (binomN == 1) {
        if (count == 0)
            return ( 1 - lambda );
        else
            return ( lambda );
    }

    // negative binomial
    else if (binomN < 0) {
        boost::math::negative_binomial_distribution<> nbin(binomN, lambda);
        return (boost::math::pdf(nbin, count));
    }

    // binomial
    else {
        boost::math::binomial_distribution<> bin(binomN, lambda);
        return (boost::math::pdf(bin, count));
    }
}
//--------------------------------------------------------------------------

double zrcpp (double r, int detectfn, NumericVector par)
{
    if (detectfn == 14) {  // hazard halfnormal
        return (exp(-r*r / 2 / par(1) / par(1)));    
    }
    else {
        if (detectfn == 15) {  // hazard hazard rate
            return (1 - exp(- pow(r /par(1), - par(2))));
        }
        else if (detectfn == 16) {  // hazard exponential
            return (exp(-r / par(1)));
        }
        else if (detectfn == 17) {  // hazard annular normal
            return (exp(-(r-par(2))*(r-par(2)) / 
                    2 / par(1)/ par(1)));
        }
        else if (detectfn == 18) {  // hazard cumulative gamma
            boost::math::gamma_distribution<> gam(par(2),par(1)/par(2));
            return (boost::math::cdf(complement(gam,r))); 
        }
        else if (detectfn == 19) {  // hazard variable power
            return (exp(- pow(r /par(1), par(2))));
        }
        else 
            return (R_NaN);  //Rcpp::stop("unknown or invalid detection function in gxy"));
    }
}

//--------------------------------------------------------------------------

double randomtime (double p)
    // return random event time for event with probability p 
{
    double minprob = 1e-5;
    double lambda;
    double random_U;
    
    if (p < minprob)
        return(huge);                        // ignore trivial p/lambda 
    else if (p >= 1.0)
        return (-unif_rand());                  // trick to spread P=1 
    else {
        lambda   = -log(1-p);                // rate parameter 
        random_U = unif_rand();
        if (random_U <= 0)                   // trap for zero 
            return(huge);
        else
            return (-log(random_U)/lambda);   // random exponential e.g. Ripley 1987 Algorithm 3.2 p 55 
    }
}
//----------------------------------------------------------------

void probsort (
        const int n, 
        std::vector<trap_animal> &tran)
    // Sort using Shell algorithm see Press et al 1989 p 257
    // tran is an array of trap_animal records
    
{
    double aln2i = 1.442695022;
    double tiny  = 1.0e-5;
    int nn,m,lognb2,l,k,j,i;
    trap_animal t;
    lognb2 = trunc(log(n)*aln2i+tiny);
    m = n;
    for (nn=1; nn<=lognb2; nn++)
    {
        m = m / 2;
        k = n-m;
        for (j=1; j<=k; j++)
        {
            i = j;
            lab1:    l = i+m;
            if (tran[l-1].time < tran[i-1].time)
            {
                t = tran[i-1];
                tran[i-1] = tran[l-1];
                tran[l-1] = t;
                i = i-m;
                if (i >= 1)  goto lab1;
            }
        }
    }
}    // end of probsort 
//----------------------------------------------------------------

// random count from different distributions 

double rcount (int binomN, double lambda, const double Tsk) {
    
    // Poisson 
    if (binomN == 0)
        return (R::rpois(lambda * Tsk) );
    
    // negative binomial 
    else if (binomN < 0) {
        // must use 'size, prob' parameters 
        // prob = size / (size + mu) 
        binomN = abs(binomN);
        return (R::rnbinom(binomN, binomN / (binomN+ (lambda * Tsk))) );
    }
    
    else { 
        if (fabs(Tsk-1) > 1e-10)             
            lambda = 1 - pow(1-lambda, Tsk);   // 2012-12-18 
        
        // Bernoulli 
        if (binomN == 1) {
            if (unif_rand() < lambda)
                return (1);
            else
                return (0);
        }
        // binomial 
        else
            return (R::rbinom(binomN, lambda) );
    }
}
//----------------------------------------------------------------

// return probability g(r) for given detection function fn 
// used in simsecr.cpp and trapping.cpp 
// double gr (
//         const int fn,
//         const Rcpp::NumericVector gsb,
//         const rpoint xy,
//         const rpoint animal) {
//     double r;
//     fnptrC fnp;
//     fnp = getgfns(fn);
//     r = distance1 (xy, animal);
//     return (fnp(as<std::vector<double>>(gsb),r));
// }
//----------------------------------------------------------------


double hazard (double pp) {
    if (pp > (1-fuzz))  // pp close to 1.0 - approx limit 
        pp = huge;      // g0 very large (effecti inf hazard) 
    else {
        if (pp <= 0) 
            pp = 0;
        else 
            pp = -log(1-pp);
    }
    return(pp);
}
//=============================================================

// detect may take values -
// 0  multi-catch traps
// 1  binary proximity detectors
// 2  count  proximity detectors
// 8  capped

//--------------------------------------------------------------------------

bool anycapped  (const IntegerVector detect) {
    bool capped = false;
    for (int s=0; s<detect.size(); s++) {
        if (detect[s]==8)
            capped = true;
    }
    return capped;
}

//  check if we need to consider variation among individuals 
// i.e. check if detection parameters constant for given s,k 
bool anyvarying (
        const int    nc,     // number of capture histories (or groups if PIA0 has that dim) 
        const int    ss,     // number of occasions 
        const int    nk,     // number of traps 
        const int    nmix,   // number of mixture classes 
        const IntegerVector &PIA0  // lookup which g0/sigma/b combination to use for given n, S, K [naive] 
) {
    int i,n,s,k,x;
    int wxi;
    bool indiv = false;
    for (s=0; s<ss; s++) {
        for (k=0; k<nk; k++) {
            for (x=0; x<nmix; x++) {
                wxi = i4(0,s,k,x,nc,ss,nk);       
                i = PIA0[wxi];
                for (n=1; n<nc; n++) {
                    wxi = i4(n,s,k,x,nc,ss,nk);    
                    if (i != PIA0[wxi]) {
                        indiv = true; break;
                    }
                }
            }
        }
    }
    return(indiv);
}
//--------------------------------------------------------------------

bool allcapped  (const IntegerVector detect) {
    bool OK = true;
    for (int s=0; s<detect.size(); s++) {
        OK = OK && (detect[s] == 8);
    }
    return OK;
}

bool allmulti (const IntegerVector detect) {
    bool notmulti = false;
    for (int s=0; s<detect.size(); s++) {
        if (detect[s]!=0)
            notmulti = true;
    }
    return (!notmulti);
}
//--------------------------------------------------------------------

// Do parameter values for naive animals differ at all from those for other animals ?
bool anyb (const NumericMatrix &gsbval, const NumericMatrix &gsb0val) {
    bool identical = true;
    for (int i=0; i<gsbval.size(); i++) {
        if (gsbval[i] != gsb0val[i]) identical = false;
    }
    return (!identical);
}
//==============================================================================
// 

void fillngcpp(
        const int nc, 
        const int gg, 
        const IntegerVector &grp, 
        std::vector<int> &ng) {
    int g, n;
    // Count number per group (not used for CL)                
    // Assume histories sorted by group = individual           
    // CH are numbered 0 <= n < nc in C code                  
    for (g=0; g<gg; g++)
        ng[g] = 0;
    for (n=0; n<nc; n++) { 
        g = grp[n] - 1; 
        ng[g]++;
    }
}

//=============================================================


// probability of count for session s, detector k, animal i
// The argument 'g' is understood to be a cumulative hazard if binomN=0,
// a probability otherwise

double pski ( int binomN,
              int count,
              double Tski,
              double g,
              double pI) {
    
    double lambda;
    double result = 1.0;
    
    if (binomN == -1) {                              // binary proximity detectors : Bernoulli
        if (abs(Tski-1) > 1e-10) {                   // effort not unity; adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        if (count>0)                                 
            result = g*pI;  
        else 
            result = 1 - g*pI;
    }
    else if (binomN == 0) {                          // count detectors : Poisson 
        lambda = Tski * g * pI;
        if ((count < 0) || (count>0 && lambda<=0)) {         
            result = 0;
        }
        else if (count == 0) {
            result = exp(-lambda);            // routinely apply Tsk adjustment to cum. hazard 
        }
        else {
            //result = R::dpois(count, Tski * g * pI); 
            boost::math::poisson_distribution<> pois(lambda);
            result = boost::math::pdf(pois,count);
        }
    }
    else if (binomN == 1) {                          // count detectors : Binomial, size from Tsk
        result = gbinom (count, round(Tski), g*pI); 
    }
    else if (binomN > 1) {                           // count detectors : Binomial, specified size 
        if (abs(Tski-1) > 1e-10) {                   // effort not unity, adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        result = gbinom (count, binomN, g*pI);
    }
    else result = NAN; // Rcpp::stop("binomN < -1 not allowed");  // code multi -2 separately
    
    return (result);
}
//--------------------------------------------------------------------------

int par3 (int fn) {
    if ((fn==1) || (fn==3) || (fn == 5)  || (fn == 6)  || (fn == 7) || 
        (fn == 8) || (fn==10) || (fn == 11)  || (fn == 12)  || (fn == 13) || 
        (fn == 15) || (fn==17) || (fn == 18))
        return(1);
    else
        return(0);
}
//--------------------------------------------------------------------------

// https://stackoverflow.com/questions/24618370/using-rmultinom-with-rcpp
Rcpp::IntegerVector oneMultinomCall(Rcpp::NumericVector probs, int N) {
    int k = probs.size();
    Rcpp::IntegerVector ans(k);
    rmultinom(N, probs.begin(), k, ans.begin());
    return(ans);
}

//===============================================================================


