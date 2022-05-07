#include "ipsecr.h"

//===============================================================================

// hazard (fn 14:19) or g (fn 0:11)
double zcpp (
        const double r2,
        const int detectfn,
        const Rcpp::NumericVector gsbval)
{
    double temp;
    if ((detectfn == 0) || (detectfn == 14)) {        // halfnormal or hazard halfnormal
        return (gsbval(0) * std::exp(-r2 / 2 / gsbval(1)/ gsbval(1)));
    }
    else if (detectfn == 3) {                         // compound halfnormal
        temp = gsbval(0) * std::exp(- r2  / 2 / gsbval(1) / gsbval(1));
        if (round(gsbval(2)) > 1) temp = 1 - pow(1 - temp, gsbval(2));
        return (temp);
    }
    else {
        double r = std::sqrt(r2);
        if ((detectfn == 1) || (detectfn == 15)) {       // hazard rate or hazard hazard rate
            return (gsbval(0) * ( 1 - std::exp(- pow(r /gsbval(1), - gsbval(2)))));
        }
        else if ((detectfn == 2) || (detectfn == 16)) {  // exponential or hazard exponential
            return (gsbval(0) * std::exp(-r / gsbval(1)));
        }
        else if (detectfn == 4) {                        // uniform
            if (r<gsbval(1)) return (gsbval(0));
            else return (0);
        }
        else if (detectfn == 5) {                        // w exponential
            if (r<gsbval(2)) return (gsbval(0));
            else return (gsbval(0) * std::exp(-(r-gsbval(2)) / gsbval(1)));
        }
        else if ((detectfn == 6) || (detectfn == 17)) {  // annular normal or hazard annular normal
            return (gsbval(0) * std::exp(-(r-gsbval(2))*(r-gsbval(2)) / 2 /
                gsbval(1) / gsbval(1)));
        }
        else if (detectfn == 7) {                        // cumulative lognormal
            double CV2, meanlog, sdlog;
            CV2 = gsbval(2)*gsbval(2)/gsbval(1)/gsbval(1);
            meanlog = log(gsbval(1)) - log(1 + CV2)/2;
            sdlog = std::sqrt(log(1 + CV2));
            boost::math::lognormal_distribution<> ln(meanlog,sdlog);
            return (gsbval(0) * boost::math::cdf(complement(ln,r)));
        }
        else if ((detectfn == 8) || (detectfn == 18)) {  // cumulative gamma or hazard cumulative gamma
            boost::math::gamma_distribution<> gam(gsbval(2), gsbval(1)/gsbval(2));
            return (gsbval(0) * boost::math::cdf(complement(gam,r)));
        }
        else if (detectfn == 19) {  // hazard variable power
            return (gsbval(0) * std::exp(- pow(r /gsbval(1), gsbval(2))));
        }
        else (Rcpp::stop("unknown or invalid detection function"));
    }
}

// generate a population of animals distributed according to cell density

//==============================================================================

// btype is code for behavioural response 
// 0 none
// 1 individual
// 2 individual, trap-specific
// 3 trap-specific
// NOTE: behavioural responses not checked

int bswitch (
        const int btype, 
        const int N, 
        const int i, 
        const int k, 
        const std::vector<int> &caughtbefore)
{
    if (btype == 0)
        return(0);
    else if (btype == 1) 
        return(caughtbefore[i]);
    else if (btype == 2) 
        return(caughtbefore[k * (N-1) + i]);
    else if (btype == 3) 
        return(caughtbefore[k]);
    else 
        Rcpp::stop("unrecognised btype in bswitch");
    return(0);
}
//==============================================================================

// [[Rcpp::export]]
Rcpp::List CHcpp (
        const Rcpp::NumericMatrix &animals, // x-y coord
        const Rcpp::NumericMatrix &traps,   // x-y coord
        const Rcpp::NumericMatrix &Tsk,     // usage
        int   detectfn,
        int   detect,
        const Rcpp::NumericVector &gsb,
        const double        lambdak,
        const int           btype,    // code for behavioural response  0 none etc. 
        const int           Markov,   // learned vs transient behavioural response 0 learned 1 Markov 
        const Rcpp::IntegerVector &binomN   // number of trials for 'count' detector modelled with binomial 
) {
    
    //  detect may take values -
    // -1  single-catch traps
    //  0  multi-catch traps
    //  1  binary proximity detectors
    //  2  count  proximity detectors
    
    int N = animals.nrow();
    // not to be called with N < 1
    if (N<1) Rcpp::stop ("no animals in population");
    int N1 = N + (lambdak>0);   // increment if nontarget to be modelled
    int K = Tsk.nrow();
    int S = Tsk.ncol();
    int i,k,n,s;
    double d2;
    double h0;   // intermediate value of hazard
    Rcpp::NumericMatrix hik (N1,K);
    
    double p;
    int    ik;
    int    nc = 0;
    int    count = 0;
    double runif;
    double Tski = 1.0;  
    bool   before;
    
    std::vector<int> caughtbefore(N * K, 0);
    
    // return values
    Rcpp::IntegerVector caught(N1);         // caught in session 
    Rcpp::IntegerVector value (N1*S*K);     // return value array
    Rcpp::IntegerMatrix nontarget (K, S);   // return value array
    
    for (n = 0; n<N; n++) {
        for (k=0; k<K; k++) {
            d2 = d2cpp (n, k, animals, traps);
            hik(n,k) = zcpp(d2, detectfn, gsb);
            if (detectfn<13) {
                hik(n,k) = -std::log(1-hik(n,k)); 
            }
            
        }
    }
    // constant hazard for nontarget process
    if (N1>N) {
        for (k=0; k<K; k++) {
            hik(N1-1,k) = lambdak;
        }
    }
    
    //========================================================
    // 'single-catch only' declarations 
    int    nanimals;   // dynamic
    int    ntraps;     // dynamic
    int    tr_an_indx = 0;
    int    anum = 0;
    int    tnum = 0;
    int    nextcombo;
    int    finished;
    int    OK;
    double event_time;
    std::vector<int> occupied(K);
    std::vector<double> intrap(N);
    std::vector<trap_animal> tran(N * K);
    
    //========================================================
    // 'multi-catch and capped only' declarations 
    std::vector<double> h(N1 * K);        
    std::vector<double> hsum(R::imax2(N1,K)); 
    std::vector<double> cump(K+1,0);     // multi-catch only 
    std::vector<double> cumk(N1+1,0);     // capped only 
    
    //========================================================
    // MAIN LINE 
    
    Rcpp::List nullresult = Rcpp::List::create(
        Rcpp::Named("n") = 0,
        Rcpp::Named("caught") = caught,
        Rcpp::Named("value") = value,
        Rcpp::Named("nontarget") = nontarget,
        Rcpp::Named("resultcode") = 2);
    
    Rcpp::RNGScope scope;             // Rcpp initialise and finalise random seed 
    
    if ((detect < -1) || (detect > 2 && detect != 8)) {
        return(nullresult);
    }
    
    // ------------------------------------------------------------------------- 
    // MAIN LOOP 
    
    for (s=0; s<S; s++) {
        
        // --------------------------------------------------------------------- 
        // single-catch traps 
        if (detect == -1) {
            // initialise day 
            tr_an_indx = 0;
            nanimals = N;                          // only real animals
            ntraps   = K;
            for (i=0; i<N; i++) intrap[i] = 0;
            for (k=0; k<K; k++) occupied[k] = 0;
            nextcombo = 0;
            
            // make tran 
            for (i=0; i<N1; i++) {  // animals, including nontarget
                for (k=0; k<K; k++) { // traps 
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        // learned response not implemented
                        // before = bswitch (btype, N, i, k, caughtbefore);
                        // if (before)
                        //     h = hik(i,k);
                        // else
                        h0 = hik(i,k);
                        if (fabs(Tski-1) > 1e-10) {
                            h0 = Tski * h0;
                        }
                        event_time = randomtime(1-exp(-h0));
                        if (event_time <= 1) {
                            tran[tr_an_indx].time   = event_time;
                            tran[tr_an_indx].animal = i;    // 0..N1-1 
                            tran[tr_an_indx].trap   = k;    // 0..K-1 
                            tr_an_indx++;
                        }
                    }
                }
            }
            // end of make tran 
            
            if (tr_an_indx > 1) probsort (tr_an_indx, tran);
            
            while ((nextcombo < tr_an_indx) && (nanimals>0) && (ntraps>0)) {
                finished = 0;
                OK       = 0;
                while ((1-finished)*(1-OK) > 0) {      // until finished or OK 
                    if (nextcombo >= (tr_an_indx))
                        finished = 1;                  // no more to process 
                    else {
                        anum = tran[nextcombo].animal;
                        tnum = tran[nextcombo].trap;
                        OK = (1-occupied[tnum]) * (1-intrap[anum]); // not occupied and not intrap 
                        nextcombo++;
                    }
                }
                if (finished==0) {
                    // Record this capture 
                    occupied[tnum] = 1;
                    ntraps--;
                    if (anum<N) {
                        intrap[anum]   = tnum+1;         // trap = k+1 
                        nanimals--;
                    }
                    else {
                        nontarget(tnum,s) = 1;
                    }
                }
            }
            
            for (i=0; i<N; i++) {
                if (intrap[i]>0) {
                    if (caught[i]==0) {                    // first capture of this animal 
                        nc++;
                        caught[i] = nc;                    // nc-th animal to be captured 
                    }
                    value[i3(s, intrap[i]-1, caught[i]-1, S, K)] = 1;  
                }
            }
        }
        
        // -------------------------------------------------------------------------- 
        // multi-catch trap; only one site per occasion 
        else if (detect == 0) {
            for (i=0; i<N; i++) {
                hsum[i] = 0;
                for (k=0; k<K; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        // before = bswitch (btype, N, i, k, caughtbefore);
                        // if (before)
                        //     h[k * N + i] = Tski * -log(1-gik(i,k));
                        // else
                        h[k * N + i] = Tski * hik(i,k);
                        hsum[i] += h[k * N + i];
                    }
                }
                
                for (k=0; k<K; k++) {
                    cump[k+1] = cump[k] + h[k * N + i]/hsum[i];
                }
                if (R::runif(0,1) < (1-exp(-hsum[i]))) {
                    if (caught[i]==0)  {        // first capture of this animal 
                        nc++;
                        caught[i] = nc;
                    }
                    // find trap with probability proportional to p
                    // searches cumulative distribution of p  
                    runif = R::runif(0,1);
                    k = 0;
                    // while ((runif > cump[k]) && (k<K)) k++;   // bug fix 2019-10-04
                    while ((runif > cump[k+1]) && (k<K)) k++;
                    value[i3(s, k, caught[i]-1, S, K)] = 1;  
                }
            }
        }
        
        // -------------------------------------------------------------------------- 
        // capped proximity detectors; only one detection per site per occasion 
        else if (detect == 8) {
            for (k=0; k<K; k++) {
                Tski = Tsk(k,s);
                if (fabs(Tski) > 1e-10) {
                    hsum[k] = 0;
                    for (i=0; i<N1; i++) {
                        // currently no behavioural response 2022-04-21
                        // before = bswitch (btype, N, i, k, caughtbefore);
                        // if (before)
                        //     h[k * N + i] = Tski * hik(i,k);
                        // else
                        h[k * N1 + i] = Tski * hik(i,k);
                        hsum[k] += h[k * N1 + i];
                    }
                }
            }
            
            // work in progress
            for (k=0; k<K; k++) {
                Tski = Tsk(k,s);
                if (fabs(Tski) > 1e-10) {
                    for(i=0; i<N1; i++) {  // includes nontarget
                        cumk[i+1] = cumk[i] + h[k * N1 + i]/hsum[k];
                    }
                    if (R::runif(0,1) < (1-exp(-hsum[k]))) {
                        // find animal with probability proportional to p
                        // searches cumulative distribution of p  
                        runif = R::runif(0,1);
                        i = 0;
                        while ((runif > cumk[i+1]) && (i<N1)) i++;
                        if (i<N) {
                            // Rprintf("trapped animal i %4d \n", i);
                            if (caught[i]==0)  {        // first capture of this animal 
                                nc++;
                                caught[i] = nc;
                            }
                            value[i3(s, k, caught[i]-1, S, K)] = 1;  
                        }
                        else {
                            nontarget(k,s) = 1;
                        }
                    }
                }
            }
        }    
        // -------------------------------------------------------------------------------- 
        // the 'proximity' group of detectors 1:2 - proximity, count 
        else if ((detect >= 1) && (detect <= 2)) {
            for (i=0; i<N; i++) {
                for (k=0; k<K; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        // before = bswitch (btype, N, i, k, caughtbefore);
                        // if (before)
                        //     p = gik(i,k);
                        // else
                        h0 = hik(i,k);
                        // if (p < -0.1) { 
                        //     return(nullresult);
                        // }  
                        if (h0>0) {
                            if (detect == 1) {
                                if (fabs(Tski-1) > 1e-10) {
                                    h0 = h0 * Tski;
                                }
                                p = 1 - std::exp(-h0);
                                count = R::runif(0,1) < p;  // binary proximity 
                            }
                            else if (detect == 2) {             // count proximity 
                                p = 1 - std::exp(-h0);
                                if (binomN[s] == 1)
                                    count = rcount(round(Tski), p, 1);
                                else
                                    count = rcount(binomN[s], p, Tski);
                            }
                            if (count>0) {
                                if (caught[i]==0) {              // first capture of this animal 
                                    nc++;
                                    caught[i] = nc;
                                }
                                value[i3(s, k, caught[i]-1, S, K)] = count;
                            }
                        }
                    }
                }
            }
        }
        
        
        if ((btype > 0) && (s < (S-1))) {
            // update record of 'previous-capture' status 
            if (btype == 1) {
                for (i=0; i<N; i++) {
                    if (Markov) 
                        caughtbefore[i] = 0;
                    for (k=0; k<K; k++)
                        caughtbefore[i] = R::imax2 (value[i3(s, k, i, S, K)], caughtbefore[i]);
                }
            }
            else if (btype == 2) {
                for (i=0; i<N; i++) {
                    for (k=0; k<K; k++) {
                        ik = k * (N-1) + i;
                        if (Markov) 
                            caughtbefore[ik] = value[i3(s, k, i, S, K)];
                        else 
                            caughtbefore[ik] = R::imax2 (value[i3(s, k, i, S, K)], 
                                caughtbefore[ik]);
                    }
                }
            }
            else {
                for (k=0;k<K;k++) {
                    if (Markov) 
                        caughtbefore[k] = 0;
                    for (i=0; i<N; i++) 
                        caughtbefore[k] = R::imax2 (value[i3(s, k, i, S, K)], caughtbefore[k]);
                }
            }
        }
        
    }   // loop over s 
    
    return (Rcpp::List::create(
            Rcpp::Named("n") = nc, 
            Rcpp::Named("caught") = caught,
            Rcpp::Named("value") = value,
            Rcpp::Named("nontarget") = nontarget,
            Rcpp::Named("resultcode") = 0));
    
}


//==============================================================================