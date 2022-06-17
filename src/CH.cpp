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

// int bswitch (
//         const int btype, 
//         const int N, 
//         const int i, 
//         const int k, 
//         const std::vector<int> &caughtbefore)
// {
//     if (btype == 0)
//         return(0);
//     else if (btype == 1) 
//         return(caughtbefore[i]);
//     else if (btype == 2) 
//         return(caughtbefore[k * (N-1) + i]);
//     else if (btype == 3) 
//         return(caughtbefore[k]);
//     else 
//         Rcpp::stop("unrecognised btype in bswitch");
//     return(0);
// }
//==============================================================================

// [[Rcpp::export]]
Rcpp::List CHcpp (
        const Rcpp::NumericMatrix &animals, // x-y coord
        const Rcpp::NumericMatrix &traps,   // x-y coord
        const Rcpp::NumericMatrix &Tsk,     // usage
        const Rcpp::NumericVector &gsb,
        const Rcpp::NumericVector &NT,
        const int                 detectfn,
        const int                 detectorcode,
        const int                 nontargetcode,
        const int                 btype,    // code for behavioural response  0 none etc. 
        const int                 Markov,   // learned vs transient behavioural response 0 learned 1 Markov 
        const Rcpp::IntegerVector &binomN   // number of trials for 'count' detector modelled with binomial 
) {
    
    //  detectorcode may take values -
    // -1  single-catch traps
    //  0  multi-catch traps
    //  1  binary proximity detectors
    //  2  count  proximity detectors
    //  8  capped binary proximity detectors
    
    int N = animals.nrow();
    // not to be called with N < 1
    if (N<1) Rcpp::stop ("no animals in population");
    int N1 = N + (nontargetcode==1);   // increment for exclusive nontargets
    int K = Tsk.nrow();
    int S = Tsk.ncol();
    int i,k,j,n,s;
    int isk;
    double d2;
    double h0;   // intermediate value of hazard
    Rcpp::NumericMatrix hik (N1,K);
    
    double p;
    int    nc = 0;
    int    count = 0;
    double Tski = 1.0;  
    // int    ik;
    // bool   before;
    
    std::vector<int> caughtbefore(N * K, 0);
    
    // return values
    Rcpp::IntegerVector CH (N*S*K);        // return value array
    Rcpp::IntegerMatrix nontarget (K, S);  // return nontarget array
    
    //========================================================
    for (n = 0; n<N; n++) {
        for (k=0; k<K; k++) {
            d2 = d2cpp (n, k, animals, traps);
            hik(n,k) = zcpp(d2, detectfn, gsb);
            if (detectfn<13) {
                hik(n,k) = -std::log(1-hik(n,k)); 
            }
            
        }
    }
    // notional animal N1 for exclusive nontarget process
    // used only for single-catch traps
    if (N1>N) {
        for (k=0; k<K; k++) {
            hik(N1-1,k) = NT[k];
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
    std::vector<int> occupied(K);        // single, multi
    std::vector<double> intrap(N);
    std::vector<trap_animal> tran(N1 * K);
    double maxt;
    
    //========================================================
    // 'multi' and 'proximity' and 'capped' declarations 
    std::vector<double> inttime(K,1.0);
    double dettime = 1.0;
    double tmptime;
    
    //========================================================
    // MAIN LINE 
    
    Rcpp::List nullresult = Rcpp::List::create(
        // Rcpp::Named("n") = 0,
        // Rcpp::Named("caught") = caught,
        Rcpp::Named("CH") = CH,
        Rcpp::Named("nontarget") = nontarget,
        Rcpp::Named("resultcode") = 2);
    
    Rcpp::RNGScope scope;             // Rcpp initialise and finalise random seed 
    
    if ((detectorcode < -1) || (detectorcode > 2 && detectorcode != 8)) {
        return(nullresult);
    }
    
    // ------------------------------------------------------------------------- 
    // MAIN LOOP 
    
    for (s=0; s<S; s++) {
        
        // --------------------------------------------------------------------- 
        // interference events, excluding dependent 
        if (nontargetcode > 0 && nontargetcode != 5) {
            for (k=0; k<K; k++) {
                if (fabs(Tsk(k,s))>1e-10 && NT[k]>0) {
                    // random time of interference
                    inttime[k] = R::rexp(1/(NT[k] * Tsk(k,s)));  // scale not rate
                }
                else {
                    inttime[k] = 1;
                }
            }
        }
        
        // --------------------------------------------------------------------- 
        // single-catch traps 
        if (detectorcode == -1) {
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
                        if (nontargetcode == 2)
                            maxt = inttime[k];
                        else
                            maxt = 1.0;
                        if (event_time <= maxt) {
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
                        nontarget(tnum,s) = 1;         // exclusive nontarget
                    }
                }
            }
            
            for (i=0; i<N; i++) {
                if (intrap[i]>0) {
                    CH[i3(i,s,intrap[i]-1, N, S)] = 1;  
                }
            }
        }
        
        // -------------------------------------------------------------------------- 
        
        // multi-catch trap; only one site per occasion 
        if (detectorcode == 0) {
            for (i=0; i<N; i++) {
                dettime = 1.0;
                for (k=0; k<K; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        h0 = hik(i,k);
                        if (h0>0) {
                            tmptime = R::rexp(1/(h0 * Tski));
                            if (tmptime < dettime && 
                                (tmptime < inttime[k] || nontargetcode > 2)) {
                                dettime = tmptime;
                                tnum = k;
                            }
                        }
                    }
                }
                if (dettime < 1.0) {
                    CH[i3(i, s, tnum, N, S)] = 1;
                }
            }
        }
        // binary capped
        else if (detectorcode == 8) {
            for (k=0; k<K; k++) {
                occupied[k] = 0;
                dettime = 1.0;
                Tski = Tsk(k,s);
                if (fabs(Tski) > 1e-10) {
                    for (i=0; i<N; i++) {
                        h0 = hik(i,k);
                        if (h0>0) {
                            tmptime = R::rexp(1 / (h0*Tski));
                            if (tmptime < dettime) {
                                anum = i;
                                dettime = tmptime;
                            }
                        }
                    }
                    if ((dettime < 1.0) && 
                        (dettime < inttime[k] || nontargetcode > 2)) {
                        occupied[k] = 1;
                        CH[i3(anum, s, k, N, S)] = 1;
                    }
                }
            }
        }
    
        // the 'proximity' group of detectors: 1 proximity, 2 count
        else if (detectorcode == 1 || detectorcode == 2) {
            for (k=0; k<K; k++) {
                occupied[k] = 0;
                Tski = Tsk(k,s);
                if (fabs(Tski) > 1e-10) {
                    for (i=0; i<N; i++) {
                        h0 = hik(i,k);
                        if (h0>0) {
                            if (detectorcode == 1) {    // binary proximity 
                                // random time of detection
                                dettime = R::rexp(1/ (h0*Tski));
                                count = (dettime < 1) && 
                                    (dettime < inttime[k] || nontargetcode > 2);
                                if (count>0) {
                                    occupied[k] = 1;
                                }
                            }
                            else if (detectorcode == 2) {             // count proximity 
                                if (binomN[s] == 0) {                 // Poisson 2022-06-15
                                   count = R::rpois(h0 * Tski);
                                }
                                else {
                                    p = 1 - std::exp(-h0);
                                    if (binomN[s] == 1)
                                        count = rcount(round(Tski), p, 1);
                                    else
                                        count = rcount(binomN[s], p, Tski);
                                }
                                
                                if (count > 0 && nontargetcode>0 && nontargetcode != 5) {
                                    for (i=0; i<count; i++) {
                                        // random time of detection
                                        // NOT QUITE RIGHT BECAUSE ALREADY TIME<1
                                        dettime = R::rexp(1/h0);
                                        // discount detections after interference 
                                        if (dettime > inttime[k]) count--;
                                    }
                                }                                
                            }
                            if (count>0) {
                                CH[i3(i, s, k, N, S)] = count;
                            }
                        }
                    }
                }
            }
        }
        
        if (nontargetcode>0) {
            
            // dependent nontarget
            // reclassify fraction NT[k] of detections to 'not identified'
            if (nontargetcode == 5 && detectorcode == 2) {
                for (k=0; k<K; k++) {
                    for (i=0; i<N; i++) {
                        isk = i3(i, s, k, N, S);
                        for(j=1; j<=CH[isk]; j++) {
                            if (R::runif(0,1) < NT[k]) {
                                // transfer from ID to nontarget
                                nontarget(k,s)++;
                                CH[isk]--;
                            }
                        }
                    }
                }
            }
            
            // non-exclusive interference
            else if (nontargetcode > 1) {
                for (k=0; k<K; k++) {
                    if (inttime[k] < 1) {
                        nontarget(k,s) = 1;
                        if (nontargetcode == 3) {
                            // any detections erased
                            for (i=0; i<nc; i++) {
                                CH[i3(i,s,k,N,S)] = 0;
                            }  
                        }
                    }
                }
            }
            
            // capped detectors, exclusive interference
            else if (nontargetcode==1 && detectorcode== 8) {
                for (k=0; k<K; k++) {
                    if (inttime[k] < 1 && !occupied[k]) {
                        nontarget(k,s) = 1;
                    }
                }
            }
        }
        
    }   // loop over s 
    
    return (Rcpp::List::create(
            // Rcpp::Named("n") = nc, 
            // Rcpp::Named("caught") = caught,
            Rcpp::Named("CH") = CH,
            Rcpp::Named("nontarget") = nontarget,
            Rcpp::Named("resultcode") = 0));
    
}


//==============================================================================
