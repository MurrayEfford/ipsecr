// Simulate capture histories from model
// ******************* indicates weakness


#include <Rcpp.h>
#include "ipsecr.h"
using namespace Rcpp;

//==============================================================================

// [[Rcpp::export]]
List simdetectpointcpp (
        const int           detect,      // detector -1 single, 0 multi, 1 proximity, 2 count,... 
        const int           N, 
        const int           cc,
        const NumericVector &gk0, 
        const NumericVector &gk, 
        const NumericVector &hk0, 
        const NumericVector &hk, 
        const IntegerVector &PIA0,       // lookup which g0/sigma/b combination to use for given g, S, K [naive animal] 
        const IntegerVector &PIA1,       // lookup which g0/sigma/b combination to use for given n, S, K  [caught before] 
        const NumericMatrix &Tsk,        // ss x kk array of 0/1 usage codes or effort 
        const int           btype,       // code for behavioural response  0 none etc. 
        const int           Markov,      // learned vs transient behavioural response 0 learned 1 Markov 
        const IntegerVector &binomN      // number of trials for 'count' detector modelled with binomial 
)
{
    //  detect may take values -
    // -1  single-catch traps
    //  0  multi-catch traps
    //  1  binary proximity detectors
    //  2  count  proximity detectors
    
    int    kk = Tsk.nrow();            // number of detectors 
    int    ss = Tsk.ncol();            // number of occasions
    
    double p;
    int    i,k,s;
    int    ik;
    int    nc = 0;
    int    count = 0;
    double runif;
    int    wxi = 0;
    int    c = 0;
    double Tski = 1.0;  
    bool before;
    
    std::vector<int> caughtbefore(N * kk, 0);
    std::vector<int> x(N, 0);          // mixture class of animal i 
    
    // return values
    IntegerVector caught(N);           // caught in session 
    IntegerVector value (N*ss*kk);     // return value array
    
    //========================================================
    // 'single-catch only' declarations 
    int    nanimals;
    int    ntraps;
    int    tr_an_indx = 0;
    int    anum = 0;
    int    tnum = 0;
    int    nextcombo;
    int    finished;
    int    OK;
    double event_time;
    std::vector<int> occupied(kk);
    std::vector<double> intrap(N);
    std::vector<trap_animal> tran(N * kk);
    
    //========================================================
    // 'multi-catch only' declarations 
    std::vector<double> h(N * kk);        // multi-catch only 
    std::vector<double> hsum(N);          // multi-catch only 
    std::vector<double> cump(kk+1,0);     // multi-catch only 
    
    //========================================================
    // MAIN LINE 
    
    List nullresult = List::create(Named("n") = 0,
                                   Named("caught") = caught,
                                   Named("value") = value,
                                   Named("resultcode") = 2);
    
    RNGScope scope;             // Rcpp initialise and finalise random seed 
    
    if ((detect < -1) || (detect > 2)) {
        return(nullresult);
    }
    
    // ------------------------------------------------------------------------- 
    // MAIN LOOP 
    
    for (s=0; s<ss; s++) {
        
        // --------------------------------------------------------------------- 
        // single-catch traps 
        if (detect == -1) {
            // initialise day 
            tr_an_indx = 0;
            nanimals = N;
            ntraps   = kk;
            for (i=0; i<N; i++) intrap[i] = 0;
            for (k=0; k<kk; k++) occupied[k] = 0;
            nextcombo = 0;
            
            // make tran 
            for (i=0; i<N; i++) {  // animals 
                for (k=0; k<kk; k++) { // traps 
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        wxi =  i4(i, s, k, x[i], N, ss, kk);
                        if (before)
                            c = PIA1[wxi] - 1;
                        else 
                            c = PIA0[wxi] - 1;
                        if (c >= 0) {    // ignore unused detectors 
                            if (before)
                                p = gk[i3(c, k, i, cc, kk)];
                            else
                                p = gk0[i3(c, k, i, cc, kk)];
                            
                            if (fabs(Tski-1) > 1e-10) 
                                p = 1 - pow(1-p, Tski);    // ************************ better use hazard?
                            event_time = randomtime(p);
                            if (event_time <= 1) {
                                tran[tr_an_indx].time   = event_time;
                                tran[tr_an_indx].animal = i;    // 0..N-1 
                                tran[tr_an_indx].trap   = k;    // 0..kk-1 
                                tr_an_indx++;
                            }
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
                    intrap[anum]   = tnum+1;         // trap = k+1 
                    nanimals--;
                    ntraps--;
                }
            }
            
            for (i=0; i<N; i++) {
                if (intrap[i]>0) {
                    if (caught[i]==0) {                    // first capture of this animal 
                        nc++;
                        caught[i] = nc;                    // nc-th animal to be captured 
                    }
                    value[i3(s, intrap[i]-1, caught[i]-1, ss, kk)] = 1;  
                }
            }
        }
        
        // -------------------------------------------------------------------------- 
        // multi-catch trap; only one site per occasion 
        else if (detect == 0) {
            for (i=0; i<N; i++) {
                hsum[i] = 0;
                for (k=0; k<kk; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        wxi =  i4(i, s, k, x[i], N, ss, kk);
                        if (before)
                            c = PIA1[wxi] - 1;
                        else 
                            c = PIA0[wxi] - 1;
                        if (c >= 0) {    // ignore unused detectors 
                            if (before)
                                h[k * N + i] = Tski * hk[i3(c, k, i, cc, kk)];
                            else
                                h[k * N + i] = Tski * hk0[i3(c, k, i, cc, kk)];
                            hsum[i] += h[k * N + i];
                        }
                    }
                }
                
                for (k=0; k<kk; k++) {
                    cump[k+1] = cump[k] + h[k * N + i]/hsum[i];
                }
                if (unif_rand() < (1-exp(-hsum[i]))) {
                    if (caught[i]==0)  {        // first capture of this animal 
                        nc++;
                        caught[i] = nc;
                    }
                    // find trap with probability proportional to p
                    // searches cumulative distribution of p  
                    runif = unif_rand();
                    k = 0;
                    // while ((runif > cump[k]) && (k<kk)) k++;   // bug fix 2019-10-04
                    while ((runif > cump[k+1]) && (k<kk)) k++;
                    value[i3(s, k, caught[i]-1, ss, kk)] = 1;  
                }
            }
        }
        
        // -------------------------------------------------------------------------------- 
        // the 'proximity' group of detectors 1:2 - proximity, count 
        else if ((detect >= 1) && (detect <= 2)) {
            for (i=0; i<N; i++) {
                for (k=0; k<kk; k++) {
                    Tski = Tsk(k,s);
                    if (fabs(Tski) > 1e-10) {
                        before = bswitch (btype, N, i, k, caughtbefore);
                        wxi =  i4(i, s, k, x[i], N, ss, kk);
                        if (before)
                            c = PIA1[wxi] - 1;
                        else 
                            c = PIA0[wxi] - 1;
                        if (c >= 0) {    // ignore unused detectors 
                            if (before)
                                p = gk[i3(c, k, i, cc, kk)];
                            else
                                p = gk0[i3(c, k, i, cc, kk)];
                            if (p < -0.1) { 
                                return(nullresult);
                            }  
                            if (p>0) {
                                if (detect == 1) {
                                    if (fabs(Tski-1) > 1e-10)
                                        p = 1 - pow(1-p, Tski);
                                    count = unif_rand() < p;           // binary proximity 
                                }
                                else if (detect == 2) {             // count proximity 
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
                                    value[i3(s, k, caught[i]-1, ss, kk)] = count;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if ((btype > 0) && (s < (ss-1))) {
            // update record of 'previous-capture' status 
            if (btype == 1) {
                for (i=0; i<N; i++) {
                    if (Markov) 
                        caughtbefore[i] = 0;
                    for (k=0; k<kk; k++)
                        caughtbefore[i] = R::imax2 (value[i3(s, k, i, ss, kk)], caughtbefore[i]);
                }
            }
            else if (btype == 2) {
                for (i=0; i<N; i++) {
                    for (k=0; k<kk; k++) {
                        ik = k * (N-1) + i;
                        if (Markov) 
                            caughtbefore[ik] = value[i3(s, k, i, ss, kk)];
                        else 
                            caughtbefore[ik] = R::imax2 (value[i3(s, k, i, ss, kk)], 
                                                         caughtbefore[ik]);
                    }
                }
            }
            else {
                for (k=0;k<kk;k++) {
                    if (Markov) 
                        caughtbefore[k] = 0;
                    for (i=0; i<N; i++) 
                        caughtbefore[k] = R::imax2 (value[i3(s, k, i, ss, kk)], caughtbefore[k]);
                }
            }
        }
        
    }   // loop over s 
    
    return (List::create(Named("n") = nc, 
                         Named("caught") = caught,
                         Named("value") = value,
                         Named("resultcode") = 0));
    
}
//==============================================================================

