#include "ipsecr.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::ucube armaCHcpp( 
        const arma::mat &d,         // precomputed animal-detector distances 
        const arma::mat &Tsk,       // detector usage 
        const arma::mat &detpar,    // detection parameters N x np
        const arma::vec &NT,        // interference hazard K
        const arma::ivec &binomN,   // count distribution S (0 Poisson etc.)
        int detectfn,               // integer code for detection function
        int detectorcode,           // integer code for detector type
        int nontargetcode) {        // integer code for non-target interference type
    
    if (!(detectorcode <= 2 || detectorcode == 8)) {    
        Rcpp::stop ("detectorcode not available");
    }
    if (nontargetcode == 3 || nontargetcode == 8) {
        Rcpp::stop ("nontargetcode not implemented");
    }
        
    arma::uword N = d.n_rows;
    arma::uword K = d.n_cols;
    arma::uword S = Tsk.n_cols;
    arma::uword N1 = N + (nontargetcode>0);  // increment for interference process
    
    arma::mat hik (N,K); 
    arma::mat hiks;
    arma::uword i,j,k,s,ik;
    arma::uvec ids;
    arma::umat animal;
    arma::umat trap;
    arma::uvec indices;
    arma::umat animalv;
    arma::umat trapv;
    arma::umat sv;
    arma::umat indmat;
    arma::mat capture_time;
    
    //-------------------------------------------------------------------------
    // each N,K parameter matrix 
    // varying by individual, constant across detectors and times
    arma::mat par0 = repmat(detpar.col(0), 1, K);
    arma::mat par1 = repmat(detpar.col(1), 1, K);
    arma::mat par2;
    if (detpar.n_cols == 3) {
        par2 = repmat(detpar.col(2), 1, K);
    }
    
    if (detectfn == 0 || detectfn == 14) {  // HHN
        hik = par0 % exp(-square(d) / 2 / square(par1));
    }
    else if (detectfn == 1 || detectfn == 15) {  // HHR
        for (i=0; i<N; i++) {
            for (k=0; k<K; k++) {
                hik(i,k) = par0(i,k) * (1 - exp(- pow(d(i,k) /par1(i,k), -par2(i,k))));
            }
        }
    } 
    else if (detectfn == 2 || detectfn == 16) {  // HEX
        hik = par0 % exp(-d / par1);
    } 
    else if (detectfn == 6 || detectfn == 17) {  // HAN
        hik = par0 % exp(-square(d-par2) / 2 / square(par1));
    } 
    else if (detectfn == 8 || detectfn == 18) {  // HCG
        for (i=0; i<N; i++) {
            for (k=0; k<K; k++) {
                boost::math::gamma_distribution<> gam(par2(i,k), par1(i,k)/par2(i,k));
                hik(i,k) = par0(i,k) * boost::math::cdf(complement(gam,d(i,k)));
            }
        }
    } 
    else if (detectfn == 19) {  // HVP
        for (i=0; i<N; i++) {
            for (k=0; k<K; k++) {
                hik(i,k) = par0(i,k) * exp(- pow(d(i,k) /par1(i,k), par2(i,k)));
            }
        }
    }
    else Rcpp::stop ("detectfn not implemented");
    
    bool binomial = (detectorcode == 2) && (binomN(0) > 0);
    if (detectfn < 14 && !binomial) {
        // convert detectfn p to hazard; may be zero
        hik = -log(1-hik);
    }
    if (detectfn >= 14 && binomial) {
        // convert hazard to detectfn p t
        hik = 1 - exp(-hik);
    }
    //-------------------------------------------------------------------------
    
    // optional non-target animal with trap-specific lambda in NT
    if (nontargetcode > 0) {
        hik = arma::join_vert(hik, NT.t());
    }
    
    // output
    arma::ucube CH (N1,S,K, arma::fill::zeros);
    
    // general work
    arma::mat event_time (arma::size(hik));
    
    arma::umat animal0 (arma::size(hik));
    arma::umat trap0 (arma::size(hik));
    arma::uvec an_avail (hik.n_rows);   // N or N+1
    arma::uvec tr_avail (K);
    
    for(i=0; i<animal0.n_rows; i++) animal0.row(i).fill(i);
    for(k=0; k<K; k++) trap0.col(k).fill(k);
    
    for (s=0; s<S; s++) {
        // counts (2)
        if (detectorcode == 2) {       
            // independent interference not allowed
            for(i=0; i<N; i++) {
                for(k=0; k<K; k++) {
                    if (binomN(s) == 0)
                        CH(i,s,k) = R::rpois(Tsk(k,s) * hik(i,k));   
                    else if (binomN(s) > 1)
                        CH(i,s,k) = R::rbinom(Tsk(k,s), hik(i,k));   
                    else
                        CH(i,s,k) = R::rbinom(binomN(s), hik(i,k));  
                    // dependent interference
                    if (nontargetcode == 5) {   
                        for(j = 1; j <= CH(i, s, k); j++) {
                            if (R::runif(0,1) < NT(k)) {
                                // transfer from ID to nontarget
                                CH(i, s, k)--;
                                CH(N, s, k)++;
                            }
                        }
                    }
                }
            }
        }
        // single (-1), multi (0), proximity (1), capped (8) 
        else {
            animal = animal0;
            trap = trap0;

            // scale hazard by effort (Tsk)
            hiks = hik.each_row() % Tsk.col(s).t();
            
            // protect against 1/0
            hiks.replace(0, arma::datum::eps);
            
            // random exponential latent capture times
            event_time.randu();
            event_time = -log(event_time) / hiks;
            
            // alternative random exponential (safer)
            // for (i=0; i<N1; i++)
            //     for (k=0; k<K; k++)
            //         event_time(i,k) = R::rexp(1 / hiks(i,k));
            
            // select latent captures in unit time
            ids = find(event_time<1);             
            capture_time = event_time.elem(ids);
            animal = animal.elem(ids);
            trap = trap.elem(ids);
            
            // sort by latent capture time
            indices = arma::sort_index(capture_time);
            capture_time = capture_time.elem(indices);
            animal = animal.elem(indices);
            trap = trap.elem(indices);
            
            // all animals and detectors initially available
            an_avail.ones();
            tr_avail.ones();
            
            // loop over sorted latent capture events
            for(ik = 0; ik < capture_time.n_rows; ik++) {
                i = animal(ik);
                k = trap(ik);
                
                // no capture when animal or trap not available
                if (an_avail(i) == 0 || tr_avail(k) == 0) {
                    capture_time(ik) = -1;   
                }
                else {
                    // single (-1) and multi (0) traps remove animals
                    if (i < N && detectorcode < 1) {
                        an_avail(i) = 0;
                    }
                    
                    // single (-1), capped (8) and excl. interference remove traps
                    if (detectorcode == -1 || detectorcode == 8  ||             
                        (i == N && nontargetcode == 1)) {   
                        tr_avail(k) = 0;   
                    }
                }
            }
            
            // store captures
            ids = find(capture_time>0);
            animalv = animal.elem(ids);   // column vector
            trapv = trap.elem(ids);
            sv = arma::umat (arma::size(trapv), arma::fill::value(s));    
            indmat = arma::join_horiz(animalv, sv, trapv ).t();
            indices = arma::sub2ind(arma::size(CH), indmat);
            CH.elem(indices).ones();
        }
    }
    
    return(CH);
    
}
