region.N.ipsecr <- function (object, region = NULL, spacing = NULL, session = NULL,
    se.N = TRUE, alpha = 0.05, loginterval = TRUE,
    keep.region = FALSE, nlowerbound = TRUE, RN.method = 'poisson',
    pooled.RN = FALSE, ncores = NULL, ...) {
    
    object$fit <- list(par = object$beta)
    object$CL <- FALSE
    class(object) <- 'secr'
    region.N(object, region, spacing, session,
        group = NULL, se.N, alpha, loginterval,
        keep.region, nlowerbound, RN.method,
        pooled.RN, ncores, ...) 
}
    