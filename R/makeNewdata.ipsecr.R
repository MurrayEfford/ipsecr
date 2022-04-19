############################################################################################
## package 'ipsecr'
## makeNewData.ipsecr
## 2022-03-31
############################################################################################

makeNewData.ipsecr <- function (object, all.levels = FALSE, ...) {
        
    autovars <- c('x','y','x2','y2','xy',
                  't','T','tcov')
    capthist <- object$capthist
    mask <- object$mask
    vars <- object$vars
    timecov <- object$timecov

    nocc <- ncol (capthist)

    onesession <- function(session) {
        findvars <- function (basevars, cov) {
            ## function to add covariates to a list
            ## cov should be dataframe or list of dataframes, one per session (R > 1),
            if (!is.data.frame(cov)) cov <- cov[[session]] ## assume multisession list
            if (is.null(cov) | (length(cov)==0) | (length(sessvars)==0)) return(basevars)
            else {
                found <- ''
                for (v in sessvars) {
                    if (v %in% names(cov)) {
                        vals <- cov[,v]
                        if (is.character(vals)) vals <- factor(vals)
                        basevars[[v]] <- if (is.factor(vals))
                            factor(levels(vals), levels = levels(vals))
                        else
                            unique(vals)
                        
                        found <- c(found, v)
                    }
                }
                sessvars <<- sessvars[!(sessvars %in% found)]
                return(basevars)
            }
        }
        
        sessvars <- vars
        
        basevars <- list(session = 1)
        
        for (v in sessvars) {
            if (v=='x')  basevars$x <- 0     # mean attr(mask,'meanSD')[1,'x']
            if (v=='y')  basevars$y <- 0     # mean attr(mask,'meanSD')[1,'y']
            if (v=='x2') basevars$x2 <- 0   # mean attr(mask,'meanSD')[1,'x']
            if (v=='y2') basevars$y2 <- 0   # mean attr(mask,'meanSD')[1,'y']
            if (v=='xy') basevars$xy <- 0   # mean attr(mask,'meanSD')[1,'x']
            if (v=='T')  basevars$T <- 0   
            
            if (v=='t')  basevars$t <- factor(1:nocc)

            if (v=='tcov') {
                timecov <- object$timecov
                if (is.factor(timecov)) {
                    basevars$tcov <- unique(timecov)
                }
                else
                    basevars$tcov <- 0        # ideally use mean or standardize?
            }
        }
        ## all autovars should now have been dealt with
        basevars <- findvars (basevars, timecov)
        basevars <- findvars (basevars, covariates(traps(capthist)))
        basevars <- findvars (basevars, covariates(mask))
        
        ## revert to first level (original default)
        for (v in names(basevars)) {
            if (!all.levels) {
                basevars[[v]] <- basevars[[v]][1] 
            }
        }
        expand.grid(basevars)
    }
    
    newdata <- onesession(1)
    newdata
    
}
############################################################################################

