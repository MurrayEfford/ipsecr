## summary.ipsecr.R
## 2022-04-01, 2022-06-14

summary.ipsecr <- function (object, newdata = NULL, alpha = 0.05, ...) {
    
    out <- vector('list')
    
    out$versiontime <- paste0(object$version, ', run ', object$starttime, ', elapsed ', round(object$proctime,2), ' s')
    
    ###################
    ## Data description
    
    trp <- traps(object$capthist)
    out$traps <- data.frame (Detector = detector(trp)[1],
        Number = nrow(trp),
        Spacing = spacing(trp))
    if (is.null(usage(trp)))
        out$traps$UsagePct <- 100
    else
        out$traps$UsagePct <- 100 * sum(usage(trp))/length(usage(trp))
    if (length(detector(trp))>1)
        out$detector <- detector(trp)
    out$capthist <- summary(object$capthist, terse = TRUE, moves = TRUE)
    
    out$mask <- data.frame(Cells = nrow(object$mask), Spacing = spacing(object$mask))
    if (length(maskarea(object$mask))==0)
        out$mask <- cbind(out$mask, Length = masklength(object$mask))
    else
        out$mask <- cbind(out$mask, Area = maskarea(object$mask))
    
    ####################
    ## Model description
    
    out$modeldetails <- data.frame(
        detectfn = object$detectfn,
        fixed = fixed.string(object$fixed),
        distribution = object$details$distribution
    )
    
    ####################
    ## Completion code
    codetext <- paste (
        c('OK', 
            'reached maximum number of boxes, ', 
            'reached maximum number of simulations, '),
        c("",  
            object$details$max.nbox, 
            object$details$max.nsim)
    )
    out$completioncode <- paste0("Completion code ", object$code, ": ", codetext[object$code])
    
    out$coef <- coef(object)
    
    out$predicted <- predict (object, newdata, type = "response", alpha = alpha, 
        se.fit = object$details$var.nsim>1)
    
    ## remove distracting row names
    for (i in 1:length(out)) {
        if (is.data.frame(out[[i]]))
            if (nrow(out[[i]])==1 & (!names(out)[i] %in% c('coef', 'derived', 'predicted')))
                rownames (out[[i]]) <- ''
    }
    class(out) <- c("summary.ipsecr")
    out
}
############################################################################################

print.summary.ipsecr <- function (x, ...) {
    class(x) <- NULL
    print(x)
}
############################################################################################

