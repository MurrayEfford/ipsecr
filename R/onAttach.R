###############################################################################
## package 'ipsecr'
## onAttach.R
## last changed 2022-04-01
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- paste0(packageVersion('ipsecr'), .localstuff$packageType)
    packageStartupMessage( "This is ipsecr ", version,
                           ". For overview type ?ipsecr" )
}

## .onLoad is preferred if actions are required for single functions 
## that may be called without attaching package