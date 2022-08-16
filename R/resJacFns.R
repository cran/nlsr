# resJacFns -- functions to provide residuals and Jacobians
#   including approximations
# ?? need to define output of these!!
#
#'  jand 
#' 
#' approximate Jacobian via numDeriv::jacobian
#' 
#' @param pars a named numeric vector of parameters to the model
#' @param resfn a function to compute a vector of residuals
#' @param bdmsk ?? do we need it? Default is \code{NULL}
#' @param resbest If supplied, a vector of the residuals at the parameters
#'       \code{pars} to save re-evaluation.
#' @param ndstep A tolerance used to alter parameters to compute numerical
#'       approximations to derivatives. Default \code{1e-7}.
#' @param \dots Extra information needed to compute the residuals
#' 
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
jand <- function(pars, resfn=NULL, bdmsk=NULL, resbest=NULL, ndstep=1e-7, ...){
   # ndstep not currently used
   if (is.null(resfn)) stop("You must supply a residual function")
   # numDeriv jacobian approximation
   jacmat<-numDeriv::jacobian(resfn, pars,...) # don't specify other numDeriv parameters here yet ??
   attr(jacmat,"gradient") <- jacmat # to satisfy needs of nlxb(), nlfb() working
} # end jand

#'  jafwd
#' 
#' approximate Jacobian via forward differences
#' 
#' @param pars a named numeric vector of parameters to the model
#' @param resfn a function to compute a vector of residuals
#' @param bdmsk ?? do we need it? Default is \code{NULL}
#' @param resbest If supplied, a vector of the residuals at the parameters
#'       \code{pars} to save re-evaluation.
#' @param ndstep A tolerance used to alter parameters to compute numerical
#'       approximations to derivatives. Default \code{1e-7}.
#' @param \dots Extra information needed to compute the residuals
#' 
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
jafwd <- function(pars, resfn=NULL, bdmsk=NULL, resbest=NULL, ndstep=1e-7, ...){
  # Forward difference jacobian approximation
   if (is.null(resfn)) stop("You must supply a residual function")
   npar<-length(pars)
   nres<-length(resbest)
   jacmat<-matrix(0, ncol=npar, nrow=nres)
   tpars<-pars # need to set because values used for each dimension
   for (j in 1:npar){
      if (is.null(bdmsk) || (bdmsk[j]!=0)) {
         step<-ndstep*(abs(pars[j])+ndstep)
         tpars[j]<-pars[j]+step
         jacmat[,j]<-(resfn(tpars,...)-resbest)/step
         tpars[j]<-pars[j]
      }
   }
   attr(jacmat,"gradient") <- jacmat # to satisfy needs of nlxb(), nlfb() working
   jacmat
} # end jafwd

#'  jaback
#' 
#' approximate Jacobian via forward differences
#' 
#' @param pars a named numeric vector of parameters to the model
#' @param resfn a function to compute a vector of residuals
#' @param bdmsk ?? do we need it? Default is \code{NULL}
#' @param resbest If supplied, a vector of the residuals at the parameters
#'       \code{pars} to save re-evaluation.
#' @param ndstep A tolerance used to alter parameters to compute numerical
#'       approximations to derivatives. Default \code{1e-7}.
#' @param \dots Extra information needed to compute the residuals
#' 
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
jaback <- function(pars, resfn=NULL, bdmsk=NULL, resbest=NULL, ndstep=1e-7, ...){
   # Backward difference jacobian approximation
   if (is.null(resfn)) stop("You must supply a residual function")
   npar<-length(pars)
   nres<-length(resbest)
   jacmat<-matrix(0, ncol=npar, nrow=nres)
   tpars<-pars
   for (j in 1:npar){
      if (is.null(bdmsk) || (bdmsk[j]!=0)) {
         step<-ndstep*(abs(pars[j])+ndstep)
         tpars[j]<-tpars[j]-step
         jacmat[,j]<-(resbest-resfn(tpars,...))/step
         tpars[j]<-pars[j]
      }
   }
   attr(jacmat,"gradient") <- jacmat # to satisfy needs of nlxb(), nlfb() working
   jacmat
} # end jaback

#'  jacentral
#' 
#' approximate Jacobian via central differences. Note this needs two
#'    evaluations per parameter, but generally gives much better approximation
#'    of the derivatives.
#' 
#' @param pars a named numeric vector of parameters to the model
#' @param resfn a function to compute a vector of residuals
#' @param bdmsk ?? do we need it? Default is \code{NULL}
#' @param resbest If supplied, a vector of the residuals at the parameters
#'       \code{pars} to save re-evaluation.
#' @param ndstep A tolerance used to alter parameters to compute numerical
#'       approximations to derivatives. Default \code{1e-7}.
#' @param \dots Extra information needed to compute the residuals
#' 
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
jacentral <- function(pars, resfn=NULL, bdmsk=NULL, resbest=NULL, ndstep=1e-7, ...){
   # Central difference jacobian approximation
   if (is.null(resfn)) stop("You must supply a residual function")
   npar<-length(pars)
   nres<-length(resbest)
   jacmat<-matrix(0, ncol=npar, nrow=length(resbest))
   tparf<-pars
   tparm<-pars
   for (j in 1:npar){
      if (is.null(bdmsk) || (bdmsk[j]!=0)) {
         step<-ndstep*(abs(pars[j])+ndstep)
         tparf[j]<-tparf[j]+step
         tparm[j]<-tparm[j]-step
         jacmat[,j]<-0.5*(resfn(tparf,...)-resfn(tparm,...))/step
         tparf[j]<-pars[j]
         tparm[j]<-pars[j]
      }
   }
   attr(jacmat,"gradient") <- jacmat # to satisfy needs of nlxb(), nlfb() working
   jacmat
} # end jacentral


#' resgr
#' 
#' Computes the gradient of the sum of squares function for nonlinear least
#' squares where \code{resfn} and \code{jacfn} supply the residuals and Jacobian
#' 
#' ?? does it work with approximate Jacobian functions
#' 
#' @usage resgr(prm, resfn, jacfn, ...)
#' 
#' @param prm a numeric vector of parameters to the model
#' @param resfn a function to compute a vector of residuals
#' @param jacfn  a function to compute the Jacobian of the sum of squares 
#' @param ... Extra information needed to compute the residuals
#' 
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' 
#' @export
resgr <- function(prm, resfn, jacfn, ...) {
    # computes the gradient 2 * J' %*% res for residuals (res)
    # and jacobian (Jac) defined by resfn and jacfn at
    # parameters prm and extra variables defined in the
    # dot-arguments J C Nash 2012-4-26
    res <- resfn(prm, ...)  # ?? try()
    Jac <- jacfn(prm, ...)
    grj <- 2 * as.numeric(crossprod(Jac, res))
    attr(res, "Jacobian") <- Jac 
    attr(res, "gradient") <- grj
    res
}


#'  resss
#' 
#' compute the sum of squares from \code{resfn} at parameters \code{prm}
#' 
#' @usage resss(prm, resfn, ...)
#' 
#' @param prm a named numeric vector of parameters to the model
#' @param resfn a function to compute a vector of residuals
#' @param ... Extra information needed to compute the residuals
#' 
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
resss <- function(prm, resfn, ...) {
    # computes sumsquares function from residuals defined by
    # resfn at parameters prm and extra variables defined in
    # the dot-arguments J C Nash 2012-4-26
    resids <- resfn(prm, ...)  # ?? try()
    ss <- as.numeric(crossprod(resids))
}

