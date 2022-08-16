#' nlsr-package
#' @aliases nlsr
#' 
#' Tools for solving nonlinear least squares problems
#'  
#' The package provides some tools related to using the Nash variant
#' of Marquardt's algorithm for nonlinear least squares. Jacobians can
#' usually be developed by automatic or symbolic derivatives.
#' 
#' @details 
#' This package includes methods for solving nonlinear least squares problems
#' specified by a modeling expression and given a starting vector of named 
#' paramters. Note: You must provide an expression of the form
#'     lhs ~ rhsexpression
#' so that the residual expression 
#'     rhsexpression - lhs
#' can be computed. The expression can be enclosed in quotes, and this seems to give
#' fewer difficulties with R. Data variables must already be defined, either within 
#' the parent environment or else in the dot-arguments. Other symbolic elements in
#' the modeling expression must be standard functions or else parameters that are 
#' named in the start vector.
#' 
#' The main functions in \code{nlsr} are:
#' 
#' nlfb   Nash variant of the Marquardt procedure for nonlinear least squares,
#'    with bounds constraints, using a residual and optionally Jacobian
#'    described as \code{R} functions.
#'    
#' nlxb  Nash variant of the Marquardt procedure for nonlinear least squares,
#'   	with bounds constraints, using an expression to describe the residual via
#'   	an \code{R} modeling expression. The Jacobian is computed via symbolic
#'   	differentiation.
#'   	
#' wrapnlsr  Uses \code{nlxb} to solve nonlinear least squares then calls 
#'    \code{nls()} to create an object of type nls.
#'    
#' model2rjfun  returns a function with header \code{function(prm)}, which
#'    evaluates the residuals (and if jacobian is TRUE the Jacobian matrix) 
#'    of the model at \code{prm}.  The residuals are defined to be the 
#'    right hand side of \code{modelformula} minus the left hand side.
#'    
#' model2ssgrfun  returns a function with header \code{function(prm)}, which
#'     evaluates the sum of squared residuals (and if gradient is \code{TRUE} the
#'     gradient vector) of the model at \code{prm}.
#'     
#' modelexpr  returns the expression used to calculate the vector of 
#'         residuals (and possibly the Jacobian) used in the previous functions.
#'         
#' @author John C Nash and Duncan Murdoch
#' 
#' @references 
#' 
#' Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
#'     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
#'     of Physics Publications
#'
#'   Nash, J. C. (2014) _Nonlinear Parameter Optimization Using R Tools._
#'     Wiley 
#'  
#' @keywords nls 
#' @keywords nonlinear&least&squares
#' 
nlsr <- nlsr.package <- function(){
   cat("nlsr package: a set of tools for solving nonlinear least squares problems\n")
   return(1)
}
#' 
#' 
#' 
#' nlsr.control
#' 
#' set and provide defaults of controls for package \code{nlsr}
#' 
#' @param control A list of controls. 
#'        If missing, the defaults are provided. See below.
#'        If a named control is provided,  e.g., via a call
#'        ctrllist<- nlsr.control(japprox="jand"), then
#'        that value is substituted for the default of the control
#'        in the FULL list of controls that is returned.
#' 
#' NOTE: at 2022-6-17 there is NO CHECK FOR VALIDITY
#' 
#' The set of possible controls to set is as follows, and is returned.
#'
#' @return 
#'            femax      INTEGER limit on the number of evaluations of residual function
#'                           Default  10000.
#'                           
#'            japprox    CHARACTER name of the Jacobian approximation to use
#'                           Default NULL, since we try to use analytic gradient
#'                  
#'            jemax      INTEGER limit on the number of evaluations of the Jacobian
#'                           Default 5000
#'                           
#'            lamda      REAL initial value of the Marquardt parameter
#'                           Default 0.0001 
#'            Note: mis-spelling as in JNWMS, kept as historical serendipity.
#'            
#'            lamdec     REAL multiplier used to REDUCE \code{lambda} (0 < \code{lamdec} < \code{laminc})
#'                           Default 4, so \code{lamda <- lamda * (lamdec/laminc)}
#'                           
#'            laminc     REAL multiplier to INCREASE \code{lambda} (1 < \code{laminc}
#'                           Default 10
#'                           
#'            ndstep     REAL  stepsize for numerical Jacobian approximation
#'                           Default 1e-7
#'                           
#'            numjac     LOGICAL  If TRUE use a numerical approximation to Jacobian.
#'                           Default FALSE
#'                           
#'            offset   REAL A value used to test for numerical equality, i.e. \code{a} and
#'                   \code{b} are taken equal if \code{(a + offset) == (b + offset)}
#'                           Default 100.
#'                           
#'            phi      REAL Factor used to add unit Marquardt stabilization matrix in Nash
#'                         modification of LM method.  Default 1.  choices
#'                           
#'            psi      REAL Factor used to add scaled Marquardt stabilization matrix in Nash
#'                         modification of LM method.  Default 0. 
#'                           
#'            rofftest LOGICAL If TRUE, perform (safeguarded) relative offset convergence test
#'                           Default TRUE
#'                           
#'            smallsstest LOGICAL If TRUE tests sum of squares and terminates if very small ??describe
#'                           Default TRUE
#'                           
#'            stepredn  REAL Factor used to reduce the stepsize in a Gauss-Newton algorithm (Hartley's
#'			method). 0 means NO backtrack. Default 0. 
#'                           
#'            watch      LOGICAL to provide a pause at the end of each iteration for user to monitor
#'                       progress.   Default FALSE
#'                           
#'            prtlvl     INTEGER The higher the value, the more intermediate output is provided.
#'                       Default 0
#'                           
#'            scaleOffset   a positive constant to be added to the denominator sum-of-squares in
#'                          the relative offset convergence criteria. Default 0.
#'                           
#' 
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
nlsr.control <- function(control) {
  # ?? NEED TO PUT IN CHECKS of validity
  # Defaults
  defctrl=list(femax=10000,
              japprox=NULL,
              jemax=5000,
              lamda=0.0001,
	      laminc=10, 
	      lamdec=4, 
	      ndstep=1e-7,
	      numjac=FALSE,
	      offset=100, 
	      phi=1.0, 
	      psi=0.0,
	      rofftest = TRUE, 
	      smallsstest = TRUE, 
	      stepredn=0.0,
	      watch=FALSE,
              prtlvl=0,
              scaleOffset=1.0)
  if (missing(control) || is.null(control)){control <- defctrl}
  else {
      namctl <- names(control)
      namdef <- names(defctrl)
      if (! all(namctl %in% namdef) ) stop("unrecognized control names")
      for (onename in namctl) {defctrl[onename]<-control[onename]}
      control<-defctrl
  }
  control # return the control list
}

#' summary.nlsr
#' 
#' prepare display result of \code{nlsr} computations - NOT compact output
#' 
#' The set of possible controls to set is as follows
#'
#' @param object  an object of class \code{nlsr}                           
#' @param ... additional data needed to evaluate the modeling functions
#'                          
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' 
#' @importFrom stats pt
#' 
#' @export
summary.nlsr <- function(object, ...) {
    smalltol <- .Machine$double.eps * 1000
    options(digits = 5) # 7 is default
    resname <- deparse(substitute(object))
    JJ <- object$jacobian
    res <- object$resid # ?? are these weighted or not?
    coeff <- object$coefficients
    pnames<-names(coeff)
    npar <- length(coeff)
    w <- object$weights
    nobs <- if (!is.null(w)) sum(w > 0) else length(res)
    lo <- object$lower
    if (is.null(lo)) lo <- rep( -Inf, npar)
    up <- object$upper
    if (is.null(up)) up <- rep( Inf, npar)
    mi <- object$maskidx
    mt <- rep(" ",npar) # start with all "unmasked"
    mt[mi] <- "M" # Put in the masks
    bdmsk <- rep(1, npar) # bounds and masks indicator ?? should it be 1L
    bdmsk[mi] <- 0 # masked
    ct <- rep(" ",npar) # start with all "free"
    for (i in seq_along(coeff)){
       if (lo[[i]] - coeff[[i]] > 0) {
          ct[[i]] <- "-" # lower bound violation
          if (bdmsk[[i]] == 1) bdmsk[[i]] <- -3
       } else { 
          if (coeff[[i]] - lo[[i]] < smalltol*(abs(coeff[[i]])+smalltol) ) {
             ct[[i]] <- "L" # "at" lower bound
             if (bdmsk[[i]] != 0) bdmsk[[i]] <- -3 # leave mask indication intact
          }
       }
       if (coeff[[i]] - up[[i]] > 0) {
          ct[[i]] <- "+" # upper bound violation
          if (bdmsk[[i]] == 1) bdmsk[[i]] <- -1
       } else { 
          if (up[[i]] - coeff[[i]] < smalltol*(abs(coeff[[i]])+smalltol) ) {
             ct[[i]] <- "U" # "at" upper bound
             if (bdmsk[[i]] != 0) bdmsk[[i]] <- -1 # leave mask indication intact
          }
       }
    }
    notmask<-which(bdmsk != 0)
    ss <- object$ssquares
    rdf <- nobs - length(notmask) # New way 220801
    if (rdf <= 0) {
          if (rdf < 0) { stop(paste("Inadmissible degrees of freedom =",rdf,sep='')) }
          else { resvar <- Inf }
    } else {
       resvar <- ss/rdf
    }
    if (length(notmask) < 1) {
       warning("No unmasked parameters")
       return(NULL)
    }
    Jwork<-JJ[, notmask]
    if (any(is.na(Jwork))){ stop("Bad working Jacobian") }
    else {
      dec <- svd(Jwork)
      U <- dec$u
      V <- dec$v
      Sd <- dec$d
    }
    Sdout <- SEs <- rep(NA, npar) # ?? Inf or NA. Set so it is defined and right length
    # ?? Fix here for example sqres that fails.
    if (min(Sd) <= smalltol * max(Sd)) { # singular
       XtXinv <- matrix(NA, nrow=npar, ncol=npar)
    } else {
       Sinv <- 1/Sd
       if (length(notmask) > 1)  { # 220809 use 
           VS <- crossprod(t(V), diag(Sinv))
       } else {
           VS <- V/Sinv
       }
       XtXinv <- crossprod(t(VS))
       SEs[notmask] <- sqrt(diag(XtXinv) * resvar) # SEs of masked params still NA
    }
    Sdout[notmask] <- Sd
##    cat("CHECK XtXinv:")
##    print(XtXinv)
##    dimnames(XtXinv) <- list(pnames, pnames) # ?? don't need this?
    gr <- crossprod(JJ, res)
    tstat<-rep(NA, npar)
    tstat[notmask] <- coeff[notmask]/SEs[notmask]
    tval <- coeff/SEs
    param <- cbind(coeff, SEs, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <-
        list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

# Note: We don't return formula because we may be doing nlfb summary 
#   i.e., resfn and jacfn approach  ?? But could we??
    ans <- list(residuals = res, sigma = sqrt(resvar),  
                df = c(npar, rdf), cov.unscaled = XtXinv,
                param = param, resname=resname, ssquares=ss, nobs=nobs, 
                ct=ct, mt=mt, Sd=Sdout, gr=gr, jeval=object$jeval,feval=object$feval)
                ## parameters = param)# never documented, for back-compatibility
#    attr(ans,"pkgname") <- "nlsr"
#    class(ans) <- "summary.nlsr" # ?? does this make sense
    ans
}


#' coef.nlsr
#' 
#' prepare and display result of \code{nlsr} computations
#' 
#' The set of possible controls to set is as follows
#'
#' @param object  an object of class \code{nlsr}                           
#' @param ... additional data needed to evaluate the modeling functions
#'                           Default FALSE
#'                          
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
# coef() function
coef.nlsr <- function(object, ...) {
       out <- object$coefficients
       attr(out,"pkgname")<-"nlsr"
##       invisible(out) # ?? is this needed, and if so why
       out # JN 170109
}

#' print.nlsr
#' 
#' prepare and display result of \code{nlsr} computations
#' 
#' The set of possible controls to set is as follows
#'
#' @param x  an object of class \code{nlsr}                           
#' @param ... additional data needed to evaluate the modeling functions
#'                           Default FALSE
#'                          
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
## ?? probably need to fix??
print.nlsr <- function(x, ...) {
  xname<-deparse(substitute(x))
  cat(class(x)," object:",xname,"\n") # ?? gets wrong resname -- just x
  if (inherits(x,"try-error") || is.null(x$coefficients)) {
     cat("Object has try-error or missing parameters\n")
     return(invisible(x))
  }
  xx<-summary(x) # calls summary to get information
  with(xx, { 
    pname<-dimnames(param)[[1]] # param is augmented coefficients with SEs and tstats
    npar <- dim(param)[1] # ?? previously length(coeff) 
    cat("residual sumsquares = ",ssquares," on ",nobs,"observations\n")
    cat("    after ",jeval,"   Jacobian and ",feval,"function evaluations\n")
    cat("  name     ","      coeff    ","     SE   ","   tstat  ",
        "   pval  ","   gradient  "," JSingval  ","\n")
    SEs <- param[,2]
    tstat <- param[,3]
    pval <- param[,4]
    for (i in seq_along(param[,1])){
      tmpname<-pname[i]
      if (is.null(tmpname)) {tmpname <- paste("p_",i,sep='')}
      cat(format(tmpname, width=10)," ")
      cat(format(param[[i]], digits=6, width=12))
      cat(ct[[i]],mt[[i]]," ")
      cat(format(SEs[[i]], digits=4, width=9)," ")
      cat(format(tstat[[i]], digits=4, width=9)," ")
      cat(format(pval[[i]], digits=4, width=9)," ")
      cat(format(gr[[i]], digits=4, width=10)," ")
      cat(format(Sd[[i]], digits=4, width=10)," ")
      cat("\n")
    }
  }) # remember to close with()
# return(NULL) 
  invisible(x)
}

#' pshort
#' 
#' 1 line result display of \code{nlsr} computations
#' 
#' @param x  an object of class \code{nlsr} 
#'                          
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
## ?? probably need to fix??
pshort <- function(x) {
  xname<-deparse(substitute(x))
  cat(xname," -- ss=",x$ssquares,":")
  prm<-x$coefficients
  pnam<-names(prm)
  for (i in 1:length(prm)){ cat(" ",pnam[i],"=",prm[i]) }
  cat(";",x$feval,"res/",x$jeval,"jac\n")
  invisible(x)
}


#' resid.nlsr  
#' 
#' prepare and display result of \code{nlsr} computations
#' 
#' @param object  an object of class \code{nlsr}                           
#' @param ...     additional data needed to evaluate the modeling functions
#'                          
#' @author  J C Nash nashjc _at_ uottawa.ca
#' 
#' ### remove _at_export to try to overcome NAMESPACE issue
#'  
resid.nlsr <- function(object, ...){ # ?? are these weighted?
  resids <- object$resid
  resids # so the function prints
}

#' residuals.nlsr
#' 
#' prepare and display result of \code{nlsr} computations
#' 
#' @param object  an object of class \code{nlsr}                           
#' @param ...  additional data needed to evaluate the modeling functions
#'                          
#' @author  J C Nash nashjc _at_ uottawa.ca
#' @export
#' 
## ?? can we do a generic resids?? ?? do we want dot args
## residuals.nlsr <- function(object, ...){ resid(object, ...) }
residuals.nlsr <- function(object, ...){
  resids <- object$resid
  resids # so the function prints
}


#' predict.nlsr
#' 
#' prepare and display predictions from an \code{nlsr} model
#'
#' @param object  an object of class \code{nlsr}                           
#' @param newdata a dataframe containing columns that match the original dataframe
#'      used to estimate the nonlinear model in the \code{nlsr} object
#' @param ... additional data needed to evaluate the modeling functions
#'                           Default FALSE
#'                          
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
predict.nlsr <- function(object=NULL, newdata=list(), ...) { 
#  This ONLY works if we have used nlxb. Do we want to add class 'nlxb'??
    if( is.null(object) ) stop("predict.nlsr REQUIRES an nlsr solution object")
    form <- object$formula
    if (is.null(form)) stop("nlsr.predict works only if formula is defined")
# ?? give more output
#
#  we assume a formula of style y~something, and use the something
#  In some ways need to check this more carefully
    env4coefs <- list2env(as.list(object$coefficients))
    preds <- eval(form[[3]], as.list(newdata), env4coefs)
    class(preds)<- "predict.nlsr"
    attr(preds,"pkgname") <- "nlsr"
    preds
}

#' fitted.nlsr
#' 
#' prepare and display fits of \code{nlsr} computations
#' 
#' @param object  an object of class \code{nlsr}
#' @param data  a data frame with the date for which fits are wanted
#'           ?? how do we guarantee it is the data used to fit the model
#'           Or do we need a different approach?
#' @param ... additional data needed to evaluate the modeling functions
#'                           Default FALSE
#'                          
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
fitted.nlsr <- function(object=NULL, data=parent.frame(), ...) { 
  #  This ONLY works if we have used nlxb. Do we want to add class 'nlxb'??
  # ?? is parent.frame() right?
  if( is.null(object) ) stop("predict.nlsr REQUIRES an nlsr solution object")
  form <- object$formula
  if (is.null(form)) stop("fitted.nlsr works only if formula is defined")
  # ?? give more output
  #
  #  we assume a formula of style y~something, and use the something
  #  In some ways need to check this more carefully
  env4coefs <- list2env(as.list(object$coefficients))
  fits <- eval(form[[3]], as.list(data), env4coefs)
  class(fits)<- "predict.nlsr" # ?? do we want all this
  attr(fits,"pkgname") <- "nlsr" # ?? do we want all this
  fits
}

#' rawres.nlsr
#' 
#' prepare and display raw residuals of \code{nlsr} computations
#' NOTE: we use model - data i.e., rhs - lhs
#' 
#' @param object  an object of class \code{nlsr}
#' @param data  a data frame with the date for which fits are wanted
#'           ?? how do we guarantee it is the data used to fit the model
#'           Or do we need a different approach?
#' @param ... additional data needed to evaluate the modeling functions
#'                           Default FALSE
#'                          
#' @return
#'    A vector of the raw residuals
#'                      
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
rawres.nlsr <- function(object=NULL, data=parent.frame(), ...) { 
  #  This ONLY works if we have used nlxb. Do we want to add class 'nlxb'??
  # ?? is parent.frame() right?
  if( is.null(object) ) stop("predict.nlsr REQUIRES an nlsr solution object")
  form <- object$formula
  if (is.null(form)) stop("fitted.nlsr works only if formula is defined")
  # ?? give more output
  #
  #  we assume a formula of style y~something, and use the something
  #  In some ways need to check this more carefully
  env4coefs <- list2env(as.list(object$coefficients))
  rawres <- eval(form[[3]]-form[[1]], as.list(data), env4coefs)
  class(rawres)<- "rawres.nlsr" # ?? do we want all this
  attr(fits,"pkgname") <- "nlsr" # ?? do we want all this
  fits
}

#' pctrl
#' 
#' Compact display of specified \code{control} vector
#' 
#' @param control  a control list
#'                          
#' @return
#'    none (??may want to change that)
#'                      
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
pctrl <- function(control) { 
   nprow <- 4 # default to 4 per row. May change later.??
   nc <- length(control) # note some items may be multiple
   namctl<-names(control)
   for (ii in 1:nc) {
      item <- control[[ii]]
      if (length(item) > 1) item <- "vec/list"
      cat("  ",namctl[[ii]],"=",item)
      if ((ii %% nprow) == 0) cat("\n")
   }
   if ((ii %% nprow) != 0) cat("\n")
   invisible(1)
}

#' pnlm0
#' 
#' Compact display of specified \code{minpack.lm} object \code{x}
#' 
#' @param x  a minpack.lm nlsLM() or nls.lm() result object
#'                          
#' @return
#'    none 
#'                      
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @importFrom stats coef
#' @importFrom stats deviance
#' @export
# pnlm0.R -- a 1 line summary of an nls() result object
pnlm0 <- function(x){
  xname <- deparse(substitute(x))
  ss<-deviance(x)
  prm<-coef(x)
  pnam <- names(prm)   
  cat(xname," -- ss=",ss,":")
  for (i in 1:length(prm)){cat(" ",pnam[i],"=",prm[i])}
  cat(";",x$convInfo$finIter," itns\n")
  invisible(1)
}

#' pnls0
#' 
#' Compact display of specified \code{nls} object \code{x}
#' 
#' @param x  an nls() result object
#'                          
#' @return
#'    none 
#'                      
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
# pnls0.R -- a 1 line summary of an nls() result object
pnls0 <- function(x){
  xname <- deparse(substitute(x))
  ss<-deviance(x)
  prm<-coef(x)
  pnam <- names(prm)   
  cat(xname," -- ss=",ss,":")
  for (i in 1:length(prm)){cat(" ",pnam[i],"=",prm[i])}
  cat(";",x$convInfo$finIter," itns\n")
  invisible(1)
}

#' nvec
#' 
#' Compact display of specified \code{nvec} named vector
#' 
#' @param vec  a named vector of parameters
#'                          
#' @return
#'    none (??may want to change that)
#'                      
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' @export
nvec <- function(vec){
   vnam<-deparse(substitute(vec))
   cat(vnam,":")
   pnam <- names(vec)
   for (i in 1:length(pnam)){ cat(pnam[i],"=",as.numeric(vec[i])," ")}
   cat("\n")
   invisible(1)
}
