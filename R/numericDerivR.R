#' numericDerivR: numerically evaluates the gradient of an expression. All in R
#' 
#' This version is all in R to replace the C version in package stats
#' @param expr expression or call to be differentiated. Should evaluate to a numeric vector.
#' @param theta	character vector of names of numeric variables used in expr.
#' @param rho environment containing all the variables needed to evaluate expr.
#' @param dir numeric vector of directions, typically with values in -1, 1 to use for the finite 
#'    differences; will be recycled to the length of theta.
#' @param eps	a positive number, to be used as unit step size hh for the approximate 
#'     numerical derivative (f(x+h)-f(x))/h (f(x+h)-f(x))/h or the central version, see central.
#' @param central logical indicating if central divided differences should be computed, 
#'  i.e., (f(x+h) - f(x-h)) / 2h (f(x+h)-f(x-h))/2h. These are typically more accurate but 
#'  need more evaluations of f()f().
#' @return 
#' The value of eval(expr, envir = rho) plus a matrix attribute "gradient". The columns of 
#' this matrix are the derivatives of the value with respect to the variables listed in theta.
#' @examples 
#' ex <- expression(a/(1+b*exp(-tt*c)) - weed)
#' weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
#'               38.558, 50.156, 62.948, 75.995, 91.972)
#' tt <- 1:12
#' a <- 200; b <- 50; c <- 0.3
#' dhobb <- numericDerivR(ex, theta=c("a", "b", "c"))
#' print(dhobb)
#' # exf <- ~ a/(1+b*exp(-tt*c)) - weed
#' # Note that a formula doesn't work
#' # dh1 <-  try(numericDerivR(exf, theta=c("a", "b", "c")))
#' @export
numericDerivR <- function(expr, theta, rho = parent.frame(), dir = 1,
                 eps = .Machine$double.eps ^ (1/if(central) 3 else 2), central = FALSE)
## Note: this expr must be set up as a call to work properly according to JN??
## ?? we set eps conditional on central. But central set AFTER eps. Is this OK.
{   ## cat("numericDeriv-Alt\n")
    dir <- rep_len(dir, length(theta))
    stopifnot(is.finite(eps), eps > 0)
##    rho1 <- new.env(FALSE, rho, 0)
    if (!is.character(theta) ) {stop("'theta' should be of type character")}
    if (is.null(rho)) {
            stop("use of NULL environment is defunct")
            #        rho <- R_BaseEnv;
    } else {
          if(! is.environment(rho)) {stop("'rho' should be an environment")}
          #    int nprot = 3;
    }
    if( ! ((length(dir) == length(theta) ) & (is.numeric(dir) ) ) )
              {stop("'dir' is not a numeric vector of the correct length") }
    if(is.na(central)) { stop("'central' is NA, but must be TRUE or FALSE") }
    res0 <- eval(expr, rho) # the base residuals. ?? C has a check for REAL ANS=res0
#    cat("res0:"); print(res0)
##    if (any(is.infinite(res0)) ) {stop("residuals cannot be evaluated at base point")}
## 220715 -- let it fail!
    ##  CHECK_FN_VAL(res, ans);  ?? how to do this. Is it necessary?
    nt <- length(theta) # number of parameters
    mr <- length(res0) # number of residuals
    JJ <- matrix(NA, nrow=mr, ncol=nt) # Initialize the Jacobian
    for (j in 1:nt){
       origPar<-get(theta[j],rho) # This is parameter by NAME, not index
       xx <- abs(origPar)
       delta <- if (xx == 0.0) {eps} else { xx*eps }
       ## JN: I prefer eps*(xx + eps)  which is simpler ?? Should we suggest / use a control switch
       prmx<-origPar+delta*dir[j]
       assign(theta[j],prmx,rho)
       res1 <- eval(expr, rho) # new residuals (forward step)
#       cat("res1:"); print(res1)
       if (central) { # compute backward step resids for central diff
          prmb <- origPar - dir[j]*delta
          assign(theta[j], prmb, envir=rho) # may be able to make more efficient later??
          resb <- eval(expr, rho)
          JJ[, j] <- dir[j]*(res1-resb)/(2*delta) # vectorized
       } else { ## forward diff
          JJ[, j] <- dir[j]*(res1-res0)/delta
       }  # end forward diff
       assign(theta[j],origPar, rho) # reset value ??
    } # end loop over the parameters
    attr(res0, "gradient") <- JJ
    return(res0)
}
