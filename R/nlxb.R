#' nlxb: nonlinear least squares modeling by formula
#' 
#' A simplified and hopefully robust alternative to finding
#' the nonlinear least squares minimizer that causes
#' 'formula' to give a minimal residual sum of squares.
#' 
#' nlxb is particularly intended to allow for the
#' resolution of very ill-conditioned or else near
#' zero-residual problems for which the regular nls()
#' function is ill-suited. 
#' 
#' This variant uses a qr solution without forming the sum
#' of squares and cross products t(J)%*%J
#' 
#' Neither this function nor \code{nlfb} have provision for parameter
#' scaling (as in the \code{parscale} control of \code{optim} and
#' package \code{optimx}). This would be more tedious than difficult to
#' introduce, but does not seem to be a priority feature to add.
#' 
#' @param formula The modeling formula. Looks like 'y~b1/(1+b2*exp(-b3*T))' 
#' @param start MUST be a named vector with all elements present
#'     e.g., start=c(b1=200, b2=50, b3=0.3) 
#' @param trace TRUE for console output during execution
#' @param data a data frame containing data for variables
#'     used in the formula that are NOT the parameters. This
#'     data may also be defined in the parent frame i.e.,
#'     'global' to this function 
#' @param lower a vector of lower bounds on the parameters. 
#'         If a single number, this will be applied to all parameters
#'         Default \code{NULL}.
#' @param upper a vector of upper bounds on the parameters. If a single number, 
#'     this will be applied to all parameters. Default \code{NULL}.
## #' @param masked a character vector of quoted parameter names. These parameters 
## #'    will NOT be altered by the algorithm. Masks may also be defined by setting 
## #'    lower and upper bounds equal for the parameters to be fixed. Note that the 
## #'    starting parameter value must also be the same as the lower and upper bound 
## #'    value.
#' @param weights A vector of fixed weights. The objective function that will be 
#'    minimized is the sum of squares where each residual is multiplied by the 
#'    square root of the corresponding weight. Default \code{NULL} implies 
#'    unit weights.
#' @param control a list of control parameters. See nlsr.control().
#' @param ... additional data needed to evaluate the modeling functions
#' 
#' @details 
#'   There are many controls, and some of them are important for \code{nlxb}.
#'   In particular, if the derivatives needed for developing the Jacobian are
#'   NOT in the derivatives table, then we must supply code elsewhere as 
#'   specified by the control \code{japprox}. This was originally just for
#'   numerical approximations, with the character strings "jafwd", "jaback",
#'   "jacentral" and "jand" leading to the use of a forward, backward, central
#'   or package \code{numDeriv} approximation. However, it is also possible to
#'   use code embedded in the residual function created using the \code{formula}.
#'   This is particularly useful for \code{selfStart} models, and we use the
#'   character string "SSJac" to point to such Jacobian code. Note how the
#'   starting parameter vector is found using the \code{getInitial} function
#'   from the \code{stats} package as in an example below.
#' 
#' @returns 
#' 
#'   list of solution elements
#'   
#'   resid    = weighted residuals at the proposed solution
#'   jacobian = Jacobian matrix at the proposed solution
#'   feval    = residual function evaluations used to reach solution from starting parameters
#'   jeval    = Jacobian function (or approximation) evaluations used
#'   coefficients = a named vector of proposed solution parameters
#'   ssquares = weighted sum of squared residuals (often the deviance)
#'   lower    = lower bounds on parameters
#'   upper    = upper bounds on parameters
#'   maskidx  = vector if indices of fixed (masked) parameters
#'   weights  = specified weights on observations
#'   formula  = the modeling formula
#'   resfn    = the residual function (unweighted) based on the formula
#' 
#' @author  J C Nash 2014-7-16   nashjc _at_ uottawa.ca
#' 
#' @importFrom stats as.formula
#' @importFrom stats setNames
#' 
#' @examples 
#' library(nlsr)
#' weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
#'         38.558, 50.156, 62.948, 75.995, 91.972)
#' tt <- 1:12
#' weeddf <- data.frame(tt, weed)
#' frm <- 
#' wmodu <- weed ~ b1/(1+b2*exp(-b3*tt)) # Unscaled
#' ## nls from unit start FAILS
#' start1<-c(b1=1, b2=1, b3=1)
#' hunls1 <- try(nls(wmodu, data=weeddf, start=start1, trace=FALSE))
#' if (! inherits(hunls1, "try-error")) print(hunls1) ## else cat("Failure -- try-error\n")
#' ## nlxb from unit start
#' hunlx1 <- try(nlxb(wmodu, data=weeddf, start=c(b1=1, b2=1, b3=1), trace=FALSE))
#' if (! inherits(hunlx1, "try-error")) print(hunlx1)
#' 
#' st2h<-c(b1=185, b2=10, b3=.3)
#' #' hunls2 <- try(nls(wmodu, data=weeddf, start=st2h, trace=FALSE))
#' if (! inherits(hunls1, "try-error")) print(hunls1) ## else cat("Failure -- try-error\n")
#' ## nlxb from unit start
#' hunlx1 <- try(nlxb(wmodu, data=weeddf, start=st2h, trace=FALSE))
#' if (! inherits(hunlx1, "try-error")) print(hunlx1)
#' 
#' # Functional models need to use a Jacobian approximation or external calculation.
#' # For example, the SSlogis() selfStart model from \code{stats} package.
#' 
#' # nls() needs NO starting value
#' hSSnls<-try(nls(weed~SSlogis(tt, Asym, xmid, scal), data=weeddf))
#' summary(hSSnls)
#' # We need to get the start for nlxb explicitly
#' stSS <- getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
#' hSSnlx<-try(nlxb(weed~SSlogis(tt, Asym, xmid, scal), data=weeddf, start=stSS))
#' hSSnlx
#' 
#' # nls() can only bound parameters with algorithm="port"
#' #   and minpack.lm is unreliable in imposing bounds, but nlsr copes fine.
#' lo<-c(0, 0, 0)
#' up<-c(190, 10, 2) # Note: start must be admissible.
#' bnls0<-try(nls(wmodu, data=weeddf, start=st2h,
#'          lower=lo, upper=up)) # should complain and fail
#'  
#' bnls<-try(nls(wmodu, data=weeddf, start=st2h,
#'          lower=lo, upper=up, algorith="port"))
#' summary(bnls)
#' bnlx<-try(nlxb(wmodu, data=weeddf, start=st2h, lower=lo, upper=up))
#' bnlx
#' 
#' # nlxb() can also MASK (fix) parameters. The mechanism of maskidx from nls
#' # is NO LONGER USED. Instead we set upper and lower parameters equal for
#' # the masked parameters. The start value MUST be equal to this fixed value.
#' lo<-c(190, 0, 0) # mask first parameter
#' up<-c(190, 10, 2)
#' strt <- c(b1=190, b2=1, b3=1)
#' mnlx<-try(nlxb(wmodu, start=strt, data=weeddf, 
#'          lower=lo, upper=up))
#' mnlx
#' mnls<-try(nls(wmodu, data=weeddf, start=strt,
#'          lower=lo, upper=up, algorith="port"))
#' summary(mnls)
#'
#' # Try first parameter masked and see if we get SEs 
#' lo<-c(200, 0, 0) # mask first parameter
#' up<-c(100, 10, 5)
#' strt <- c(b1=200, b2=1, b3=1)
#' mnlx<-try(nlxb(wmodu, start=strt, data=weeddf, 
#'          lower=lo, upper=up))
#' mnlx
#' mnls<-try(nls(wmodu, data=weeddf, start=strt,
#'          lower=lo, upper=up, algorith="port"))
#' summary(mnls) 
#' 
#' 
#' @export
nlxb <- function(formula, data=parent.frame(), start, trace = FALSE,  lower = NULL,
                 upper = NULL, weights=NULL, control=list(), ...) {
    if("subset" %in% names(list(...)))
        stop("subset NOT used in nlxb(). Please explicitly subset your data.")
    ctrl <- nlsr.control(control)
    if (missing(start)) stop("A start vector MUST be provided")
    pnames <- names(start)
    start <- as.numeric(start) # ensure we convert (e.g., if matrix)
    names(start) <- pnames ## as.numeric strips names, so this is needed
    # bounds
    npar <- length(start)  # number of parameters
    if (is.null(lower)) lower<- -Inf
    if (is.null(upper)) upper<- Inf
    if (length(lower) == 1) 
        lower <- rep(lower, npar) # expand to full dimension
    if (length(upper) == 1) 
        upper <- rep(upper, npar)
# more tests on bounds
    if (length(lower) != npar) 
        stop("Wrong length: lower")
    if (length(upper) != npar) 
        stop("Wrong length: upper")
    if (any(start < lower) || any(start > upper)) {
        if(trace) {  cat("start:"); print(start)
                     cat("\n");cat("lower:"); print(lower)
                     cat("\n"); cat("upper:"); print(upper) 
                  }
        stop("Infeasible start")
    }
    if (trace && ctrl$prtlvl>1) {
        cat("formula: ")
        print(formula)
        cat("lower:"); print(lower); cat("upper:"); print(upper) 
    }
##      maxlamda <- 1e+60) ## dropped 130709 ??why?
##    epstol <- (.Machine$double.eps) * ctrl$offset # ??161018 - not used elsewhere
    if (trace && ctrl$prtlvl>2) print(str(ctrl))
# Note spelling of lamda -- a throwback to Ag Can 1974 and way to see if folk are copying code.
# First get all the variable names:
#    vn <- all.vars(parse(text = formula))
# ??? need to fix -- why?, what is wrong
    vn <- all.vars(formula)
    # Then see which ones are parameters (get their positions
    # in the set xx
    pnum <- start  # may simplify later??
    pnames <- names(pnum)
    bdmsk <- rep(1, npar)  # set all params free for now
    maskidx <- which(lower==upper)
    if (length(maskidx) > 0 && trace) {
        cat("The following parameters are masked:")
        print(pnames[maskidx])
    }
    bdmsk[maskidx] <- 0  # fixed parameters
    if (trace && ctrl$prtlvl>1) { # diagnostic printout
        cat("Finished masks check\n")
        parpos <- match(pnames, vn) # ?? check this is right??
        datvar <- vn[-parpos]  # NOT the parameters
        cat("datvar:")
        print(datvar)
        for (i in 1:length(datvar)) {
            dvn <- datvar[[i]]
            cat("Data variable ", dvn, ":")
            if (is.null(data)) { 
                print(eval(parse(text = dvn)))
            } else {
                print(with(data, eval(parse(text = dvn)))) 
            }
        }
    }
    if (is.null(ctrl$japprox)) {
       trjfn<-model2rjfun(formula, pnum, data=data) 
    } else {
      if(ctrl$japprox == "SSJac") {
         trjfn<-SSmod2rjfun(formula, pnum, data=data)
      } else {
         trjfn<-model2rjfun(formula, pnum, data=data, jacobian=FALSE)
      }
    } 
    if (trace && ctrl$prtlvl>2) { cat("str(trjfn):\n"); print(str(trjfn))}
    if (trace && ctrl$prtlvl>2){ cat("ctrl:\n"); print(ctrl) }
    if (is.null(ctrl$japprox)) { # call regularly when no jacobian approx
       resfb <- nlfb(start=pnum, resfn=trjfn, jacfn=trjfn, trace=trace, 
            data=data, lower=lower, upper=upper, 
            weights=weights, ctrlcopy=TRUE, control=ctrl)
    } else { # need approx
      if(!is.character(ctrl$japprox)) stop("Non-character ctrl$japprox")
      if(ctrl$japprox == "SSJac") { # Jacobian from selfStart model
        if (trace) { cat('Using Jacobian code from selfStart model \n') }
       resfb <- nlfb(start=pnum, resfn=trjfn, jacfn=trjfn, trace=trace, 
            data=data, lower=lower, upper=upper, 
            weights=weights, ctrlcopy=TRUE, control=ctrl)
      } else {
      if (trace) { cat('Using approximation ',as.character(ctrl$japprox),' 
		 is.character(ctrl$japprox)=', is.character(ctrl$japprox),'\n') }
       resfb <- nlfb(start=pnum, resfn=trjfn, jacfn=ctrl$japprox, trace=trace, 
            data=data, lower=lower, upper=upper, 
            weights=weights, ctrlcopy=TRUE, control=ctrl)
      }
    }
    resfb$formula <- formula # 190805 to add formula to solution (??for predict?)
    resfb$resfn <- trjfn # May be lacking Jacobian
# ?? should there be any ... arguments
    pnum <- as.vector(resfb$coefficients)
    names(pnum) <- pnames # Make sure names re-attached. ??Is this needed??
##    resfb$coefficients <- pnum ## commented 190821
    result <- resfb
##    attr(result, "pkgname") <- "nlsr"
    class(result) <- "nlsr" ## needed for print method
    result
}
