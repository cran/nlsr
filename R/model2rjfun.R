#' model2rjfun
#'
#' @description 
#' These functions create functions to evaluate residuals or sums of squares at 
#' particular parameter locations.
#'
#' @aliases model2ssgrfun modelexpr SSmod2rjfun
#' @usage 
#' model2rjfun(modelformula, pvec, data = NULL, jacobian = TRUE, testresult = TRUE, ...)
#' SSmod2rjfun(modelformula, pvec, data = NULL, jacobian = TRUE, testresult = TRUE, ...)
#' model2ssgrfun(modelformula, pvec, data = NULL, gradient = TRUE, 
#'               testresult = TRUE, ...)
#' modelexpr(fun)
#'
#' @param modelformula  A formula describing a nonlinear regression model.
#' @param pvec A vector of parameters.
#' @param data A dataframe, list or environment holding data used in the calculation.
#' @param jacobian Whether to compute the Jacobian matrix.
#' @param gradient Whether to compute the gradient vector.
#' @param testresult Whether to test the function by evaluating it at \code{pvec}.
#' @param fun A function produced by one of \code{model2rjfun} or \code{model2ssgrfun}.
#' @param \dots Dot arguments, that is, arguments that may be supplied by \code{name = value}
#'    to supply information needed to compute specific quantities in the model.
#'    
#'  @details 
#'  If \code{pvec} does not have names, the parameters will have names
#'  generated in the form \samp{p_<n>}, e.g. \code{p_1, p_2}.  Names that appear in
#'  \code{pvec} will be taken to be parameters of the model.  
#'  
#'  The \code{data} argument may be a dataframe, list or environment, or \code{NULL}.
#'  If it is not an environment, one will be constructed using the components
#'  of \code{data} with parent environment set to be
#'  the environment of \code{modelformula}.  
#'  
#'  \code{SSmod2rjfun} returns a function with header \code{function(prm)}, which
#'  evaluates the residuals (and if \code{jacobian} is \code{TRUE} the
#'  Jacobian matrix) of the selfStart model (the rhs is used) at \code{prm}.  
#'  The residuals are defined to be the right hand side of \code{modelformula} 
#'  minus the left hand side. Note that the selfStart model used in the model
#'  formula must be available (i.e., loaded). If this function is called from 
#'  \code{nlxb()} then the \code{control} element \code{japprox} must be 
#   \code{"SSJac"}.
#'  
#'  @section value
#'  \code{model2rjfun} returns a function with header \code{function(prm)}, which
#'  evaluates the residuals (and if \code{jacobian} is \code{TRUE} the
#'  Jacobian matrix) of the model at \code{prm}.  The residuals are defined to be
#'  the right hand side of \code{modelformula} minus the left hand side.
#'  
#'  \code{model2ssgrfun} returns a function with header \code{function(prm)}, which
#'  evaluates the sum of squared residuals (and if \code{gradient} is \code{TRUE} the
#'  gradient vector) of the model at \code{prm}. 
#'  
#'  \code{modelexpr} returns the expression used to calculate the vector of 
#'  residuals (and possibly the Jacobian) used in the previous functions.
#'  
#' @author John Nash and Duncan Murdoch
#'
#' @importFrom stats deriv
#'      
#' @importFrom utils str
#'
#' @seealso \code{\link{nls}}
#' 
#' @examples 
#' y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 
#'         50.156, 62.948, 75.995, 91.972) 
#'         
#' tt <- seq_along(y)  # for testing
#' 
#' mydata <- data.frame(y = y, tt = tt)
#' f <- y ~ b1/(1 + b2 * exp(-1 * b3 * tt))
#' p <- c(b1 = 1, b2 = 1, b3 = 1)
#' rjfn <- model2rjfun(f, p, data = mydata)
#' rjfn(p)
#' myexp <- modelexpr(rjfn)
#' cat("myexp:"); print(myexp)
#' 
#' ssgrfn <- model2ssgrfun(f, p, data = mydata)
#' ssgrfn(p)
#' 
#' @section keyword:
#'    nonlinear 
#'
#'
#' @export
#' 
model2rjfun <- function(modelformula, pvec, data = NULL, jacobian = TRUE, 
                          testresult = TRUE, ...) {
# add dots 161127, but they have NOT been tested ??
##    cat("model2rjfun: modelformula = ")
##    print(modelformula)
##    print(class(modelformula))

    stopifnot(inherits(modelformula, "formula"))

    if (length(modelformula) == 2) { ## ?? Explain this in vignette!
        residexpr <- modelformula[[2]]
    } else if (length(modelformula) == 3) {
        residexpr <- call("-", modelformula[[3]], modelformula[[2]])
    } else stop("Unrecognized formula")
    
    if (is.null(names(pvec))){
          cat("apparently no names in pvec\n")
	  names(pvec) <- paste0("p_", seq_along(pvec))
    }
    
    if (jacobian)
        residexpr <- deriv(residexpr, names(pvec))
## SHOULD TRY:??
##	residexpr <- fnDeriv(residexpr, names(pvec))
#    data <- list(...)	
    if (is.null(data)) {
	data <- environment(modelformula) # this will handle variables in the parent frame
        }
        ## ?? Don't yet handle variable in dot args. But no dot args here.
    else if (is.list(data))
	data <- list2env(data, parent = environment(modelformula))
    else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")
    rjfun <- function(prm) {
        if (is.null(names(prm))) 
	    names(prm) <- names(pvec)
	  localdata <- list2env(as.list(prm), parent = data)
	  eval(residexpr, envir = localdata)  
        # Saves Jacobian matrix as "gradient" attribute (consistent with deriv())
    }
    
    if (testresult) {
	resids <- rjfun(pvec)
	if (any(bad <- !is.finite(resids))) 
	    stop("residuals contain ", unique(resids[bad]))
	if (jacobian && any(bad <- !is.finite(attr(resids, "gradient"))))
	    stop("Jacobian contains ", unique(attr(resids, "gradient")[bad]))
	rm(resids, bad)  # Don't want to capture these in the environment of rjfun
    }
    rjfun
}

#' @export
#' 
model2ssgrfun <- function(modelformula, pvec, data = NULL, gradient = TRUE, 
                          testresult = TRUE, ...) {
    rjfun <- model2rjfun(modelformula, pvec, data = data, jacobian = gradient, 
                         testresult = testresult, ...)
			
    function(prm) {
    	resids <- rjfun(prm)
	    ss <- as.numeric(crossprod(resids))
	    if (gradient) {
	        jacval <- attr(resids, "gradient")
	        grval <- 2*as.numeric(crossprod(jacval, resids))
	        attr(ss, "gradient") <- grval
	    }
	    attr(ss, "resids") <- resids
	    ss
    }
}

#' @export
#' 
modelexpr <- function(fun) {
    env <- environment(fun)
    if (exists("rjfun", env))
	env <- environment(env$rjfun)
    env$residexpr
}


#' @export
#' 

SSmod2rjfun <- function(modelformula, pvec, data = NULL, jacobian = TRUE, 
                        testresult = TRUE, ...) {
  stopifnot(inherits(modelformula, "formula"))
  
  if (length(modelformula) == 2) { ## ?? Explain this in vignette!
    residexpr <- modelformula[[2]]
    rhs <- residexpr # for SS case
  } else if (length(modelformula) == 3) {
    residexpr <- call("-", modelformula[[3]], modelformula[[2]])
    rhs <- modelformula[[3]] # for SS case
  } else stop("Unrecognized formula")

  if (is.null(names(pvec))){
    cat("apparently no names in pvec\n")
    names(pvec) <- paste0("p_", seq_along(pvec))
  }
  
  if (jacobian)
  if (is.null(data)) {
    data <- environment(modelformula) # this will handle variables in the parent frame
  }
  ## ?? Don't yet handle variable in dot args. But no dot args here.
  else if (is.list(data))
    data <- list2env(data, parent = environment(modelformula))
  else if (!is.environment(data))
    stop("'data' must be a dataframe, list, or environment")
  rjfun <- function(prm) {
    if (is.null(names(prm))) 
      names(prm) <- names(pvec)
    localdata <- list2env(as.list(prm), parent = data)
    eval(residexpr, envir = localdata)  
    # Saves Jacobian matrix as "gradient" attribute (consistent with deriv())
  }
  
  if (testresult) {
    resids <- rjfun(pvec)
    if (any(bad <- !is.finite(resids))) 
      stop("residuals contain ", unique(resids[bad]))
    if (jacobian && any(bad <- !is.finite(attr(resids, "gradient"))))
      stop("Jacobian contains ", unique(attr(resids, "gradient")[bad]))
    rm(resids, bad)  # Don't want to capture these in the environment of rjfun
  }
  rjfun
}
