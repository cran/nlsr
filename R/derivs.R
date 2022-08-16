# R-based replacement for deriv() function
#' dex
#' 
#' Calculate expression for derivative calculations.
#' Converts input to an expression suitable for use
#' in \code{\link{nlsDeriv}} and related functions.
#' 
#' @usage
#'   dex(x, do_substitute = NA, verbose = FALSE)
#'   
#' @param x An expression represented in a variety of ways.
#'      See Details.
#'      
#' @param  do_substitute Whether to use the expression passed as \code{x}, or
#'      to evaluate it and use its value.  
#'      
#' @param  verbose Print messages describing the process.
#'      
#' @section Details:
#'    If \code{do_substitute} is \code{NA}, the following 
#'      rules are used:
#'      
#'      An attempt is made to evaluate \code{x}.  If that fails, the expression is used.
#'      
#'      If the evaluation succeeds and the value is a character vector, it is parsed.
#'      
#'      If the value is not a character vector and the expression is a single name, the value is used.
#'      
#'      Otherwise, the expression is used.
#'      
#'      Once the expression is determined it may be simplified,
#'      by extracting the language object from a length-one
#'      expression vector, or the right-hand-side from a 
#'      formula.
#'      
#'      Normally a warning will be issued if \code{x} is a formula
#'      containing a left-hand side.  To suppress this, 
#'      wrap the formula in \code{expression()}, or pass it 
#'      as a character string to be parsed.
#'      
#' @section value:
#'      An expression or language object suitable as input
#'      to \code{\link{nlsDeriv}} and related functions.
#'      
#' @author Duncan Murdoch
#'  
#' @examples 
#'    aa <- dex(~ x^2)
#'    aa
#'    str(aa)
#'    bb <- dex(expression(x^2))
#'    bb
#'    str(bb)
#'    cc <- dex("x^2")
#'    cc
#'    str(cc)
#'    
#' @export
#'    
dex <- function(x, do_substitute = NA, verbose = FALSE) {
  expr <- substitute(x)
  if (is.na(do_substitute) || !do_substitute)
    value <- tryCatch(x, error = function(e) e)
  if (is.na(do_substitute)) {
    if (verbose)
      message("Determining whether to use expression or value...")
    if (inherits(value, "error"))
      do_substitute <- TRUE
    else if (is.character(value))
      do_substitute <- FALSE
    else if (is.name(expr))
      do_substitute <- FALSE
    else
      do_substitute <- TRUE
  } else if (!do_substitute && inherits(value, "error"))
    stop(conditionMessage(value), call. = FALSE)
  if (verbose) {
    if (do_substitute) 
      message("Using expression.")
    else message("Using value.")
  }  
  if (!do_substitute) {
    if (is.character(value)) {
      if (verbose) 
        message("Parsing value.")
      expr <- parse(text=value)
    } else
      expr <- value
  }
  if (is.expression(expr) && length(expr) == 1)
    expr <- expr[[1]]
  else if (is.call(expr) && as.character(expr[[1]]) == "~") {
    if (length(expr) == 3) {
      warning("Left hand side of formula will be ignored.")
      expr <- expr[[3]]
    } else
      expr <- expr[[2]]
  }
  expr
}

#' sysDerivs
#' 
#' creates a new environment whose parent is emptyenv
#' 
#' @usage  sysDerivs()
#'  
#' @export
sysDerivs <- new.env(parent = emptyenv())

#' sysSimplifications
#' 
#' creates a new environment whose parent is emptyenv
#' 
#' @export
sysSimplifications <- new.env(parent = emptyenv())

#' newDeriv
#' 
#' Define a new derivative expression
#' 
#' @usage newDeriv (expr, deriv, derivEnv = sysDerivs) 
#'   
#' @param expr An expression represented in a variety of ways.
#' @param deriv   An expression giving the derivative of the function call in \code{expr}.
#' @param derivEnv The environment in which to evaluate (??) these
#'      
#' @export
newDeriv <- function(expr, deriv, derivEnv = sysDerivs) {
    if (missing(expr))
    	return(ls(derivEnv)) # ?? should we not throw an error?
    expr <- substitute(expr) # Get parse tree of expr 
    if (!is.call(expr))
    	stop("expr must be a call to a function")
    fn <- as.character(expr[[1]])
    if (missing(deriv)) 
    	return(derivEnv[[fn]])
    deriv <- substitute(deriv)
    args <- expr[-1]
    argnames <- names(args)
    if (is.null(argnames))
    	argnames <- rep("", length(args))
    required <- which(argnames == "")
    for (i in required)
    	argnames[i] <- as.character(args[[i]])
    value <- list(expr = expr, argnames = argnames, 
                  required = required, deriv = deriv)
    if (!is.null(oldval <- derivEnv[[fn]]) && !identical(value, oldval))
      warning(gettextf("changed derivative for %s", dQuote(fn)))
    assign(fn, value, envir = derivEnv)
    invisible(value) # returns but does not print ?? why not?
} 

#' newSimplification
#' 
#' Define a new simplification expression
#' 
#' @usage newSimplification(expr, test, simplification, do_eval = FALSE, 
#'                        simpEnv = sysSimplifications)
#'   
#' @param expr An expression represented in a variety of ways.
#' @param test ?? what is this
#' @param simplification ?? what is this
#' @param do_eval ??evaluate the result
#' @param simpEnv the environment in which the simplification is carried out
#'      
#' @export
newSimplification <- function(expr, test, simplification, 
                              do_eval = FALSE, simpEnv = sysSimplifications) {
    if (missing(expr))
    	return(ls(simpEnv))
    expr <- substitute(expr)
    if (!is.call(expr))
    	stop("expr must be a call to a function")
    fn <- as.character(expr[[1]])
    nargs <- length(expr) - 1L
    simps <- simpEnv[[fn]]
    if (missing(test)) {
        if (nargs <= length(simps))	
    	    return(simps[[nargs]])
    	else
    	    return(NULL)
    }
    test <- substitute(test)
    simplification <- substitute(simplification)
    
    args <- expr[-1]
    if (!is.null(names(args)))
    	stop("expr should not have named arguments")
    if (!all(sapply(args, is.name)))
    	stop("expr should have simple names as arguments")
    argnames <- sapply(args, as.character)
    if (any(duplicated(argnames)))
    	stop("expr names should be unique")
    	
    if (is.null(simps)) simps <- list()
    if (nargs <= length(simps)) 
    	simpn <- simps[[nargs]]
    else
    	simpn <- list()
    simpn <- c(simpn, list(list(expr = expr, argnames = argnames, test = test, 
                                simplification = simplification, do_eval = do_eval)))
    simps[[nargs]] <- simpn
    assign(fn, simps, envir = simpEnv)
}
    	

#' nlsDeriv
#'   Functions to take symbolic derivatives.
#' 
#' @aliases codeDeriv fnDeriv
#'   
#' @description 
#'   
#' Compute derivatives of simple expressions symbolically, allowing user-specified derivatives.
#' 
#' @usage  nlsDeriv(expr, name, derivEnv = sysDerivs, do_substitute = FALSE, verbose = FALSE, ...)
#'            
#' @usage  codeDeriv(expr, namevec, hessian = FALSE, derivEnv = sysDerivs, 
#'           do_substitute = FALSE, verbose = FALSE, ...) 
#'           
#' @usage  fnDeriv(expr, namevec, args = all.vars(expr), env = environment(expr), 
#'           do_substitute = FALSE, verbose = FALSE, ...)
#'           
#' @param expr  An expression represented in a variety of ways. See Details.
#'  
#' @param name  The name of the variable with respect to which the derivative will be computed.
#'  
#' @param derivEnv  The environment in which derivatives are stored.
#'  
#' @param do_substitute   If \code{TRUE}, use \code{\link{substitute}} to get the expression passed as
#'        \code{expr}, otherwise evaluate it.
#'        
#' @param verbose  If \code{TRUE}, then diagnostic output will be printed as derivatives
#'      and simplifications are recognized.
#'      
#' @param ...   Additional parameters which will be passed to \code{codeDeriv}
#'       from \code{fnDeriv}, and to \code{nlsSimplify} from 
#'       \code{nlsDeriv} and \code{codeDeriv}.
#'       
#' @param namevec    Character vector giving the variable names with respect to 
#'    which the derivatives will be taken.
#'      
#' @param hessian    Logical indicator of whether the 2nd derivatives should also be computed.
#'  
#' @param args    Desired arguments for the function.  See Details below.
#'  
#' @param env The environment to be attached to the created function.  
#'        If \code{NULL}, the caller's frame is used.
#'       
#' @details
#' {   
#'  Functions \code{nlsDeriv} and \code{codeDeriv} are designed as replacements 
#'  for the \pkg{stats} package functions \code{\link{D}} and \code{\link{deriv}}
#'  respectively, though the argument lists do not match exactly.
#'  
#'  The \code{nlsDeriv} function computes a symbolic derivative of an expression
#'  or language object.  Known derivatives are stored in
#'  \code{derivEnv}; the default \code{sysDerivs} contains expressions for
#'  all of the derivatives recognized by \code{\link{deriv}}, but in
#'  addition allows differentiation with respect to any parameter
#'  where it makes sense.  It also allows the derivative of \code{abs}
#'  and \code{sign}, using an arbitrary choice of 0 at the discontinuities.
#'  
#'  The \code{codeDeriv} function computes
#'  an expression for efficient calculation of the expression value together
#'  with its gradient and optionally the Hessian matrix.
#'  
#'  The \code{fnDeriv} function wraps the \code{codeDeriv} result
#'  in a function.  If the \code{args} are given as a character
#'  vector (the default), the arguments will have those names,
#'  with no default values.  Alternatively, a custom argument list with default values can
#'  be created using \code{\link{alist}}; see the example below.
#'  
#'  The \code{expr} argument will be converted to a
#'  language object using \code{\link{dex}} (but note
#'  the different default for \code{do_substitute}).  
#'  Normally it should be a formula with no left
#'  hand side, e.g. \code{ ~ x^2 }, or an expression vector
#'  e.g. \code{ expression(x, x^2, x^3) }, or a language
#'  object e.g. \code{quote(x^2)}.  In \code{codeDeriv} and
#'  \code{fnDeriv} the expression vector must be of length 1.
#'  
#'  The \code{newDeriv} function is used to define a new derivative.
#'  The \code{expr} argument should match the header of the function as a
#'  call to it (e.g. as in the help pages), and the \code{deriv} argument
#'  should be an expression giving the derivative, including calls to
#'  \code{D(arg)}, which will not be evaluated, but will be substituted
#'  with partial derivatives of that argument with respect to \code{name}.
#'  See the examples below.  
#'  
#'  If \code{expr} or \code{deriv} is missing in a call to
#'  \code{newDeriv()}, it will return the currently saved derivative
#'  record from \code{derivEnv}.  If \code{name} is missing in a call to
#'  \code{nlsDeriv} with a function call, it will print a message describing
#'  the derivative formula and return \code{NULL}.
#'  
#'  To handle functions which act differently if a parameter is
#'  missing, code the default value of that parameter to \code{.MissingVal},
#'  and give a derivative that is conditional on \code{missing()}
#'  applied to that parameter.  See the derivatives of \code{"-"} and \code{"+"} 
#'  in the file \code{derivs.R} for an example.
#'  }
#'  
#'  @section value
#'  {
#'  If \code{expr} is an expression vector, \code{nlsDeriv} and \code{nlsSimplify}
#'  return expression vectors containing the response.  
#'  For formulas or language objects, a language object is returned.
#'  
#'  \code{codeDeriv} always returns a language object.
#'  
#'  \code{fnDeriv} returns a closure (i.e. a function).
#'  
#'  \code{nlsDeriv} returns the symbolic derivative of the expression.
#'  
#'  \code{newDeriv} with \code{expr} and \code{deriv} specified is
#'  called for the side effect of recording the derivative in \code{derivEnv}.
#'  If \code{expr} is missing, it will return the list of names of functions
#'  for which derivatives are recorded.  If \code{deriv} is missing, it
#'  will return its record for the specified function.
#'  }
#'  @section note
#'  #'  \code{newDeriv(expr, deriv, ...)} will issue a warning
#'  if a different definition for the derivative exists
#'  in the derivative table.
#'  
#'  @author Duncan Murdoch
#'  
#'  @seealso \code{\link{deriv}}
### dropped link to nlsSimplify
#'  
#'  @examples 
#'  newDeriv()
#'  newDeriv(sin(x))
#'  nlsDeriv(~ sin(x+y), "x")
#'
#'  f <- function(x) x^2
#'  newDeriv(f(x), 2*x*D(x))
#'  nlsDeriv(~ f(abs(x)), "x")
#'  
#'  nlsDeriv(~ pnorm(x, sd=2, log = TRUE), "x")
#'  fnDeriv(~ pnorm(x, sd = sd, log = TRUE), "x")
#'  f <- fnDeriv(~ pnorm(x, sd = sd, log = TRUE), "x", args = alist(x =, sd = 2))
#'  f
#'  f(1)
#'  100*(f(1.01) - f(1))  # Should be close to the gradient
#'  
#'        # The attached gradient attribute (from f(1.01)) is
#'        # meaningless after the subtraction.
#'        
#'  # Multiple point example
#'  xvals <- c(1, 3, 4.123)
#'  print(f(xvals))
#'  # Getting a hessian matrix
#'  f2 <- ~ (x-2)^3*y - y^2
#'  mydf2 <- fnDeriv(f2, c("x","y"), hessian=TRUE)
#'  # display the resulting function
#'  print(mydf2)
#'  x <- c(1, 2)
#'  y <- c(0.5, 0.1)
#'  evalmydf2 <- mydf2(x, y)
#'  print(evalmydf2)
#'  # the first index of the hessian attribute is the point at which we want the hessian
#'  hmat1 <- as.matrix(attr(evalmydf2,"hessian")[1,,])
#'  print(hmat1)
#'  hmat2 <- as.matrix(attr(evalmydf2,"hessian")[2,,])
#'  print(hmat2)
#'  
#'  @section keyword 
#'    math, nonlinear
#'    
# This is a more general version of D()
#' @export
#' 
nlsDeriv <- function(expr, name, derivEnv = sysDerivs, do_substitute = FALSE, verbose = FALSE, ...) {
    Recurse <- function(expr) {
    	if (is.call(expr)) {
    	    if (as.character(expr[[1]]) == "D")
    	    	expr <- nlsDeriv(expr[[2]], name, derivEnv, do_substitute = FALSE, verbose = verbose, ...)
    	    else
    	    	for (i in seq_along(expr)[-1])
    	    	    expr[[i]] <- Recurse(expr[[i]])
    	}
    	expr
    }
    expr <- dex(expr, do_substitute = do_substitute, verbose = verbose)
    if (is.expression(expr))
    	return(as.expression(lapply(expr, nlsDeriv, name = name, derivEnv = derivEnv, do_substitute = FALSE, verbose = verbose, ...)))
    else if (is.numeric(expr) || is.logical(expr))
    	return(0)
    else if (is.call(expr)) {
    	fn <- as.character(expr[[1]])
	if (fn == "expression")
	    return(as.expression(lapply(as.list(expr)[-1], nlsDeriv, name = name, derivEnv = derivEnv, do_substitute = FALSE, verbose = verbose, ...)))
    	model <- derivEnv[[fn]]
    	if (is.null(model))
    	    stop("no derivative known for '", fn, "'")
 	if (missing(name) || verbose) {
 	    message(paste("Expr:", deparse(expr)))
 	    message(if (missing(name)) "Pattern for" else "Using pattern")
 	    message(paste("  ", deparse(model$expr), collapse = "\n"))
 	    message("is")
 	    message(paste("  ", deparse(model$deriv), collapse = "\n"))
 	    if (missing(name))
 	    	return(invisible(NULL))
 	}
        args <- expr[-1]
        argnames <- names(args)
        if (is.null(argnames)) 
            argnames <- rep("", length(args))
        modelnames <- model$argnames
        argnum <- pmatch(argnames, modelnames)
        if (any(bad <- is.na(argnum) & argnames != ""))
            stop("Argument names not matched: ", paste(argnames[bad], collapse = ", "))
        unused <- setdiff(seq_along(modelnames), argnum)
        nonamecount <- sum(is.na(argnum))
        length(unused) <- nonamecount
        argnum[which(is.na(argnum))] <- unused
        default <- setdiff(seq_along(modelnames), argnum)
        if (length(bad <- setdiff(model$required, argnum)))
            stop("Missing required arguments: ", paste(modelnames[bad], collapse = ", "))
        
        # Now do the substitutions
        subst <- list()
        subst[argnum] <- as.list(args)
        subst[default] <- as.list(model$expr[-1])[default]
        names(subst) <- modelnames
        result <- do.call(substitute, list(model$deriv, subst))
        result <- Recurse(result)
        nlsSimplify(result, verbose = verbose, ...)
    } else if (is.name(expr))
        return( as.numeric(as.character(expr) == name) )
}

# This is a more general version of deriv(), since it allows user specified 
# derivatives and simplifications
#' @export
codeDeriv <- function(expr, namevec, 
       hessian = FALSE, derivEnv = sysDerivs, 
       do_substitute = FALSE, verbose = FALSE, ...) {
  expr <- dex(expr, do_substitute = do_substitute, verbose = verbose)
  expr <- as.expression(expr)
  if (length(expr) > 1)
    stop("Only single expressions allowed")
  exprs <- as.list(expr)
  n <- length(namevec)
  length(exprs) <- n + 1L
  for (i in seq_len(n))
    exprs[[i + 1]] <- nlsDeriv(expr[[1]], namevec[i], derivEnv = derivEnv, do_substitute = FALSE, verbose = verbose, ...)
  names(exprs) <- c(".value", namevec)
  if (hessian) {
    m <- length(exprs)
	  length(exprs) <- m + n*(n+1)/2
	  for (i in seq_len(n))
	    for (j in seq_len(n-i+1) + i-1) {
		    m <- m + 1
		    exprs[[m]] <- nlsDeriv(exprs[[i + 1]], namevec[j], derivEnv = derivEnv, 
		                    do_substitute = FALSE, verbose = verbose, ...)
	    }
  }
  exprs <- as.expression(exprs)
  subexprs <- findSubexprs(exprs)
  m <- length(subexprs)
  final <- subexprs[[m]]
  subexprs[[m]] <- substitute(.value <- expr, list(expr = final[[".value"]]))
  subexprs[[m+1]] <- substitute(.grad <- array(0, c(length(.value), namelen), list(NULL, namevec)),
				  list(namevec = namevec, namelen = length(namevec)))
  if (hessian)
    subexprs[[m+2]] <- substitute(.hessian <- array(0, c(length(.value), namelen, namelen), 
					list(NULL, namevec, namevec)), 
					list(namelen = length(namevec), namevec = namevec))
  m <- length(subexprs)
  for (i in seq_len(n))
    subexprs[[m+i]] <- substitute(.grad[, name] <- expr,
			          list(name = namevec[i], expr = final[[namevec[i]]]))
  m <- length(subexprs)
  h <- 0
  if (hessian) {
    for (i in seq_len(n))
      for (j in seq_len(n-i+1) + i-1) {
		    h <- h + 1
		    if (i == j)
		      subexprs[[m + h]] <- substitute(.hessian[, i, i] <- expr,
				    list(i = namevec[i], expr = final[[1 + n + h]]))
		    else
		      subexprs[[m + h]] <- substitute(.hessian[, i, j] <- .hessian[, j, i] <- expr,
				    list(i = namevec[i], j = namevec[j], expr = final[[1 + n + h]]))
	    }	
	  h <- h + 1
	  subexprs[[m + h]] <- quote(attr(.value, "hessian") <- .hessian)
  }
  m <- length(subexprs)
  subexprs[[m+1]] <- quote(attr(.value, "gradient") <- .grad)
  subexprs[[m+2]] <- quote(.value)
  subexprs
}

#' @export
#' 
fnDeriv <- function(expr, namevec, args = all.vars(expr), env = environment(expr),
                    do_substitute = FALSE, verbose = FALSE, ...) {
  fn <- function() NULL
  expr <- dex(expr, do_substitute = do_substitute, verbose = verbose)
  body(fn) <- codeDeriv(expr, namevec, do_substitute = FALSE, verbose = verbose, ...)
  if (is.character(args)) {
    formals <- rep(list(bquote()), length = length(args))
    names(formals) <- args
    args <- formals
  }
  formals(fn) <- args
  if (is.null(env))
    environment(fn) <- parent.frame()
  else
    environment(fn) <- as.environment(env)
  fn
}

## #' @name isFALSE
## #' @export
## Ignore this as it causes trouble. Assume R > 3.5
## if (getRversion() < "3.5.0") {
##   isFALSE <- function(x) identical(FALSE, x)
## } else 
##  isFALSE <- isFALSE

#' isZERO
#' 
#' Test if argument is zero
#' 
#' @param x object to be tested
#' 
#' @export
isZERO <- function(x) is.numeric(x) && length(x) == 1 && x == 0

#' isONE
#' 
#' Test if argument is one
#' 
#' @param x object to be tested
#' 
#' @export
isONE  <- function(x) is.numeric(x) && length(x) == 1 && x == 1

#' isMINUSONE
#' 
#' @param x object to be tested
#' 
#' @export
isMINUSONE <- function(x) is.numeric(x) && length(x) == 1 && x == -1

#' isCALL
#'  
#' Test if argument is a call
#' 
#' @param x object to be tested
#' @param name ??need to document better -- is it a character string?
#' 
#' @export
isCALL <- function(x, name) is.call(x) && as.character(x[[1]]) == name


#' nlsSimplify
#' 
#' Try to simplify an expression re: nonlinear least squares
#' 
#' @usage  nlsSimplify(expr, simpEnv = sysSimplifications, verbose = FALSE)
#'
#' @param expr  An expression represented in a variety of ways. See Details.
#'  
#' @param simpEnv  The environment in which simplifications are stored.
#'  
#' @param verbose  If \code{TRUE}, then diagnostic output will be printed as derivatives
#'      and simplifications are recognized.
#'      
#' @export
#' @importFrom digest digest      
nlsSimplify <- function(expr, simpEnv = sysSimplifications, verbose = FALSE) {
    
    if (is.expression(expr))
    	return(as.expression(lapply(expr, nlsSimplify, simpEnv, verbose = verbose)))    

    if (is.call(expr)) {    	
	for (i in seq_along(expr)[-1])
	    expr[[i]] <- nlsSimplify(expr[[i]], simpEnv, verbose = verbose)
	fn <- as.character(expr[[1]])
	nargs <- length(expr) - 1
	while (!identical(simpEnv, emptyenv())) {
	    simps <- simpEnv[[fn]]
	    if (nargs > length(simps))
		return(expr)
	    simpn <- simps[[nargs]]
	    for (i in seq_along(simpn)) {  	
		argnames <- simpn[[i]]$argnames
		substitutions <- lapply(seq_along(argnames)+1L, function(i) expr[[i]])
		names(substitutions) <- argnames
		test <- simpn[[i]]$test
		if (with(substitutions, eval(test))) {
		    if (verbose) {
	    		message(paste("Simplifying", deparse(expr)))		    	
	    		message("Applying simplification:")
		    	print(simpn[[i]])
		    }
		    simplification <- do.call(substitute, list(simpn[[i]]$simplification, substitutions))
	   	    if (simpn[[i]]$do_eval)
	    		simplification <- eval(simplification)
	    	    return(simplification)
	        }
	     }
	     simpEnv <- parent.env(simpEnv)
	}
    }
    expr
}

#' findSubexprs
#' 
#' Try to find the sub-expressions in \code{expr} ??
#' 
#' @usage findSubexprs(expr, simplify = FALSE, tag = ".expr", verbose = FALSE, ...) 
#'
#' @param expr  An expression represented in a variety of ways. See Details.
#' @param simplify  The environment in which simplifications are stored.
#' @param tag   to be attached to the returned object(s)??
#' @param verbose  If \code{TRUE}, then diagnostic output will be printed as derivatives
#'      and simplifications are recognized.
#' @param ... Additional arguments
#'      
#' @export
findSubexprs <- function(expr, simplify = FALSE, tag = ".expr", verbose = FALSE, ...) {
    digests <- new.env(parent = emptyenv())
    subexprs <- list()
    subcount <- 0
    
    record <- function(index) {
        if (simplify)
	         expr[[index]] <<- subexpr <- nlsSimplify(expr[[index]], 
	                     verbose = verbose, ...)
      	else
     	    subexpr <- expr[[index]]
        	if (is.call(subexpr)) {
	           digest <- digest(subexpr)
	           for (i in seq_along(subexpr))
		         record(c(index,i))
	           prev <- digests[[digest]]
	           if (is.null(prev)) assign(digest, index, envir = digests)
        	   else if (is.numeric(prev))  { # the index where we last saw this
      	        subcount <<- subcount + 1
		            name <- as.name(paste0(tag, subcount))
		            assign(digest, name, envir = digests)
        	   }
	        }
     }
    
     edit <- function(index) {
	     subexpr <- expr[[index]]
	     if (is.call(subexpr)) {
	       digest <- digest(subexpr)	    
	       for (i in seq_along(subexpr)) edit(c(index,i))
    	   prev <- digests[[digest]]
	       if (is.name(prev)) {
		       num <- as.integer(substring(as.character(prev), nchar(tag)+1L))
		       subexprs[[num]] <<- call("<-", prev, expr[[index]])
		       expr[[index]] <<- prev
	       } 
	     }
     }
     for (i in seq_along(expr)) record(i)
     for (i in seq_along(expr)) edit(i)
     result <- quote({})
     result[seq_along(subexprs)+1] <- subexprs
     result[[length(result)+1]] <- expr
     result
}
    
# These are the derivatives supported by deriv()

newDeriv(log(x, base = exp(1)), 
         D(x)/(x*log(base)) - (1/base)*log(x, base)/log(base)*D(base))
newDeriv(exp(x), exp(x)*D(x))
newDeriv(sin(x), cos(x)*D(x))
newDeriv(cos(x), -sin(x)*D(x))
newDeriv(tan(x), 1/cos(x)^2*D(x))
newDeriv(sinh(x), cosh(x)*D(x))
newDeriv(cosh(x), sinh(x)*D(x))
newDeriv(sqrt(x), D(x)/2/sqrt(x))
newDeriv(pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE),
  (if (lower.tail && !log.p) 
  	dnorm((q-mean)/sd)*(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2) 
  else if (!lower.tail && !log.p) 
  	-dnorm((q-mean)/sd)*(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2) 
  else if (lower.tail && log.p) 
  	(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2)*exp(dnorm((q-mean)/sd, log = TRUE) 
  		- pnorm((q-mean)/sd, log = TRUE))
  else if (!lower.tail && log.p)	
  	-(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2)*exp(dnorm((q-mean)/sd, log = TRUE) 
  		- pnorm((q-mean)/sd, lower.tail = FALSE, log = TRUE)))
  + stop("cannot take derivative wrt 'lower.tail' or 'log.p'")*(D(lower.tail) + D(log.p)))
newDeriv(dnorm(x, mean = 0, sd = 1, log = FALSE),
  (if (!log) 
  	dnorm((x-mean)/sd)/sd^2*((D(mean)-D(x))*(x-mean)/sd + D(sd)*((x-mean)^2/sd^2 - 1)) 
  else if (log)
        ((D(mean)-D(x))*(x-mean)/sd + D(sd)*((x-mean)^2/sd^2 - 1))/sd) 
  + stop("cannot take derivative wrt 'log'")*D(log))
newDeriv(asin(x), D(x)/sqrt(1+x^2))
newDeriv(acos(x), -D(x)/sqrt(1+x^2))
newDeriv(atan(x), D(x)/(1+x^2))
newDeriv(gamma(x), gamma(x)*digamma(x)*D(x))
newDeriv(lgamma(x), digamma(x)*D(x))
newDeriv(digamma(x), trigamma(x)*D(x))
newDeriv(trigamma(x), psigamma(x, 2L)*D(x))
newDeriv(psigamma(x, deriv = 0L), 
          psigamma(x, deriv + 1L)*D(x) 
        + stop("cannot take derivative wrt 'deriv'")*D(deriv))

newDeriv(x*y, x*D(y) + D(x)*y)
newDeriv(x/y, D(x)/y - x*D(y)/y^2)
newDeriv(x^y, y*x^(y-1)*D(x) + x^y*log(x)*D(y))
newDeriv((x), D(x))
# Need to be careful with unary + or -
newDeriv(`+`(x, y = .MissingVal), if (missing(y)) D(x) else D(x) + D(y))
newDeriv(`-`(x, y = .MissingVal), if (missing(y)) -D(x) else D(x) - D(y))

# These are new

newDeriv(abs(x), sign(x)*D(x))
newDeriv(sign(x), 0)

newDeriv(`~`(x, y = .MissingVal), if (missing(y)) D(x) else D(y))

# Now, the simplifications

newSimplification(+a, TRUE, a)
newSimplification(-a, is.numeric(a), -a, do_eval = TRUE)
newSimplification(-a, isCALL(a, "-") && length(a) == 2, quote(a)[[2]], do_eval = TRUE)

newSimplification(exp(a), isCALL(a, "log") && length(a) == 2, quote(a)[[2]], do_eval = TRUE)
newSimplification(exp(a), is.numeric(a), exp(a), do_eval = TRUE)

newSimplification(log(a), isCALL(a, "exp") && length(a) == 2, quote(a)[[2]], do_eval = TRUE)
newSimplification(log(a), is.numeric(a), log(a), do_eval = TRUE)

newSimplification(!a, isTRUE(a), FALSE)
newSimplification(!a, isFALSE(a), TRUE)

newSimplification((a), TRUE, a)

newSimplification(a + b, isZERO(b), a)
newSimplification(a + b, isZERO(a), b)
newSimplification(a + b, identical(a, b), nlsSimplify(quote(2*a)), do_eval = TRUE)
newSimplification(a + b, is.numeric(a) && is.numeric(b), a+b, do_eval = TRUE)
# Add these to support our error scheme, don't test for stop() everywhere.
newSimplification(a + b, isCALL(a, "stop"), a)
newSimplification(a + b, isCALL(b, "stop"), b)

newSimplification(a - b, isZERO(b), a)
newSimplification(a - b, isZERO(a), nlsSimplify(quote(-b)), do_eval = TRUE)
newSimplification(a - b, identical(a, b), 0)
newSimplification(a - b, is.numeric(a) && is.numeric(b), a - b, do_eval = TRUE)

newSimplification(a * b, isZERO(a), 0)
newSimplification(a * b, isZERO(b), 0)
newSimplification(a * b, isONE(a), b)
newSimplification(a * b, isONE(b), a)
newSimplification(a * b, isMINUSONE(a), nlsSimplify(quote(-b)), do_eval = TRUE)
newSimplification(a * b, isMINUSONE(b), nlsSimplify(quote(-a)), do_eval = TRUE)
newSimplification(a * b, is.numeric(a) && is.numeric(b), a * b, do_eval = TRUE)

newSimplification(a / b, isONE(b), a)
newSimplification(a / b, isMINUSONE(b), nlsSimplify(quote(-a)), do_eval = TRUE)
newSimplification(a / b, isZERO(a), 0)
newSimplification(a / b, is.numeric(a) && is.numeric(b), a / b, do_eval = TRUE)

newSimplification(a ^ b, isONE(b), a)
newSimplification(a ^ b, is.numeric(a) && is.numeric(b), a ^ b, do_eval = TRUE)

newSimplification(log(a, base), isCALL(a, "exp"), nlsSimplify(call("/", quote(a)[[2]], quote(log(base)))), do_eval = TRUE)

newSimplification(a && b, isFALSE(a) || isFALSE(b), FALSE)
newSimplification(a && b, isTRUE(a), b)
newSimplification(a && b, isTRUE(b), a)

newSimplification(a || b, isTRUE(a) || isTRUE(b), TRUE)
newSimplification(a || b, isFALSE(a), b)
newSimplification(a || b, isFALSE(b), a)

newSimplification(if (cond) a, isTRUE(cond), a)
newSimplification(if (cond) a, isFALSE(cond), NULL)

newSimplification(if (cond) a else b, isTRUE(cond), a)
newSimplification(if (cond) a else b, isFALSE(cond), b)
newSimplification(if (cond) a else b, identical(a, b), a)

# This one is used to fix up the unary -
newSimplification(missing(a), identical(a, quote(.MissingVal)), TRUE)
newSimplification(missing(a), !identical(a, quote(.MissingVal)), FALSE)
