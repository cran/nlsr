## ----derivexp, eval=TRUE------------------------------------------------------
partialderiv <- D(expression(a * (x ^ b)),"b")
print(partialderiv)

## ----hiddenexprob1, echo=FALSE------------------------------------------------
# Here we set up an example problem with data
# Define independent variable
t0 <- 0:19
t0a<-t0
t0a[1]<-1e-20 # very small value
# Drop first value in vectors
t0t<-t0[-1]
y1 <- 4 * (t0^0.25)
y1t<-y1[-1]
n <- length(t0)
fuzz <- rnorm(n)
range <- max(y1)-min(y1)
## add some "error" to the dependent variable
y1q <- y1 + 0.2*range*fuzz
edta <- data.frame(t0=t0, t0a=t0a, y1=y1, y1q=y1q)
edtat <- data.frame(t0t=t0t, y1t=y1t)

start1 <- c(a=1, b=1)
try(nlsy0t0ax <- nls(formula=y1~a*(t0a^b), start=start1, data=edta, control=nls.control(maxiter=10000)))
library(nlsr)
nlsry1t0a <- nlxb(formula=y1~a*(t0a^b), start=start1, data=edta)
library(minpack.lm)
nlsLMy1t0a <- nlsLM(formula=y1~a*(t0a^b), start=start1, data=edta)

## ----nlsstr-------------------------------------------------------------------
str(nlsy0t0ax)

## ----nlsrstr------------------------------------------------------------------
str(nlsry1t0a)

## ----pred1--------------------------------------------------------------------
nudta <- colMeans(edta)
predict(nlsy0t0ax, newdata=nudta)
predict(nlsLMy1t0a, newdata=nudta)
predict(nlsry1t0a, newdata=nudta)

## ----exprob1------------------------------------------------------------------
# Here we set up an example problem with data
# Define independent variable
t0 <- 0:19
t0a<-t0
t0a[1]<-1e-20 # very small value
# Drop first value in vectors
t0t<-t0[-1]
y1 <- 4 * (t0^0.25)
y1t<-y1[-1]
n <- length(t0)
fuzz <- rnorm(n)
range <- max(y1)-min(y1)
## add some "error" to the dependent variable
y1q <- y1 + 0.2*range*fuzz
edta <- data.frame(t0=t0, t0a=t0a, y1=y1, y1q=y1q)
edtat <- data.frame(t0t=t0t, y1t=y1t)

## ----extrynls-----------------------------------------------------------------
cprint <- function(obj){
   # print object if it exists
  sobj<-deparse(substitute(obj))
  if (exists(sobj)) {
      print(obj)
  } else {
      cat(sobj," does not exist\n")
  }
#  return(NULL)
}
start1 <- c(a=1, b=1)
try(nlsy0t0 <- nls(formula=y1~a*(t0^b), start=start1, data=edta))
cprint(nlsy0t0)
# Since this fails to converge, let us increase the maximum iterations
try(nlsy0t0x <- nls(formula=y1~a*(t0^b), start=start1, data=edta,
                    control=nls.control(maxiter=10000)))
cprint(nlsy0t0x)
try(nlsy0t0ax <- nls(formula=y1~a*(t0a^b), start=start1, data=edta, 
                     control=nls.control(maxiter=10000)))
cprint(nlsy0t0ax)
try(nlsy0t0t <- nls(formula=y1t~a*(t0t^b), start=start1, data=edtat))
cprint(nlsy0t0t)

## ----extry1nlsr---------------------------------------------------------------
library(nlsr)
nlsry1t0 <- try(nlxb(formula=y1~a*(t0^b), start=start1, data=edta))
cprint(nlsry1t0)
nlsry1t0a <- nlxb(formula=y1~a*(t0a^b), start=start1, data=edta)
cprint(nlsry1t0a)
nlsry1t0t <- nlxb(formula=y1t~a*(t0t^b), start=start1, data=edtat)
cprint(nlsry1t0t)

## ----extry1minlm--------------------------------------------------------------
library(minpack.lm)
nlsLMy1t0 <- nlsLM(formula=y1~a*(t0^b), start=start1, data=edta)
nlsLMy1t0
nlsLMy1t0a <- nlsLM(formula=y1~a*(t0a^b), start=start1, data=edta)
nlsLMy1t0a
nlsLMy1t0t <- nlsLM(formula=y1t~a*(t0t^b), start=start1, data=edtat)
nlsLMy1t0t

## ----brownden-----------------------------------------------------------------
m <- 20
t <- seq(1, m) / 5
y <- rep(0,m)
library(nlsr)
library(minpack.lm)

bddata <- data.frame(t=t, y=y)
bdform <- y ~ ((x1 + t * x2 - exp(t))^2 + (x3 + x4 * sin(t) - cos(t))^2)
prm0 <- c(x1=25, x2=5, x3=-5, x4=-1)
fbd <-model2ssgrfun(bdform, prm0, bddata)
cat("initial sumsquares=",as.numeric(crossprod(fbd(prm0))),"\n")
nlsrbd <- nlxb(bdform, start=prm0, data=bddata, trace=FALSE)
nlsrbd

nlsbd10k <- nls(bdform, start=prm0, data=bddata, trace=FALSE, 
                control=nls.control(maxiter=10000))
nlsbd10k

nlsLMbd10k <- nlsLM(bdform, start=prm0, data=bddata, trace=FALSE, 
                    control=nls.lm.control(maxiter=10000, maxfev=10000))
nlsLMbd10k

## ----browndenpred-------------------------------------------------------------
ndata <- data.frame(t=c(5,6), y=c(0,0))
predict(nlsLMbd10k, newdata=ndata)

# now nls
predict(nlsbd10k, newdata=ndata)

# now nlsr
predict(nlsrbd, newdata=ndata)



## ----brownden2----------------------------------------------------------------
bdf2 <-  (x1 + t * x2 - exp(t))^2 ~ - (x3 + x4 * sin(t) - cos(t))^2

nlsbd2 <- try(nls(bdf2, start=prm0, data=bddata, trace=FALSE))
nlsbd2

nlsLMbd2 <- try(nlsLM(bdf2, start=prm0, data=bddata, trace=FALSE, 
                      control=nls.lm.control(maxiter=10000, maxfev=10000)))
nlsLMbd2

nlsrbd2 <- nlxb(bdf2, start=prm0, data=bddata, trace=FALSE)
##summary(nlsrbd2)
nlsrbd2

## ----brownden2pred------------------------------------------------------------
ndata <- data.frame(t=c(5,6), y=c(0,0))
predict(nlsrbd2, newdata=ndata)

## ----bdcheck, eval=FALSE------------------------------------------------------
#  #' Brown and Dennis Function
#  #'
#  #' Test function 16 from the More', Garbow and Hillstrom paper.
#  #'
#  #' The objective function is the sum of \code{m} functions, each of \code{n}
#  #' parameters.
#  #'
#  #' \itemize{
#  #'   \item Dimensions: Number of parameters \code{n = 4}, number of summand
#  #'   functions \code{m >= n}.
#  #'   \item Minima: \code{f = 85822.2} if \code{m = 20}.
#  #' }
#  #'
#  #' @param m Number of summand functions in the objective function. Should be
#  #'   equal to or greater than 4.
#  #' @return A list containing:
#  #' \itemize{
#  #'   \item \code{fn} Objective function which calculates the value given input
#  #'   parameter vector.
#  #'   \item \code{gr} Gradient function which calculates the gradient vector
#  #'   given input parameter vector.
#  #'   \item \code{fg} A function which, given the parameter vector, calculates
#  #'   both the objective value and gradient, returning a list with members
#  #'   \code{fn} and \code{gr}, respectively.
#  #'   \item \code{x0} Standard starting point.
#  #' }
#  #' @references
#  #' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
#  #' Testing unconstrained optimization software.
#  #' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
#  #' \url{https://doi.org/10.1145/355934.355936}
#  #'
#  #' Brown, K. M., & Dennis, J. E. (1971).
#  #' \emph{New computational algorithms for minimizing a sum of squares of
#  #' nonlinear functions} (Report No. 71-6).
#  #' New Haven, CT: Department of Computer Science, Yale University.
#  #'
#  #' @examples
#  #' # Use 10 summand functions
#  #' fun <- brown_den(m = 10)
#  #' # Optimize using the standard starting point
#  #' x0 <- fun$x0
#  #' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
#  #' "L-BFGS-B")
#  #' # Use your own starting point
#  #' res <- stats::optim(c(0.1, 0.2, 0.3, 0.4), fun$fn, fun$gr, method =
#  #' "L-BFGS-B")
#  #'
#  #' # Use 20 summand functions
#  #' fun20 <- brown_den(m = 20)
#  #' res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
#  #' @export
#  #`
#  brown_den <- function(m = 20) {
#    list(
#      fn = function(par) {
#        x1 <- par[1]
#        x2 <- par[2]
#        x3 <- par[3]
#        x4 <- par[4]
#  
#        ti <- (1:m) * 0.2
#        l <- x1 + ti * x2 - exp(ti)
#        r <- x3 + x4 * sin(ti) - cos(ti)
#        f <- l * l + r * r
#        sum(f * f)
#      },
#      gr = function(par) {
#        x1 <- par[1]
#        x2 <- par[2]
#        x3 <- par[3]
#        x4 <- par[4]
#  
#        ti <- (1:m) * 0.2
#        sinti <- sin(ti)
#        l <- x1 + ti * x2 - exp(ti)
#        r <- x3 + x4 * sinti - cos(ti)
#        f <- l * l + r * r
#        lf4 <- 4 * l * f
#        rf4 <- 4 * r * f
#        c(
#          sum(lf4),
#          sum(lf4 * ti),
#          sum(rf4),
#          sum(rf4 * sinti)
#        )
#      },
#      fg = function(par) {
#        x1 <- par[1]
#        x2 <- par[2]
#        x3 <- par[3]
#        x4 <- par[4]
#  
#        ti <- (1:m) * 0.2
#        sinti <- sin(ti)
#        l <- x1 + ti * x2 - exp(ti)
#        r <- x3 + x4 * sinti - cos(ti)
#        f <- l * l + r * r
#        lf4 <- 4 * l * f
#        rf4 <- 4 * r * f
#  
#        fsum <- sum(f * f)
#        grad <- c(
#          sum(lf4),
#          sum(lf4 * ti),
#          sum(rf4),
#          sum(rf4 * sinti)
#        )
#  
#        list(
#          fn = fsum,
#          gr = grad
#        )
#      },
#      x0 = c(25, 5, -5, 1)
#    )
#  }
#  mbd <- brown_den(m=20)
#  mbd
#  mbd$fg(mbd$x0)
#  bdsolnm <- optim(mbd$x0, mbd$fn, control=list(trace=0))
#  bdsolnm
#  bdsolbfgs <- optim(mbd$x0, mbd$fn, method="BFGS", control=list(trace=0))
#  bdsolbfgs
#  
#  library(optimx)
#  methlist <- c("Nelder-Mead","BFGS","Rvmmin","L-BFGS-B","Rcgmin","ucminf")
#  
#  solo <- opm(mbd$x0, mbd$fn, mbd$gr, method=methlist, control=list(trace=0))
#  summary(solo, order=value)
#  
#  ## A failure above is generally because a package in the 'methlist' is not installed.
#  

