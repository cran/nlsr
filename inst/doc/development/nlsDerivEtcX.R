library(nlsr)
y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558,
       50.156, 62.948, 75.995, 91.972)
tt <- seq_along(y)  # for testing
mydata <- data.frame(y = y, tt = tt)
frm <- y ~ b1/(1 + b2 * exp(-1 * b3 * tt))
st <- c(b1=2, b2=1, b3=1) # a default starting vector (named!)
nst<-names(st)
rfrm<-as.character(frm[3])
# Try different inputs
try(nlsDeriv(frm[3], "b1")) # fails 
try(nlsDeriv(expression(frm[3]), "b1"))# fails 
try(nlsDeriv(as.expression(frm[3]), "b1"))# fails 
try(nlsDeriv(as.character(frm[3]), "b1")) # works
try(nlsDeriv(as.character(rfrm), "b1")) # works
crfrm <- "b1/(1 + b2 * exp(-1 * b3 * tt))"
identical(crfrm, rfrm)
tnlsD1 <- nlsDeriv(crfrm, "b1")
tnlsD1
tnlsD2 <- nlsDeriv(crfrm, "b2")
tnlsD2
tnlsD3 <- nlsDeriv(crfrm, "b3")
tnlsD3
fnDeriv(crfrm, nst)
fnDeriv(crfrm, nst, verbose=TRUE)
# ?? how to use do_substitute, args, and env
codeDeriv(crfrm, nst)
codeDeriv(crfrm, nst, hessian=TRUE)
codeDeriv(crfrm, nst, verbose=TRUE)
# ?? how to use do_substitute, args, and derivEnv
# ?? Examples of ... in these functions
efrm<-expression(b1/(1 + b2 * exp(-1 * b3 * tt)))
str(efrm)
# numericDeriv and numericDerivR
tt <- 1:12
b1 <- 200; b2 <- 50; b3 <- 0.3
try(numericDerivR(efrm, nst))
# original fails. Why ??
try(numericDeriv(efrm, nst))
