## Problem in 1 parameter to ensure methods work in trivial case
require(nlsr)
cat("Check nlsr.control\n")
cat("defaults\n")
ctrl <- nlsr.control()
print(ctrl)
print(names(ctrl))
ctrl<-NULL
cat("Try control<-list(femax=4, phi=.12)\n")
control<-list(femax=4, phi=.12)
ctrl <- nlsr.control(control)
print(ctrl)
print(names(ctrl))
ctrl<-NULL
cat("Try control<-list(fellmax=4, phi=.12) [bad name]\n")
control<-list(fellmax=4, phi=.12)
ctrl <- try(nlsr.control(control))
print(ctrl)
print(names(ctrl))
