#: Last update 13 december 2022
#: Author: Nicolas Bianco (nicolas.bianco@phd.unipd.it)
#: [for any issue, feel free to mail me]

#: Example code for:
#: Smoothing volatility targeting (Bernardi and Bianchi and Bianco, 2022)

#: Required packages: 
#: Rcpp, RcppArmadillo, RcppEigen, RcppNumerical, BH,
#: ggplot2, TeachingDemos, wavethresh
#: Download dependencies
list_of_packages = c("Rcpp", "RcppArmadillo", "RcppEigen", "RcppNumerical", "BH", 
                     "ggplot2", "TeachingDemos", "wavethresh")
new_packages = list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)>0) install.packages(new_packages)

#: Load package
library(VBSV)

# Libraries
library(splines2)

# Example
set.seed(2022)

n = 500
mu = -2
rho = 0.99
eta2 = 0.1

# Simulate data without covariates, i.e. beta=NULL
sim_data = sv_sim(n,mu,rho,eta2,beta=NULL)

y = sim_data$y
h = sim_data$h

# VB without smoothin (i.e. smooth=FALSE)
# To approximate with the homoskedastic GMRF set homo=TRUE
hyper = list(A=0.1,B=0.1,s2=100,s2beta=100)
mod = VBSV(y,X=NULL,W=NULL,hyper,options=list(homo=FALSE,smooth=FALSE),Trace=1)

#: Process parameters
mod$par
#: Plot volatility (log=TRUE provide h_t, log=FALSE provides sigma^2_t)
#: Set true_val=NULL to avoid plotting the simulates process
sv_plot(mod,true_val=h,log=T)
sv_plot(mod,true_val=h,log=F)

#: Predict one step ahead
hpred = pred_h(1000,mod)
plot(density(hpred),main="predictive density",xlab=expression(h[t+1]))
plot(density(exp(hpred)),main="predictive density",xlab=expression(sigma[t+1]^2))

#: Include smoothing - set the matrix W
#: B-spline
w = 10
knots = seq(w,n-w,by=w)
W = bSpline(0:n, 
            knots=knots, 
            degree=3, 
            intercept=TRUE)
mod = VBSV(y,X=NULL,W=W,hyper,options=list(homo=FALSE,smooth=TRUE),Trace=1)

mod$par
sv_plot(mod,true_val=h,log=T)
sv_plot(mod,true_val=h,log=F)


#: Daubechies Wavelets
numLevels = 4
k = 2^numLevels - 1
W = cbind(rep(1,n+1),ZDaub(0:n,numLevels=numLevels))
mod = VBSV(y,X=NULL,W=W,hyper,options=list(homo=FALSE,smooth=TRUE),Trace=1)

mod$par
sv_plot(mod,true_val=h,log=T)
sv_plot(mod,true_val=h,log=F)
