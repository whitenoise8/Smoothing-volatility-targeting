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
set.seed(2023)

n = 400
mu = -2
rho = 0.99
eta2 = 0.1

# Simulate data without covariates, i.e. beta=0
sim_data = sv_sim(n,mu,rho,eta2,beta=0)

y = sim_data$y
h = sim_data$h

# VB without smoothing
?VBSV
mod = VBSV(y)

#: Process parameters and E(y)
mod$par
mod$mu_q_beta
#: Plot volatility (log=TRUE provide h_t, log=FALSE provides sigma^2_t)
#: Set true_val=NULL to avoid plotting the simulates process
sv_plot(mod,true_val=h,log=T)
sv_plot(mod,true_val=h,log=F)

#: Compare with MCMC stochvol algorithm looking at the accuracy measure
#: which can be seen as the "overlapping density" between MCMC and VB
library(stochvol)
mcmc = svsample(y, designmatrix = "ar0",
                draws = 20000, burnin = 10000, 
                priorphi = c(1,1), quiet = TRUE)

t = 20
#: Accuracy for the log-variance
accVarApp(VBparVec=c(mod$mu_q_h[t+1],mod$sigma2_q_h[t+1]),
          MCsample=mcmc$latent[[1]][,t],
          type='Normal',
          plotDensities=TRUE,
          parTrue=h[t])
#: Accuracy for the variance
accVarApp(VBparVec=c(mod$mu_q_h[t+1],mod$sigma2_q_h[t+1]),
          MCsample=exp(mcmc$latent[[1]][,t]),
          type='Lognormal',
          plotDensities=TRUE,
          parTrue=exp(h[t]))

#: Predict one step ahead
pred_vb = pred_h(10000,mod)
h_pred = pred_vb$h
y_pred = pred_vb$y

mcpred_h = as(predict(mcmc,steps=1)$h[[1]],"vector")
mcpred_y = as(predict(mcmc,steps=1)$y[[1]],"vector")

#: Evaluate the accuracy of the predictive densities
#: Log-scale
accVarApp(VBparVec=y_pred,
          MCsample=mcpred_y,
          plotDensities=TRUE)
#: Variance scale
accVarApp(VBparVec=exp(h_pred),
          MCsample=exp(mcpred_h),
          plotDensities=TRUE)



#: Include smoothing - set the matrix W
#: B-spline
w = 10
knots = seq(w,n-w,by=w)
W = bSpline(0:n, 
            knots=knots, 
            degree=3, 
            intercept=TRUE)
mod = VBSV(y,W=W)

mod$par
sv_plot(mod,true_val=h,log=T)
sv_plot(mod,true_val=h,log=F)


#: Daubechies Wavelets
numLevels = 4
W = WDaub(0:n,numLevels=numLevels)
mod = VBSV(y,W=W)

mod$par
sv_plot(mod,true_val=h,log=T)
sv_plot(mod,true_val=h,log=F)
