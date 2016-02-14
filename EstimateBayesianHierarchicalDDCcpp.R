#=======================================================================
# This program estimates the Dynamic Programming model of the Rust Bus Problem
# using the Bayesian Markov Chain Monte Carlo Method. 

# Based on code from Imai, S., Jain, N., & Ching, A. (2009).
# Bayesian estimation of dynamic discrete choice models. Econometrica, 77(6), 1865-1899.

# Author: Wayne Taylor

# Version 2/14/2016

#=======================================================================

set.seed(1)

library(bayesm)
library(Rcpp)
sourceCpp('bddcMCMCloop.cpp')

data = read.csv("bus_df_RE.csv")
load("simdata.Rdata")

#The order of thetas matches that of the simulated data
#theta[1]: replacement cost
#theta[2]: maintenance cost
theta1true   = 20    # replacement cost
theta2true   = .5    # maintence cost
shocks      = FALSE  # will we be averaging over nsim1 shock simulations or not? 
#NOTE: since this example uses GEV errors, there is a closed form expression so there is no need to average over shocks
#It is included only for illustrative purposes
fastConv    = FALSE  # when TRUE, starting thetas are set close to the true thetas
kernelavg   = TRUE   # when TRUE, use kernel density to average over past Ns values

nlgt = max(data$id) # number of cross-sectional units
nvar = 2            # number of variables
nz = ncol(Z)        # number of regressors in mixing distribution
R = 5000            # MCMC draws
S = 70              # States
beta0  = .75        # discount rate (set to 0 to ignore forward-looking behavior)
NsAll = 10          # Number of past observations that the kernel values are calculated over
LAll = 5            # Of NsAll, this is how many are kept. These are the LAll "closest" parameters
nsim1 = 10          # Shock simuations 

#Priors
nu=nvar+3
V=nu*diag(nvar)
Deltabar=matrix(rep(0,nz*nvar),ncol=nvar)
ADelta=.01*diag(nz)

oldVtheta = diag(nvar)
oldVthetai= diag(nvar)
oldDelta = matrix(0,nz,nvar)

stheta = .03 #scaling coefficient on Vtheta. Represents how much smaller an individual's RW step should be relative to overall heterogeneity

#Kernel parameters
weight = 1.0
hkern1 = 1
hkern2 = .02
hkerns1 = hkern1*weight
hkerns2 = hkern2*weight
hkerns = c(hkerns1,hkerns2)

#Transition matrices (assumed to be known, but these can  be estimated from the data)
ST_mat = matrix(0,S,S)
p=c(.36,.48,.16)
lp = length(p)
for(i in 1:S){
  for(j in 1:lp){
    if((i+j-1)<S)  ST_mat[i,i+j-1] = p[j]
    if((i+j-1)==S) ST_mat[i,i+j-1] = sum(p[j:lp]) #out of columns, so collapse the probabilities
  }
}
R_mat = cbind(1,matrix(0,S,S-1))

#Put Data in List() format
lgtdata = NULL
for(i in 1:nlgt){
  datai = subset(data,id==i)
  state = datai$state+1
  y = datai$choice
  lgtdata[[i]] = list(y=y,state=state)
}

#######################################################################################################
# GIBBS SAMPLER
#######################################################################################################
out = bddcMCMCloop(R,nsim1,States=1:S,lgtdata,Z,ST_mat,R_mat,beta0,
                   nu,V,Deltabar,ADelta,oldVtheta,oldVthetai,oldDelta,
                   stheta,fastConv,shocks,kernelavg,1,thetas,Vtheta,NsAll,LAll,hkerns)

print(out$Deltadraw)
matplot(out$Deltadraw,type="l")
colMeans(out$Deltadraw)
