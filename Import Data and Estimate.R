#**********************************************************************                                                         
#  This program estimates John Rust's bus engine replacement model    
#  (Rust, Econometrica 1987) using the Nested Pseudo Likelihood (NPL) 
#  algorithm in Aguirregabiria and Mira (Econometrica, 2002).

#  Original code in GAUSS by Victor Aguirregabiria see http://individual.utoronto.ca/vaguirre/software/library_procedures.html
#  R conversion by Wayne Taylor

#  Version: 12/8/2015

#  The program is intended to be used with data generated from "Rust Data Generating Process.R"
#**********************************************************************

#Support functions 
source('npl_sing.R')
source('clogit.R')

discthre=function(y,thre){
  
  # DISCTHRE - Discretization and Codification of a Variable using a prefixed vector of thresholds.
  
  # Original code in GAUSS by Victor Aguirregabiria
  # Converted to R by Wayne Taylor
  
  # Version: 12/8/2015
  
  # Input       y    - (nobs x 1) vector with observations of the continuous variable
  #             thre - (k x 1) vector of thresholds
  
  # Output      discy - (nobs x 1) vector with the codes of the discretized observations
  #                     Example: If y[i]>thre[5] and y[i]<=thre[6], then discy[i]=6
  
  numcel = length(thre)
  discy = rep(0,length(y))
  discy = discy + 1*(y <= thre[1])
  
  j=2
  while(j <= numcel){
    discy = discy + j*(y>thre[j-1])*(y<=thre[j])
    j=j+1
  }
  
  discy = discy + (numcel+1)*(y>thre[numcel])
  
  discy  
}

kernel1=function(xobs,xpred){
  
  # KERNEL1 - Kernel Estimation of a Univariate Density Function using a Gaussian Kernel 
  # The bandwith is equal to Silverman's rule of thumb
  
  # Original code in GAUSS by Victor Aguirregabiria
  # Converted to R by Wayne Taylor
  
  # Version: 12/8/2015
  
  # Input:       xobs    - (N x 1) vector of observations
  #              xpred   - (K x 1) vector of values where the pdf will be estimated
  
  # Output:      pest    - (K x 1) vector of estimates
  
  nobs  = length(xobs)
  k = length(xpred)
  
  #Silverman's rule of thumb bandwidth
  band0 = (1.364 * sd(xobs))/((pi*nobs)^(1/5))
  
  pest = colSums(dnorm((matrix(xpred,nobs,k,byrow=TRUE)-xobs)/band0))/(nobs*band0)
  
  pest
}

#*********************
# 0. Some constants
#*********************

#Names of structural parameters
namespar = c("ReplaceC","MaintenC")

#Value of discount factor
beta = 0.75

#Number of state variables
kvarx = 1

#Number of choice alternatives
jchoice = 2

#Number of cells in discretization of the state variable
ncelx = 70

#The output from the accompanying DGP file is already discretized
#and the transition matrices are known so some steps below are not
#required. However, to understand the original function set the
#variable to TRUE to discretize x and estimate the transition matrices
orig_est = TRUE

#******************************
# 1. Reading Rust's bus data
#******************************
#Datafile
filedat  = "bus_df_lin.csv"
data = read.csv(filedat)

if(!orig_est){
  aobs = data$choice #Replacement decision  
} else {
  data = data[order(data$id),]
  nobs = nrow(data)  #Cumulative mileage
  aobs = data$choice #Replacement decision  
  iobs = data$id     #Individual bus ID
  xobs = (data$state+runif(nobs))*1e6 #un-discretizing the data (to illustrate how the original code works)
  xobs = xobs/1e6
}

#********************************************************
# 2. Discretization of the decision and state variable
#********************************************************
indobsa = aobs+1 #indaobs should be 1,2,...,J

if(!orig_est){
  xval = 1:ncelx         #cost of maintaining at each state
  indobsx = data$state+1 #index starts at 1
} else{
  minx = quantile(xobs,0.01)
  maxx = quantile(xobs,0.99)
  stepx = round(1e6*(maxx-minx)/(ncelx-1))/1e6
  xthre = seq(minx+stepx,by=stepx,length.out=ncelx-1)
  xval = seq(minx,by=stepx,length.out=ncelx)
  indobsx = discthre(xobs,xthre)
  head(indobsx)
}

#****************************************
# 3. Specification of utility function
#****************************************
zmat1= cbind(0,-xval)            #if maintain (zmat1), replace cost = 0 and maintain cost = -xval
zmat2= cbind(-1,rep(0,ncelx))    #if replace (zmat2), replace cost = 1 and maintain cost = 0
zmat=cbind(zmat1,zmat2)
rm(zmat1,zmat2)

#Or it can be set up this way, the result is the same
# zmat1= cbind(1,-xval)            #Realtive to replace, if maintain (zmat1), replace cost = 1 and maintain cost = -xval
# zmat2= cbind(0,rep(0,ncelx))
# zmat=cbind(zmat1,zmat2)
# rm(zmat1,zmat2)

#****************************************************************
# 4. Estimation of transition probabilities of state variables
#****************************************************************

if(!orig_est){
  fmat1 = matrix(0,ncelx,ncelx)
  p=c(.36,.48,.16)
  lp = length(p)
  for(i in 1:ncelx){
    for(j in 1:lp){
      if((i+j-1)<ncelx)  fmat1[i,i+j-1] = p[j]
      if((i+j-1)==ncelx) fmat1[i,i+j-1] = sum(p[j:lp]) #out of columns, so collapse the probabilities
    }
  }
  fmat2 = cbind(1,matrix(0,ncelx,ncelx-1))
  fmat=cbind(fmat1,fmat2)
} else {
  # Nonparametric PDF of additional mileage
  iobs_1 = c(0,iobs[-nobs])
  xobs_1 = c(0,xobs[-nobs])
  aobs_1 = c(0,aobs[-nobs])
  dxobs = (1-aobs_1)*(xobs-xobs_1) + aobs_1*xobs
  dxobs = dxobs[iobs==iobs_1]
  mindx = 0
  maxdx = quantile(dxobs,0.999)
  numdx = 2 + round(maxdx/stepx)
  dxval = seq(0,by=stepx,length.out=numdx)
  pdfdx = kernel1(dxobs,dxval)
  pdfdx = pdfdx/sum(pdfdx)
  
  # Transition matrices
  fmat2 = kronecker(matrix(1,ncelx,1),t(c(pdfdx,rep(0,ncelx-numdx))))
  
  fmat1 = t(c(pdfdx,rep(0,ncelx-numdx)))
  j=2
  while(j<=(ncelx-1)){
    colz = ncelx - (j-1+numdx)
    if (colz>0) fmat1 = rbind(fmat1,c(rep(0,j-1),pdfdx,rep(0,colz)))
    if (colz==0) fmat1 = rbind(fmat1,c(rep(0,j-1),pdfdx))
    if (colz<0){
      buff = c(pdfdx[1:(numdx+colz-1)],sum(pdfdx[(numdx+colz):numdx]))
      fmat1 = rbind(fmat1,c(rep(0,j-1),buff))
    }
    j=j+1
  }
  fmat1 = rbind(fmat1,c(rep(0,ncelx-1),1))
  
  fmat=cbind(fmat1,fmat2)
}

#***********************************
# 5. Initial choice probabilities
#***********************************
xprob0=cbind(indobsx,indobsx^2,indobsx^3)
glmout=glm(aobs~xprob0-1,family='binomial')
coef=glmout$coef
est = exp(-(cbind(xval,xval^2,xval^3)%*%coef))
prob0 = 1/(1+est)                              #pr(replace) (columnas)|state (rows)
prob0 = cbind(1-prob0,prob0)                   #probability of maintain[,1] or replace[,2]|state 

#***************************
# 6. Structural estimation
# ***************************
out = npl_sing(indobsa,indobsx,zmat,prob0,beta,fmat,namespar)

out[[length(out)]][[1]]

# Bus data true parameters: MC = 20 RC = .5

# orig_est=FALSE
# [1,] 19.2014940
# [2,]  0.5233962

# orig_est=TRUE
# [1,] 17.8076640
# [2,]  0.5922635
