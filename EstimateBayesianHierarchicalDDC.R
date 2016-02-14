#=======================================================================
# This program estimates the Dynamic Programming model of the Rust Bus Problem
# using the Bayesian Markov Chain Monte Carlo Method. 

# Original Code from Imai, S., Jain, N., & Ching, A. (2009).
# Bayesian estimation of dynamic discrete choice models. Econometrica, 77(6), 1865-1899.

# Conversion to R by Wayne Taylor

# Version 2/7/2016

#=======================================================================

set.seed(1)

library(bayesm)
data = read.csv("bus_df_RE.csv")
load("simdata.Rdata")

#The order of thetas matches that of the simulated data
#theta[1]: replacement cost
#theta[2]: maintenance cost
theta1true   = 20    # replacement cost
theta2true   = .5    # maintence cost
shocks      = TRUE   # will we be averaging over nsim1 shock simulations or not?
fastConv    = TRUE   # when TRUE, starting thetas are set close to the true thetas
kernelavg   = FALSE   # when TRUE, use kernel density to average over past Ns values

nlgt = max(data$id) # number of cross-sectional units
nvar = 2            # number of variables
nz = ncol(Z)        # number of regressors in mixing distribution
R = 100             # MCMC draws
S = 70              # States
beta0  = .75        # discount rate (set to 0 to ignore forward-looking behavior)
NsAll = 10          # Number of past observations that the kernel values are calculated over
LAll = 5            # Of NsAll, this is how many are kept. These are the LAll "closest" parameters
icheck=rep(0,R)     # Store "check" indicators for the kernal values
kern = rep(0,R)     # Store the kernel values
nsim1 = 10          # Shock simuations 

#Storage
emax  = matrix(0,S,nlgt)      # Emax values for each state x cross-sectional unit
value = array(0,c(S,R,nlgt))  # Keep the updated emax values from each MCMC draw (used for kernel weighting)
sh    = matrix(NA,nsim1,2)    # Shocks

Deltadraw = matrix(0,R,nvar*nz)
thetadraw  = array(0,c(nlgt,nvar,R))
Vthetadraw = matrix(0,R,nvar*nvar)

reject = rep(0,R)
llike = rep(0,R)

#Initialize variables
oldthetas=matrix(0,nlgt,nvar)

#Priors
nu=nvar+3
V=nu*diag(nvar)
Deltabar=matrix(rep(0,nz*nvar),ncol=nvar)
ADelta=.01*diag(nz)

oldVtheta = diag(nvar)
oldVthetai= diag(nvar)
oldDelta = matrix(0,nz,nvar)

stheta = .2 #scaling coefficient on Vtheta. Represents how much smaller an individual's RW step should be relative to overall heterogeneity

#Kernel parameters
weight = 1.0
hkern1 = 1
hkern2 = .02
hkerns1 = hkern1*weight
hkerns2 = hkern2*weight

if(fastConv){
  scale = 1   #how far off should the starting values be from the truth?
  oldthetas = thetas*scale  #'thetas' from DGP file
}

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

#######################################################################################################
#START OF GIBBS SAMPLER
#######################################################################################################
r=1
while(r <= R){
  
  cat('r = ', r,fill=TRUE)
  
  rej = 0
  logl = 0  
  if(fastConv){
    oldVtheta = t(chol(Vtheta))
    oldVthetai=chol2inv(chol(Vtheta))
  }
  sV = stheta*oldVtheta
  root=t(chol(sV))
  
  #	Draw B_i|B-bar, V
  for (i in 1:nlgt) {
    
    datai = subset(data,id==i)
    state = datai$state+1
    y = datai$choice
    
    #Draw new thetas
    thetao = oldthetas[i,]
    thetan = thetao + root%*%rnorm(nvar)
    
    #--------------------------------------------------------------
    # Derive expected values
    #--------------------------------------------------------------
    
    futil_maint = ST_mat%*%emax[,i]
    futil_repl = R_mat%*%emax[,i]
    
    #--------------------------------------------------------------
    # Construct the likelihood
    #--------------------------------------------------------------
    
    #USING 'OLD' PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    val1 = thetao[2]*state+beta0*futil_maint[state]
    val2 = thetao[1]+beta0*futil_repl[state]
    diff = -(val1-val2) #Remember, since these are costs, higher values are worse (hence the negative)
    
    pval1 = (exp(diff)/(1+exp(diff)))
    prob = pval1*(1-y) + (1-pval1)*y
  
    pold = sum(log(prob))
    
    #USING 'NEW' PARAMETERS ("st") ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    val1 = thetan[2]*state+beta0*futil_maint[state]
    val2 = thetan[1]+beta0*futil_repl[state]
    diff = -(val1-val2)
    
    pval1 = (exp(diff)/(1+exp(diff)))
    prob = pval1*(1-y) + (1-pval1)*y
  
    pnew = sum(log(prob))
    
    #--------------------------------------------------------------------------
    # Calculate the next iteration parameters of the Metropolis- Hastings algorithm. 
    #--------------------------------------------------------------------------
    
    # heterogeneity
    logknew = -.5*(t(thetan)-Z[i,]%*%oldDelta) %*% oldVthetai %*% (thetan-t(Z[i,]%*%oldDelta))
    logkold = -.5*(t(thetao)-Z[i,]%*%oldDelta) %*% oldVthetai %*% (thetao-t(Z[i,]%*%oldDelta))
    
    # M-H step
    alpha = exp(pnew + logknew - pold - logkold)
    if(alpha=="NaN") alpha=-1
    u = runif(n=1,min=0, max=1)
    if(u < alpha) { 
      oldthetas[i,] = thetan
      logl = logl + pnew
    } else {
      logl = logl + pold
      rej = rej+1
    }
  }
  
  #	Draw B-bar and V as a multivariate regression
  out=rmultireg(oldthetas,Z,Deltabar,ADelta,nu,V)
  oldDelta=out$B
  oldVtheta=out$Sigma
  oldVthetai=chol2inv(chol(oldVtheta))
  
  Deltadraw[r,]=as.vector(oldDelta)
  Vthetadraw[r,]=as.vector(oldVtheta)
  thetadraw[,,r]=oldthetas
  llike[r]=logl
  reject[r]=rej/nlgt
  
  #================================================================
  # This part derives the Emax functions using the simulated data.
  #================================================================
  
  #Loop through each person
  for(i in 1:nlgt){
    
    #----------------------------------------------------------------
    # Derive the transition probabilities.
    #----------------------------------------------------------------
    
    futil_maint = ST_mat%*%emax[,i]
    futil_repl = R_mat%*%emax[,i]
    
    #Are you averaging over shocks?
    if(shocks){
      
      nhalf = .5*nsim1
      sh[1:nhalf,] = rlogis(nhalf*2)
      sh[(nhalf+1):nsim1,] = -sh[1:nhalf,]
      
      for(s in 1:S){        
        val1s =   thetadraw[i,2,r]*s+beta0*futil_maint[s]+sh[,1] #COST of MAINTAINING
        val2s =   thetadraw[i,1,r] + beta0*futil_repl[s]+sh[,2]  #COST of REPLACING
        
        emax[s,i] = sum(apply(cbind(val1s,val2s),1,min))/nsim1
      }
      
    } else {
       
      val1 =   thetadraw[i,2,r]*(1:S)+beta0*futil_maint  #COST of MAINTAINING
      val2 =   thetadraw[i,1,r]+beta0*futil_repl         #COST of REPLACING
      
      emax[,i] = apply(cbind(val1,val2),1,min)
    }
    
    #Store the emax values for each draw
    value[,r,i] = emax[,i]
  }

  if(r>1 & kernelavg){
    
    #---------------------------------------------------------
    # Using the kernel function, update the expected value matrix
    # --------------------------------------------------------
    
    #We check the past Ns draws (including the current one)
    Ns = min(r-1,NsAll) #this is the most we can check back over, which might be less than NsAll
    L =  min(r-1,LAll)  #this is the most we can average back over, which might be less than LAll
    
    #These are the indices we will check over
    indStart = r-Ns
    indEnd   = r-1
    ind = indStart:indEnd
    
    #----------------------------------
    # Derive the Gaussian Kernel values
    #----------------------------------
    const = (2*pi)^(-L/2) #not really needed, since only relative kernel values matter
    
    for(i in 1:nlgt){  
      #current thetas
      theta1_c   =  thetadraw[i,1,r]
      theta2_c   =  thetadraw[i,2,r]
      
      #past thetas
      theta1_p   =  thetadraw[i,1,ind]
      theta2_p   =  thetadraw[i,2,ind]
      
      rkern11 = (1/hkerns1)*exp(-.5*(((theta1_c-theta1_p)/hkerns1)^2))
      rkern22 = (1/hkerns2)*exp(-.5*(((theta2_c-theta2_p)/hkerns2)^2))
      
      kern[ind] = const*(rkern11*rkern22)
      
      icheck[] = 0 #Reset the check vector
      iperm = order(kern[ind],decreasing = TRUE) #order N(s) kernel values. Higher values are closer (preferred)
      iperm = iperm + (indStart - 1)             #shift the rankings up to the current starting point
      icheck[iperm[1:L]] = iperm[1:L]            #only check the L closest data points
      
      #=================================================================
      # Now, recalculate the Emax function. (see equation 7 on page 1883)
      #=================================================================  
      #This is the equation at the top of pg 1888
      kerncheck = kern[icheck]
      kernMult = kerncheck/sum(kerncheck)
      if(sum(kerncheck) == 0) kernMult = rep(1,length(kerncheck))
      
      emax[,i] = matrix(value[,icheck,i,drop=FALSE],nrow=S,ncol=length(kerncheck))%*%kernMult #weight the past values by the kernel density, more weight given to the closer values
    }
  }
  
  #==========================================================
  # Now, generate nsim new data.
  #==========================================================
  
  #(omitted)
  
  r = r+1
}

print(Deltadraw)
print(1-mean(reject))
matplot(Deltadraw,type="l")
oldthetas
colMeans(Deltadraw)
