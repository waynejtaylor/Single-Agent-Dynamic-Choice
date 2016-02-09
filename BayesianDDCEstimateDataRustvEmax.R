#=======================================================================
# This program estimates the Dynamic Programming model of the Rust Bus Problem
# using the Bayesian Markov Chain Monte Carlo Method. 

# Original Code from Imai, S., Jain, N., & Ching, A. (2009).
# Bayesian estimation of dynamic discrete choice models. Econometrica, 77(6), 1865-1899.

# This is based on their base program (no random effects, non-continuous states)

# Conversion to R by Wayne Taylor

# Version 2/3/2016

# Note: version "Emax" stores only the maximum value for each state, rather than the expected value for each state under each action
#=======================================================================

#=======================================================================
#IMPORTANT NOTES
#This is the first major conversion step from the original code
#It uses data generated from 'Rust Data Generating Process.R'
#To keep it simple, I am doing the following:
#-leave out random effects for now
#-only include one state (but this should be easy to generalize, since more states can be added in matrix form)

#I focus on getting the structure correct and clarifying the naming conventions
#=======================================================================
iseed = 3
set.seed(iseed)

#Based on basyesfv9.R
filedat  = "bus_df_lin.csv"
data = read.csv(filedat)
data = subset(data,id<=100) #For speed improvements, only include 100 buses

state    = data$state+1
y        = data$choice
mileage  = data$mileage

beta1true   = .5    #maintence cost
beta2true   = 20    #replacement cost
beta0true   = .75   #discount rate
shocks      = TRUE #will be averaging over nsim1 shock simulations or not?

ndata = nrow(data)
npar = 2
ngibb = 3000        # MCMC draws
accept = 0          # Counter ofof M-H acceptances
S = 70              # States
NsAll = 100         # Number of past observations that the kernel values are calculated over
LAll = 50           # Of NsAll, this is how many are kept. These are the LAll "closest" parameters
icheck=rep(0,ngibb) # Store "check" indicators for the kernal values
kern = rep(0,ngibb) # Store the kernel values
nsim1 = 10          # Shock simuations 

#Storage
emax  = rep(0,S)                 # Stores the Emax values
value = matrix(0,S,ngibb)        # At each iteration, keep the updated emax values
sh = matrix(NA,nsim1,2)          # Shocks
param = matrix(NA,ngibb,npar*2)  # Parameter coefficients and their draw M-H variances are stored here

#Output for select parameters
bayes1.out = data.frame(igibb=0,beta1=0,beta2=0)

#Kernel parameters
weight  = 1.0
hkern1 = .02
hkern2 = 1
hkerns1 = hkern1*weight
hkerns2 = hkern2*weight

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

scale = .7   #how far off should the starting values be from the truth?
#Starting values
beta1st   = beta1true*scale
beta2st   = beta2true*scale
beta0st   = beta0true
sigbeta1st= .01
sigbeta2st= .1

#'Old' values
beta1     = beta1st
beta2     = beta2st
beta0     = beta0st
sigbeta1  = sigbeta1st
sigbeta2  = sigbeta2st

#Initialize first row
param[1,1] = beta1st
param[1,2] = beta2st
param[1,3] = sigbeta1st
param[1,4] = sigbeta2st

#######################################################################################################
#START OF GIBBS SAMPLER
#######################################################################################################
igibb=1
while(igibb <= ngibb){
  
  tic = proc.time()[3]
  timeelap = tic
  
  if(igibb %% 100 == 0) cat('igibb = ', igibb,fill=TRUE)
  
  #Initiate variables
  beta1st    = param[igibb,1]
  beta2st    = param[igibb,2]
  sigbeta1st = param[igibb,3]
  sigbeta2st = param[igibb,4]
  
  #--------------------------------------------------------------
  # Derive expected values
  #--------------------------------------------------------------
  
  futil_maint = ST_mat%*%emax
  futil_repl = R_mat%*%emax

  #--------------------------------------------------------------
  # Construct the likelihood
  #--------------------------------------------------------------
  
  #USING 'OLD' PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  val1 = beta1*state+beta0*futil_maint[state]
  val2 = beta2+beta0*futil_repl[state]
  
  diff = -(val1-val2) #Remember, since these are costs, higher values are worse (hence the negative)
  pval1 = (exp(diff)/(1+exp(diff)))
  
  #This way is more clear, same answer   
  #cost_array = cbind(val1,val2)
  #cost = cost_array-apply(cost_array,1,min)
  #util = exp(-cost)
  #pchoice = util/rowSums(util) #maintain/replace
  #all.equal(pval1,pchoice[,1])
  
  #choice 0 = maintain
  #choice 1 = replace
  prob = pval1*(1-y) + (1-pval1)*y
    
  pold = sum(log(prob))
    
  #USING 'NEW' PARAMETERS ("st") ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  val1 = beta1st*state+beta0*futil_maint[state]
  val2 = beta2st+beta0*futil_repl[state]
  
  diff = -(val1-val2)
  pval1 = (exp(diff)/(1+exp(diff)))
  
  prob = pval1*(1-y) + (1-pval1)*y
    
  pnew = sum(log(prob))

  #--------------------------------------------------------------------------
  # Calculate the next iteration parameters of the Metropolis- Hastings algorithm. 
  #--------------------------------------------------------------------------
  
  phast = exp(pnew-pold)
  if(phast=="NaN") phast=-1
  hast = runif(1)
  if(hast < phast) {
    beta1    = beta1st
    beta2    = beta2st
    accept = accept+1
  }
  
  #--------------------------------------------------------------------------
  # Derive the candidate distribution for Metropolis Hastings algorithm. 
  #--------------------------------------------------------------------------
  beta1st = beta1+sigbeta1*rnorm(1)
  beta2st = beta2+sigbeta2*rnorm(1)    

  #--------------------------------------------------------
  # Update the variables and record the output
  #--------------------------------------------------------
  if(igibb<ngibb){
    param[igibb+1,1] = beta1st
    param[igibb+1,2] = beta2st
    param[igibb+1,3] = sigbeta1st
    param[igibb+1,4] = sigbeta2st
    
    bayes1.out[igibb,] = c(igibb,beta1,beta2)
  }
  
  toc = proc.time()[3]
  time1 = toc-tic
  ddtime = (toc - timeelap)/60
  
  #================================================================
  # This part derives the Emax functions using the simulated data.
  #================================================================
  
  #----------------------------------------------------------------
  # Derive the transition probabilities.
  #----------------------------------------------------------------

  #(using fuitl_maint and futil_repl from above)
  #In the original code these are recalculated, probably to perturb over the state space some more
  
  #Are you averaging over shocks?
  if(shocks){
    
    nhalf = .5*nsim1
    sh[1:nhalf,] = rlogis(nhalf*2)
    sh[(nhalf+1):nsim1,] = -sh[1:nhalf,]
    
    for(s in 1:S){
        
      val1s =   beta1*s+beta0*futil_maint[s]+sh[,1] #COST of MAINTAINING
      val2s =   beta2 + beta0*futil_repl[s]+sh[,2]  #COST of REPLACING
      
      emax[s] = sum(apply(cbind(val1s,val2s),1,min))/nsim1
    }
    
  } else {
    
    val1 =   beta1st*(1:S)+beta0*futil_maint        #COST of MAINTAINING
    val2 =   beta2st + beta0*futil_repl             #COST of REPLACING
  
    emax = apply(cbind(val1,val2),1,min) 
  }
  
  #Store the emax values for each draw
  value[,igibb] = emax

  # --------------------------------------------------------
  # Using the kernel function, update the expected value matrix
  # --------------------------------------------------------

  #We check the past Ns draws (including the current one)
  Ns = min(igibb,NsAll) #this is the most we can check back over, which might be less than NsAll
  L =  min(igibb,LAll)  #this is the most we can average back over, which might be less than LAll
  
  #These are the indices we will check over
  indStart = igibb-Ns+1
  indEnd   = igibb
  ind = indStart:indEnd
  
  #----------------------------------
  # Derive the Gaussian Kernel values
  #----------------------------------
  const = (2*pi)^(-L/2) #not really needed, since only relative kernel values matter
  
  beta1_k   =  param[ind,1]
  beta2_k   =  param[ind,2]
  sigbeta1_k = param[ind,3]
  sigbeta2_k = param[ind,4]
    
  rkern11 = (1/hkerns1)*exp(-.5*(((beta1st-beta1_k)/hkerns1)^2))
  rkern22 = (1/hkerns2)*exp(-.5*(((beta2st-beta2_k)/hkerns2)^2))
  
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
  
  emax = value[,icheck,drop=FALSE]%*%kernMult #weight the past values by the kernel density, more weight given to the closer values
  
  #==========================================================
  # Now, generate nsim new data.
  #==========================================================
  
  #(omitted)
  
  igibb = igibb+1
}

plot(bayes1.out[,2],type="l")
plot(bayes1.out[,3],type="l")
print(mean(bayes1.out[,2]))
print(mean(bayes1.out[,3]))
print(round(accept/ngibb,2))                    #Acceptance Rate
matplot(bayes1.out[(ngibb/2):ngibb,-1],type="l")#Convergence
