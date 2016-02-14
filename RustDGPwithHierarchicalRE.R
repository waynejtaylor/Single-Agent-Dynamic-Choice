library(msm)

#This is only a portion of the original file, stripped down to the essentials

#version 2/9/2016

myopic_costs=function(S, MF, params, p){
  
  "This function computes the myopic expected cost associated with each decision for each state, 
   and returns an array of state/decision costs.
    
   Takes:
    * An integer S, describing the possible states of the bus
    * A maintenance cost function MF, which takes a vector of parameters and a 'state' argument
    * A vector params, to be supplied to the maintenance cost function MF. The first element of 
      the vector is the replacement cost rc.
    * A (3x1) vector p describing the state transitions probabilities 
        
    Returns:
    * A (Sx2) array containing the maintenance and replacement costs for the N possible states of the bus"
  
  rc = params[1]
  maint_cost = rep(NA,S)
  repl_cost = rep(NA,S)
  
  for(s in 1:S){
    maint_cost[s] = MF(s,params[-1])
    repl_cost[s] = rc
  }
  
  cbind(maint_cost,repl_cost)
}

lin_cost=function(s,params) s*params[1]

choice_prob=function(cost_array){
  
  # Returns the probability of each choice, conditional on an array of state/decision costs.
    
  S = nrow(cost_array)
  cost = cost_array-apply(cost_array,1,min) #take the difference since 1) results are the same 2) more stable with exp()
  util = exp(-cost)
  pchoice = util/rowSums(util)
  
  pchoice
}

contraction_mapping=function(S, p, MF, params, beta=0.75, threshold=1e-6, suppr_output=FALSE){
  
  "Compute the non-myopic expected value of the agent for each possible decision and each possible 
  state of the bus.
  Iterate until the difference in the previously obtained expected value and the new expected value
  is smaller than the threshold.

  Takes:
  * A finite number of states S
  * A state-transition probability vector p = [p(0), p(1), p(2), ..., p(k)] of length k < S
  * A maintenance cost function MF
  * A vector params for the cost function
  * A discount factor beta (optional)
  * A convergence threshold (optional)

  Returns:
  * The converged choice probabilities for the forward-looking and myopic agents for each state, 
  conditional on 'params'"

  achieved = TRUE
  
  # Initialization of the state-transition matrices
  # ST_mat: describe the state-transition probabilities if the maintenance cost is incurred
  # RT_mat: regenerate the state to 0 if the replacement cost is incurred.
  # [a,b] = transition from state "a" to "b"
  ST_mat = matrix(0,S,S)
  lp = length(p)
  for(i in 1:S){
    for(j in 1:lp){
      if((i+j-1)<S)  ST_mat[i,i+j-1] = p[j]
      if((i+j-1)==S) ST_mat[i,i+j-1] = sum(p[j:lp]) #out of columns, so collapse the probabilities
    }
  }
  
  R_mat = cbind(1,matrix(0,S,S-1))
  
  # Initialization of the expected value (which is also the myopic decision cost of the agent).
  # Here, the forward-looking component is initialized at 0.
  k = 0
  EV = matrix(0,S,2)
  EV_myopic = EV_new = myopic_costs(S, MF, params, p)
  
  # Contraction mapping loop
  while(max(abs(EV_new-EV)) > threshold){
    # Store the former expected value
    EV = EV_new
    # Obtained the probability of maintenance and replacement from the former expected value
    pchoice = choice_prob(EV)
    # Compute the expected cost for each state: Nx1 vector
    ecost = rowSums(pchoice*EV)
    # Compute the two components of forward-looking utility: In case of maintenance, 
    # utility of future states weighted by transition probabilities. In case of replacement,
    # the future utility is the utility of state 0
    futil_maint = ST_mat%*%ecost
    futil_repl = R_mat%*%ecost
    futil = cbind(futil_maint,futil_repl)
    # Future utility is discounted by beta, and added to the myopic cost. 
    EV_new = EV_myopic + beta*futil
    k = k+1
    if(k == 1000) achieved = FALSE
  }

  if(!suppr_output){
    if(achieved){
      cat("Convergence achieved in ",k," iterations")
    } else {
      cat("CM could not converge! Mean difference = ",round(mean(EV_new-EV),2))
    }
  }

  list(CP_forward=choice_prob(EV_new),CP_myopic=choice_prob(EV_myopic))
}

bus_dgp=function(nbus,nperiod,pchoice){
  
  "Genearates a simulated bus decision process for nbus buses over nperiod decisions
  #using the converged choice probabilities pchoice
  
  Takes:
  * nbus:    The number of buses to simulate
  * nperiod: The number of decision periods
  * pchoice: (S x 1) vector of probabilities. Represents probability of 'maintaining' the bus in each state (versus replacing).
  
  Returns:
  * A dataframe of size (nbus*nperiod x 4) with columns id, choice, mileage, state"
  
  nobs = nbus*nperiod
  bus_df = data.frame(id=rep(1:nbus,times=nperiod),
                      choice=rep(NA,nobs),
                      mileage=rep(NA,nobs),
                      state=rep(NA,nobs))
  
  #Initialization
  bus_df[1:nbus,1] = 1:nbus
  bus_df[1:nbus,-1] = 0
  
  i = 1 
  while(i<=nobs){
    prev_ind = seq(from=i,by=1,length.out=nbus) #index of the "last" decisions
    next_ind = prev_ind+nbus                    #index of the "next" decisions
    
    prev_mileage = bus_df$mileage[prev_ind]
    prev_state = bus_df$state[prev_ind]  
    
    # Generating choices from choice probabilities, conditional on the state of each bus
    # Make a decision to maintain the bus using the probability in pchoice and a current state s
    choices = (runif(nbus)>pchoice[prev_state+1])*1 #Note: state "0" corresponds to row "1" in prev_state
    bus_df$choice[prev_ind] = choices
    
    if(next_ind[1]<=nobs){
      new_mileage = (1-choices)*prev_mileage + rtnorm(nbus,mu,sigma,lower,upper)
      new_state = floor(new_mileage/5000)
      
      bus_df$mileage[next_ind] = new_mileage
      bus_df$state[next_ind] = new_state  
    }
    
    i=i+nbus
  }
  
  bus_df
}

# IV. Generation of simulated data (using a linear cost specification)

# 1) True coefficients

# Cost and discount coefficients:
  rc=20
  theta11=0.5
  beta=0.75

# State transition probabilities: (Mileage_t+1-Mileage_t)~NTrunc(6000,4000), such that:
  p_x0=0.36
  p_x1=0.48
  p_x2=0.16

# Number of states
  S = 70

  p = c(p_x0, p_x1, p_x2)

  lower=0
  upper=15000
  mu=6000
  sigma=4000

#  Simulations
  nbus = 100
  nperiod = 1000

# We are simulating hierarchical data
# Both "rc" and "theta11" will be functions of a buses individual characteristics
# To give it context in this example, think of Z containing a measure of "weather harshness"

noMix = FALSE #If TRUE, Z only contains a matrix of 1's

#B=ZDelta + U

nvar = 2 #number of coefficients: rc and theta11
nz   = 2 #number of regressors in mixing distribution 
Z = cbind(rep(1,nbus),runif(nbus,.6,1))
  
Delta = matrix(c(0,25,0,.625),nrow=nz,ncol=nvar)
  #Delta = matrix(c(0,25,.5,0),nrow=nz,ncol=nvar) #For now, fix theta11 at .5
  #The way to interpret this matrix:
    #Delta[1,1]: the incremental effect of Z[i,1] on theta 1 (rc) (here it acts as an intercept)
    #Delta[2,1]: the incremental effect of Z[i,2] on theta 1 (rc)
    #Delta[1,2]: the incremental effect of Z[i,1] on theta 2 (theta11) (here it acts as an intercept)
    #Delta[2,2]: the incremental effect of Z[i,2] on theta 2 (theta11)
colMeans(Z%*%Delta) #notice the means are the same as the original values
hist(Z%*%Delta[,1]) #the disributions straddle the original values
hist(Z%*%Delta[,2])

#No heterogeneity
if(noMix){  
  nz   = 1
  Z = cbind(rep(1,nbus))
  Delta = matrix(c(20,.5),nrow=nz,ncol=nvar)
}

Vtheta=diag(c(.1^2,.01^2))

#Simulate the data
thetas = matrix(0,nbus,nvar)
bus_df = NULL
for(i in 1:nbus){
  
  thetai=t(Delta)%*%Z[i,]+as.vector(t(chol(Vtheta))%*%rnorm(nvar))
  print(thetai)
  thetas[i,] = thetai
  params = as.vector(thetai)
  
  out = contraction_mapping(S=S, p=p, MF=lin_cost, params=params, beta = beta,suppr_output=TRUE)
  pchoice=out$CP_forward[,1]
  bus_df_i = bus_dgp(1,nperiod,pchoice)
  bus_df_i$id = i
  
  bus_df = rbind(bus_df,bus_df_i)
  cat(i,fill=TRUE)
}

#Export the data
write.csv(bus_df,file="bus_df_RE.csv",row.names = FALSE)
save(Z,Delta,Vtheta,thetas,file="simdata.Rdata")