#see https://github.com/QuentinAndre/John-Rust-1987-Python/blob/master/Data%20Generation%20and%20Likelihood%20Fit.ipynb
library(ggplot2)
library(msm)

#1) Number of States and Transitions
# In this paper, the agent observes each bus at a given time t, and makes 
# a decision regarding the replacement of the engine or its maintenance as a 
# function of the observed state of the engine (i.e. its mileage).

# As in the original paper, we discretize the number of states into 
# intervals of length 5000, and define a transition density which governs 
# the evolution of mileage from one time period to another:

  # p(x_t+1|x_t,i_t,theta_3) = g(x_t+1-x_t,theta_3) if i_t=0
  #                          = g(x_t+1-0,theta_3) if i_t=1

# Here, the function g is defined by a multinomial distribution on 
# the set {0, 1, 2}. In other words, the probability that a bus'
# mileage will increase by a value x between now and the next 
# maintenance decision is:

  # theta_31 if x in c(0,5000)
  # theta_32 if x in c(5001,10000)
  # theta_33 = 1-(theta_31-theta_32) if  x > 10000

#2) Decision Criteria

# The utility function of the decision maker for a single time period is:

  # u(x_t,i,theta_1) = -c(x_t,theta_1) + e_t(0) if i = 0
  #                  = -RC-c(x_t,theta_1) + e_t(1) if i = 1

# The agent will chose the decision which maximizes its utility. 
# If the agent is forward looking (i.e. if beta != 0), he does 
# not only maximize his utility in the present time period, but 
# also his future utility discounted by beta

# Since beta cannot be estimated without knowing the 
# utility function, we will later set beta to its true value.

# 3) Cost Function 

# A functional form has to be imposed on the cost function c(x_t).

# In this recreation, we assume that the replacement cost is a constant RC, 
# and allow the maintenance cost MC(s,theta_1) to be a linear, exponential 
# or logarithmic function of the state of the bus scaled by a 
# parameter theta_1

#DATA GENERATING PROCESS--------------------------------------------

# 1) Bus Mileage

# A. Introduction
# The generation of a fake dataset can be done in two different ways:
#  -Generating monthly mileage increase from a continuous distribution, discretize the mileage, 
#   and recover the state transition parameters (theta_3) using a maximum-likelihood estimation 
#   (as done in the paper).
#  -Set the vector theta_3 to an arbitrary value, and generate discrete states transitions 
#   with corresponding probabilities.

# The advantage of the first method is that this is closer to what is actually done in the paper: 
# indeed, the focal goal of the paper is to recover RC and the cost scaling vector of parameters theta_1.

# If the decision of the agent is based on the discretized state (and not on the actual mileage), 
# then this approach is not noisier than the second one.

# B. Mileage Transition Definition
# For this reason, we will adopt this first approach, and assumes that the increase in mileage
# follows a truncated normal distribution bounded at [0, 15000]:

  # (M_t+1-M_t)~NTrunc(6000,4000)

# The parameters are arbitrarily chosen to ensure that the state transitions probabilities are plausible:

lower=0
upper=15000
mu=6000
sigma=4000

p_x0 = ptnorm(5000,mu,sigma,lower,upper)
p_x1 = ptnorm(10000,mu,sigma,lower,upper)-p_x0
p_x2 = 1 - p_x1 - p_x0
p = c(p_x0, p_x1, p_x2)

p

#p[1] chance of keeping the same mileage, p[2] chance of increasing the mileage state by 1, p[3] chance of increasing the mileage state by 2

#2) Agent's Decisions

# A. Initialization of the variables
# The replacement cost, the maintenance cost function parameters, and the discount rate
# are defined in this section. We set up the discount rate to be lower than in the paper
# to speed up the convergence of the contraction mapping algorithm

rc = 20
theta1_1 = 0.5
theta1_2 = 0.01
beta = 0.75

# B. Definition of the cost function

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

# We assume different possible forms for the maintenance cost function:
# A linear form: MF(s,theta)=theta_11*s
# A quadratic form: MF(s,theta)=theta_11*s+theta_12*s^2
# An exponential form: MF(s,theta)=exp(theta_11*s)
# A log form: MF(s,theta)=log(theta_11+theta_12*s)

lin_cost=function(s,params) s*params[1]
quad_cost=function(s,params) s*params[1]+s^2*params[2]
exp_cost=function(s,params) exp(s*params[1])
log_cost=function(s,params) log(params[1] + s*params[2])

# C. Definition of the choice probabilities, as a function of an array of costs

choice_prob=function(cost_array){
  
  # Returns the probability of each choice, conditional on an array of state/decision costs.
    
  S = nrow(cost_array)
  cost = cost_array-apply(cost_array,1,min) #take the difference since 1) results are the same 2) more stable with exp()
  util = exp(-cost)
  pchoice = util/rowSums(util)
  
  pchoice
}

# D. The Contraction Mapping Algorithm (generate the forward-looking choice probabilities)

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

# IV. Generation of simulated data (using a linear cost specification)

# 1) True coefficients

# Here, we use the same coefficients as the one described earlier:

  # Cost and discount coefficients:
    rc=20
    theta11=0.5
    theta12=0.01
    beta=0.75

  # State transition probabilities: (Mileage_t+1-Mileage_t)~NTrunc(6000,4000), such that:
    p_x0=0.36
    p_x1=0.48
    p_x2=0.16

# Finally, we assume that there are 70 discrete states for the bus' mileage 
# (i.e. no bus has a mileage greater than 345,000 in the dataset).

# 2) Choice probabilities, as a function of the bus' state:

# Using the Contraction Mapping algorithm implemented earlier, we can obtain the 
# probabibility of maintenance as a function of the bus' state for a forward-looking
# agent and a myopic agent.

params_lin = c(rc,theta11)
p = c(p_x0, p_x1, p_x2)
S = 70
out = contraction_mapping(S=S, p=p, MF=lin_cost, params=params_lin, beta = .75)
lin_forward=out$CP_forward
lin_myopic=out$CP_myopic
pchoice = lin_forward[,1]

ggdat1 = data.frame(decisionRule=c(rep("Forward-Looking",nrow(lin_forward)),
                                   rep("Myopic",nrow(lin_forward))),
                    pMaint=c(lin_forward[,1],lin_myopic[,1]),
                    State=rep(1:S,times=2))

ggplot(ggdat1,aes(y=pMaint,x=State,color=decisionRule))+geom_line(lwd=1)+theme_bw(20)+xlim(5,50)

# 3) The bus replacement dataset.

# A. Initialization

# First, we decide on a number of bus, that we initialize in state 0.
# We treat the dataset as an N_Bus x 4 array, containing:
# -The identifier of the buses
# -The choice made for the buses at this time
# -The mileage
# -The states

# B. State transition functions

# The dataset of buses transition from one period to another in the following way:
# -The agent observes the current state of the buses, and take a decision to maintain (0) or replace (1) each bus.

# -The mileage of each bus is updated a first time:
#   -Kept at its current value if the agent decided to maintain the engine
#   -Reset to 0 if the agent decided to replace the engine

# -The mileage of each bus is updated a second time : a random draw is taken from the 
#   mileage increment random variable for each bus, and the corresponding mileage is added 
#   to the mileage of the bus.

# -The state of each bus is updated to reflect its new mileage

# We define below the transition function

transition=function(bus_df, pchoice){

  #Return the updated bus dataset after one decision of our agent.
  
  #Takes:
  # bus_df : A dataframe of buses, containing the identifier of the buses, their mileages, and their current state.
  # *The columns are named: "id","choice","milage","state"
  # pchoice : The converged probabities of the agent choosing to maintain the bus (Sx1 vector)

  #Returns:
  #The updated dataset of buses, with the new decisions appended at the end of the dataframe.
  
  # Recovering the number of buses, the previous mileage and the previous states of the buses
  nbus = max(bus_array$id)
  obs = nrow(bus_array)
  ind = (obs-nbus+1):obs #index to capture the previous 'nbus' rows
  prev_mileage = bus_array$mileage[ind]
  prev_states = bus_array$state[ind]
  # Generating choices from choice probabilities, conditional on the state of each bus
  #Make a decision to maintain or replace the bus, based on a probability of choice p and a current state s
  choices = (runif(nbus)>pchoice[prev_states+1])*1

  # Generating the new mileage and state
  new_mileage = (1-choices)*prev_mileage + rtnorm(nbus,mu,sigma,lower,upper)
  new_states = floor(new_mileage/5000)
  new_array = data.frame(id=1:nbus,choice=rep(0,nbus),mileage=new_mileage,state=new_states)
  bus_array$choice[ind]=choices #put the choices in the "last" set of observations
  
  rbind(bus_array, new_array)
}

# C. Generation of the full dataset

# The complete dataset is generated by making the agent take a series of decisions.
# Here, we assume that the agent takes nbus consecutive decisions.
# The total number of data points generated will be equal to nbus*nperiod (nbus buses and nperiod time periods).

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
  while(i<nobs){
    prev_ind = seq(from=i,by=1,length.out=nbus) #index of the "last" decisions
    next_ind = prev_ind+nbus                    #index of the "next" decisions
      
    prev_mileage = bus_df$mileage[prev_ind]
    prev_state = bus_df$state[prev_ind]  
      
    # Generating choices from choice probabilities, conditional on the state of each bus
    # Make a decision to maintain the bus using the probability in pchoice and a current state s
    choices = (runif(nbus)>pchoice[prev_state+1])*1 #Note: state "0" corresponds to row "1" in prev_state
    bus_df$choice[prev_ind] = choices
    
    if(next_ind[1]<nobs){
      new_mileage = (1-choices)*prev_mileage + rtnorm(nbus,mu,sigma,lower,upper)
      new_state = floor(new_mileage/5000)
      
      bus_df$mileage[next_ind] = new_mileage
      bus_df$state[next_ind] = new_state  
    }
    
    i=i+nbus
  }
  
  bus_df
}

nbus = 1000
nperiod = 100
bus_df_lin = bus_dgp(nbus,nperiod,pchoice)

choice_freq=aggregate(choice~state,data=bus_df_lin,FUN=mean)

ggdat2 = data.frame(dataType=c(rep("Observed behavior",nrow(choice_freq)),
                               rep("True choice probabilities",length(pchoice))),
                    pMaint=c(1-choice_freq$choice,pchoice),
                    State=c(choice_freq$state,0:(S-1)))

ggplot(ggdat2,aes(y=pMaint,x=State,lty=dataType,color=dataType))+geom_line(lwd=1)+theme_bw(20)+xlim(1,25)

# C. Export of the data

# We only need the decision state and the choice, as the other variables will not be used in the estimation section.

write.csv(bus_df_lin,file="bus_df_lin.csv",row.names = FALSE)

# VI. Estimation

# 1) Dynamic Utility Logit Model

# In this section, we implement a Dynamic Utility Logit model, which will be our workbench for the estimation
# of the data. We wrap in this class the functions defined earlier (to avoid global references)
# and refactor some portions of the code (to make the model more flexible). 
# The maintenance cost function is kept outside of the scope of the class, and must be 
# supplied at the initialization of the object as an argument.

data=bus_df_lin

DynamicLogit=function(params,data,S,p,MF){

  "
  Evaluate the cost parameters underlying a bus replacement pattern by a forward-looking agent.
  
  Takes:
  * Data: a dataframe, which contains:
    -choice: the name of the column containing the dummy endogenous variable
    -state: the name of the column containing the exogenous variable 
  
  * p: The state-transition vector of exogenous variable.
      For instance, p = [0, 0.6, 0.4] means that the bus will 
      transition to the next mileage state with probability 0.6, 
      and to the second next mileage state with probability 0.4.
  
  * MF: A function passed as an argument, which is the functional 
      form for the maintenance cost. This function must accept as
      a first argument a state s, and as a second argument a vector of parameters.
  "
  
  endog = data$choice
  exog = data$state

  N=length(endog)
  S=max(exog)*2 # Assumes that the true maximum number states is twice the maximum observed state.

  # Matrices to speed up computations of the log-likelihood
  
  # A (SxN) matrix indicating the state of each observation
  state_mat=matrix(0,S,N)
  for(s in 0:(S-1)) state_mat[s+1,]=(exog==s)*1 #Note 0 is a state, sum(state_mat)==N should be true
  
  # A (2xN) matrix indicating with a dummy the decision taken by the agent for each time/bus observation (replace or maintain)
  dec_mat = rbind(t(1-endog),endog)
  
  "
  The log-likelihood of the Dynamic model is estimated in several steps.
  1) The current parameters are supplied to the contraction mapping function
  2) The function returns a matrix of decision probabilities for each state.
  3) This matrix is used to compute the loglikelihood of the observations
  4) The log-likelihood are then summed accross individuals, and returned
  "
  
  util = contraction_mapping(S=S, p=p, MF=MF, params=params, beta = .75,suppr_output = TRUE)
  pchoice = util$CP_forward
  logprob = log(t(state_mat)%*%pchoice)
  -sum(logprob*t(dec_mat))
}

# 2) Fitting the linear cost data

# A. Fitting the true linear cost function to the data

# In this section, we fit the data generated by the linear cost function and recover the parameters RC and ??. We do it for different characterizations of the cost function, and thereby illustrate the consequences of a misspecification.
# Recall that in the linear model:
  # RC=20
  # theta11=0.5

bounds = c(1e-6, Inf)
npars=2
lin_fit = optim(par=rep(.1,npars),fn=DynamicLogit,method=c("L-BFGS-B"),lower=bounds[1],upper=bounds[2],
               data=data,S=S,p=p,MF=lin_cost,control=list(fnscale=1))

# Return the parameters obtained after fitting the likelihood function to the data.
loglike =  lin_fit$value
fit_params = lin_fit$par
cat("Log-Likelihood: ",loglike,fill=TRUE)
cat("RC: ",fit_params[1],fill=TRUE)
cat("thetas: ",fit_params[-1],fill=TRUE)

#Compare with the true values
params_lin = c(20,.5)
linEst = contraction_mapping(S=70, p=p, MF=lin_cost, params=fit_params, beta = .75)
lin_forwardEst=linEst$CP_forward
lin_myopicEst=linEst$CP_myopic

gglinEst = data.frame(decisionRule=c(rep("Forward-Looking (Lin)",nrow(lin_forward)),
                                     rep("Myopic (Lin)",nrow(lin_myopic)),
                                     rep("Forward-Looking (Lin Est.)",nrow(lin_forwardEst)),
                                     rep("Myopic (Lin Est.)",nrow(lin_myopicEst))),
                    pMaint=c(lin_forward[,1],lin_myopic[,1],lin_forwardEst[,1],lin_myopicEst[,1]),
                    State=c(0:(nrow(lin_forward)-1),0:(nrow(lin_myopic)-1),0:(nrow(lin_forwardEst)-1),0:(nrow(lin_myopicEst)-1)))

ggplot(gglinEst,aes(y=pMaint,x=State,color=decisionRule))+geom_line(lwd=1)+theme_bw(20)+xlim(5,50)