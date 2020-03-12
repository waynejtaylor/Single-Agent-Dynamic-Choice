#see https://github.com/QuentinAndre/John-Rust-1987-Python/blob/master/Data%20Generation%20and%20Likelihood%20Fit.ipynb
#Goal: replicate Monte Carlo Simulation section of HMSS 1994

#SUPPORT FUNCTIONS -----------------------------------------------------

#Choice = 1 means replace
lin_cost  = function(s,params,choice) choice*params[1] + (1-choice)*s*params[2]

myopic_costs = function(S, MF, params, p){
  
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
  
  # maint_cost = rep(NA,S)
  # repl_cost = rep(NA,S)
  # 
  # for(s in 1:S){
  #   maint_cost[s] = MF(s,params,0)
  #   repl_cost[s]  = MF(s,params,1)
  # }
  
  maint_cost = MF(1:S,params,0)
  repl_cost  = MF(1:S,params,1)
  
  cbind(maint_cost,repl_cost)
}

#Choice probability function
choice_prob = function(cost_array){
  
  # Returns the probability of each choice, conditional on an array of state/decision costs.
  
  S    = nrow(cost_array)
  cost = cost_array-apply(cost_array,1,min) #take the difference since 1) results are the same 2) more stable with exp()
  util = exp(-cost)
  pchoice = util/rowSums(util)
  
  pchoice
}

#Contraction Mapping Algorithm (generate the forward-looking choice probabilities)
contraction_mapping = function(S, p, MF, params, beta=0.75, threshold=1e-6, suppr_output=FALSE){
  
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
  k  = 0
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

#DGP Function
bus_dgp = function(nbus,nperiod,pchoice,p){
  
  "Genearates a simulated bus decision process for nbus buses over nperiod decisions
  #using the converged choice probabilities pchoice
  
  Takes:
  * nbus:    The number of buses to simulate
  * nperiod: The number of decision periods
  * pchoice: (S x 1) vector of probabilities of 'replacing' the bus given each state S (versus 'maintaining').
  * p: vector of state transition probabilities  
  
  Returns:
  * A dataframe of size (nbus*nperiod x 3) with columns id, state, choice
  * Note: choice = 1 if the decision is to replace, 0 = maintain
  
  Notes:
  * The state is observed in time period t, and then the choice is made in that same period
  * The updated state will be realized in the NEXT time period t + 1"
  
  nobs = nbus*nperiod
  bus_df = data.frame(id=rep(1:nbus,times=nperiod),
                      state=rep(NA,nobs),
                      choice=rep(NA,nobs))
  
  #Initialization
  bus_df[1:nbus,"id"] = 1:nbus
  bus_df[1:nbus,"state"] = 1
  bus_df[1:nbus,"choice"] = 0
  
  i = 1 
  while(i<nobs){
    prev_ind = seq(from=i,by=1,length.out=nbus) #index of the "last" decisions
    next_ind = prev_ind+nbus                    #index of the "next" decisions
    
    prev_states = bus_df$state[prev_ind]  
    new_states  = prev_states
    
    # Generating choices, where 1 = replace, from choice probabilities, conditional on the state s of each bus
    choices = (pchoice[prev_states] > runif(nbus))*1 
    bus_df$choice[prev_ind] = choices
    
    if(next_ind[1]<nobs){
      
      #When the choice is to maintain (choice = 0), advance the mileage state, otherwise reset to zero
      #note: choices = 1 with replace, which is the minimum starting state
      new_states = (1-choices)*prev_states + choices + sample(0:2,nbus,prob = p,replace=TRUE)
      
      bus_df$state[next_ind] = new_states
      
    }
    
    i = i + nbus
  }
  
  bus_df
}

#Log-Likelihood|theta for forward looking agent
dynamicLogit = function(params,data,S,p,MF,beta=.75){
  
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
  exog  = data$state
  
  N = length(endog)
  S = max(exog)*2 # Assumes that the true maximum number states is twice the maximum observed state.
  
  # Matrices to speed up computations of the log-likelihood
  
  # A (SxN) matrix indicating the state of each observation
  state_mat=matrix(0,S,N)
  for(s in 1:S) state_mat[s,]=(exog==s)*1 #sum(state_mat)==N should be true
  
  # A (2xN) matrix indicating with a dummy the decision taken by the agent for each time/bus observation (replace or maintain)
  dec_mat = rbind(t(1-endog),endog)
  
  "
  The log-likelihood of the Dynamic model is estimated in several steps.
  1) The current parameters are supplied to the contraction mapping function
  2) The function returns a matrix of decision probabilities for each state.
  3) This matrix is used to compute the loglikelihood of the observations
  4) The log-likelihood are then summed accross individuals, and returned
  "
  
  util = contraction_mapping(S=S, p=p, MF=MF, params=params, beta=beta,suppr_output = TRUE)
  pchoice = util$CP_forward
  logprob = log(t(state_mat)%*%pchoice)
  -sum(logprob*t(dec_mat))
}

#RUST DATA GENERATING PROCESS--------------------------------------------
#see page 280 in HMSS 1994

#Transition probabilities
#p[1] chance of keeping the same mileage
#p[2] chance of increasing the mileage state by 1
#p[3] chance of increasing the mileage state by 2

p_x0 = .349
p_x1 = .639
p_x2 = 1 - p_x1 - p_x0
p = c(p_x0, p_x1, p_x2)

p

#structural parameters
theta1 = 2.0   #replacement cost
theta2 = 0.09  #maintenance cost
# #OR
# theta1 = 8.0   #replacement cost
# theta2 = 0.09  #maintenance cost

#other parameters
beta   = 0.9   #discount factor (.9 in paper)
S      = 90    #discretized states (90 in paper)
nbus    = 100
nperiod = 10

#probability of replacing|state
theta = c(theta1,theta2)

out_cm   = contraction_mapping(S=S,p=p,MF=lin_cost,params=theta,beta=beta)

cm_forward = out_cm$CP_forward
cm_myopic  = out_cm$CP_myopic

#probability of replace in each of the S states
pchoice = cm_forward[,2]
#pchoice = seq(0,1,length.out = S)

#Generate Simulated Data----
busData = bus_dgp(nbus,nperiod,pchoice,p)

# #Does pchoice look similar to empirical pr(choice|state)?
# aggregate(choice~state,busData,FUN=mean)
# pchoice

#MLE -------------------------------------------------------------
bounds = c(1e-6, Inf)
fit_mle = optim(par=rep(.1,2),fn=dynamicLogit,method=c("L-BFGS-B"),lower=bounds[1],upper=bounds[2],
                data=busData,S=S,p=p,MF=lin_cost,beta=beta,control=list(fnscale=1),hessian = FALSE)
fit_mle$par
theta

#HM CCP ------------------------------------------------
#https://notes.quantecon.org/submission/5c4906d9f68373000f919cd0
#From HMSS 1993, use conditional choice probabilities

gamma_cnst = 0.577

#Rust calculated the value function inside of a logit

#But this can be slow. 
#In CCP, we instead start with a simple estimate of choice probabilities 
#and then adjust to account for the future.
thetahat = theta

#transition probabilities (usually estimated from data)
phat = p

#choice probabilities (usually estimated from data)
pchoicehat = pchoice

#transition matrices (usually estimated from data)
#probability of moving from state column to state column, 
#given replace = 1 or maintain = 0

#when decision is to replace, we reset each state
transMat1 = matrix(c(phat,rep(0,S-length(phat))),S,S,byrow = TRUE)
#transMat1 = cbind(1,matrix(0,S,S-1))

#when the decision is to maintain, we increment according to phat
# transMat0 = matrix(0,S,S)
# for(s in 1:(S-length(phat) + 1)) transMat0[s,s:(s+2)] = phat
# #fill in last 2 manually
# transMat0[S-1,(S-1):S] = phat[1:2]
# transMat0[S,S]         = phat[1]

transMat0 = matrix(0,S,S)
lp = length(phat)
for(i in 1:S){
  for(j in 1:lp){
    if((i+j-1)<S)  transMat0[i,i+j-1] = phat[j]
    if((i+j-1)==S) transMat0[i,i+j-1] = sum(phat[j:lp]) #out of columns, so collapse the probabilities
  }
}

#Using:
#1) the payoff function (whose parameters we want to learn)
#2) the transition matrices: transMat1, transMat2
#3) the initial choice probabilities pchoicehat
#we can now calculate the bus manager's value function using the following formula:


hm_value = function(thetahat,lin_cost,pchoicehat,transMat0,transMat1){

  lambda = matrix(1,1,S)
  
  denom  = diag(S) - beta*((1-pchoicehat)%*%lambda*transMat0 + pchoicehat%*%lambda*transMat1)
  
  numer = pchoice*(-lin_cost(1:S,thetahat,1) + gamma_cnst - log(pchoice)) + 
      (1-pchoice)*(-lin_cost(1:S,thetahat,0) + gamma_cnst - log(1-pchoice))
  
  value = solve(denom) %*% numer
  
  value
}


# #alternative specification (same output)
# hm_value = function(thetahat,lin_cost,pchoicehat,transMat0,transMat1){
#     
#   pchoicehatMat = matrix(pchoicehat,S,S)
#   
#   denom = diag(S) - beta*(1-pchoicehatMat)*transMat0 - beta*pchoicehatMat*transMat1 
#     
#   # #Original: but note UTILITY is NEGATIVE of COST
#   # numer = (1-pchoicehat)*(lin_cost(1:S,thetahat,0) + gamma_cnst - log(1 - pchoicehat)) +
#   #             pchoicehat*(lin_cost(1:S,thetahat,1) + gamma_cnst - log(pchoicehat))
#     
#   numer = (1-pchoicehat)*(-lin_cost(1:S,thetahat,0) + gamma_cnst - log(1 - pchoicehat)) +
#               pchoicehat*(-lin_cost(1:S,thetahat,1) + gamma_cnst - log(pchoicehat))
#   
#   value = solve(denom) %*% numer
#     
#   value
# 
# }

#replacement probabilities
hm_prob = function(thetahat,lin_cost,pchoicehat,transMat0,transMat1){
  
  value = hm_value(thetahat,lin_cost,pchoicehat,transMat0,transMat1)
  value = value - min(value) #subtract the smallest value
  
  # delta1 = exp(lin_cost(1:S,thetahat,1) + beta*transMat1%*%value)
  # delta2 = exp(lin_cost(1:S,thetahat,0) + beta*transMat0%*%value)
  
  delta1 = exp(-lin_cost(1:S,thetahat,1) + beta*transMat1%*%value)
  delta2 = exp(-lin_cost(1:S,thetahat,0) + beta*transMat0%*%value)
  
  probhat = delta1/(delta1 + delta2)

  probhat

}

#Create LL function
nloglikeobs = function(thetahat,probhat=pchoicehat) {
  
  y = busData$choice
  X = busData$state
  
  prob = hm_prob(thetahat,lin_cost,probhat,transMat0,transMat1)
  prob = prob[X]

  LL = (1-y)*log(1-prob) + y*log(prob)
  
  -sum(LL)
 
}


fit_CCP  = optim(c(0,0),nloglikeobs,method = "BFGS",probhat=pchoicehat)
fit_CCP$par
theta

# #other methods (very similar results)
# bounds = c(1e-6, Inf)
# fit_CCP  = optim(c(0,0),nloglikeobs,method=c("L-BFGS-B"),lower=bounds[1],upper=bounds[2],probhat=pchoicehat)
# fit_CCP$par
# fit_CCP  = optim(c(0,0),nloglikeobs,method=c("Nelder-Mead"),probhat=pchoicehat)
# fit_CCP$par

fit_CCP$value
nloglikeobs(theta,pchoice)

#iterate
K = 5
thetahat = fit_CCP$par
probhat  = pchoicehat
k = 1
while(k < K) {
  
  #restimate choice probabilities
  probhat  = hm_prob(thetahat,lin_cost,probhat,transMat0,transMat1)
  
  #refit CCP
  ccp = optim(c(0,0),nloglikeobs,method = "BFGS",probhat=probhat)

  #update thetahat
  thetahat = ccp$par
  
  k = k + 1
    
}
thetahat

#HMSS 1994: CCP USING FORWARD SIMULATION ------------------

#From above, these are the estimated CCP probabilities that will be used for comparisons
probhat_CCP = hm_prob(theta,lin_cost,pchoice,transMat0,transMat1)
value_CCP   = hm_value(theta,lin_cost,pchoicehat,transMat0,transMat1)
value_CCP   = value_CCP - min(value_CCP) #subtract the smallest value
summary(beta*transMat1%*%value_CCP)
summary(beta*transMat0%*%value_CCP)
head(cbind(beta*transMat0%*%value_CCP,beta*transMat1%*%value_CCP))

#First, construct the conditional value function utility + beta*V as if theta was known
valueFS       = futureValFS(S,T=75,R=200,beta,p,pchoice,theta)
valueFS_util1 = -lin_cost(1:S,theta,1) + valueFS[,2]    #Current utility plus discounted future utilities given t0 action = replace
valueFS_util0 = -lin_cost(1:S,theta,0) + valueFS[,1]    #Current utility plus discounted future utilities given t0 action = maintain
probhat_theta = exp(valueFS_util1)/(exp(valueFS_util0) + exp(valueFS_util1))

#Note similarity in range
range(beta*transMat1%*%value_CCP - beta*transMat0%*%value_CCP)
range(valueFS[,2] - valueFS[,1])

#Next, forward simulate the future discounted utilities separating out theta and differening
valueFS_HMSSFuture      = futureValHMSS(S,T=75,R=200,beta,p,pchoice)
valueFS_HMSSFuture_diff = (-lin_cost(1:S,theta,1) - -lin_cost(1:S,theta,0) + valueFS_HMSSFuture %*% c(1,theta)) 
probhat_HMSSFuture      = exp(valueFS_HMSSFuture_diff)/(1 + exp(valueFS_HMSSFuture_diff))

#Next, finally reconstruct the function to include t = 0 utility
valueFS_HMSS      = valHMSS(S,T=75,R=200,beta,p,pchoice)
valueFS_HMSS_diff = valueFS_HMSS %*% c(1,theta)
probhat_HMSS      = exp(valueFS_HMSS_diff)/(1 + exp(valueFS_HMSS_diff))

#Next, check that you can enter states as a vector (output still as mat)
valueFS_HMSS_ME_mat      = valHMSS_ME_mat(1:S,T=75,R=200,beta,p,pchoice)
valueFS_HMSS_ME_mat_diff = valueFS_HMSS_ME_mat %*% c(1,theta)
probhat_HMSS_ME_mat      = exp(valueFS_HMSS_ME_mat_diff)/(1 + exp(valueFS_HMSS_ME_mat_diff))

#Finally, see if the output can be an array
valueFS_HMSS_ME      = valHMSS_ME(1:S,T=75,R=200,beta,p,pchoice)
valueFS_HMSS_ME      = apply(valueFS_HMSS_ME,c(1,2),FUN=mean)
valueFS_HMSS_ME_diff = valueFS_HMSS_ME %*% c(1,theta)
probhat_HMSS_ME      = exp(valueFS_HMSS_ME_diff)/(1 + exp(valueFS_HMSS_ME_diff))

matplot(cbind(
  pchoice,             #Actual pchoice|state
  probhat_CCP,         #From CCP estimates
  probhat_theta,       #From forward simulation if theta is known
  probhat_HMSSFuture,  #From forward simulation separating out theta
  probhat_HMSS,        #Reconstruct function to include t = 0
  probhat_HMSS_ME_mat,
  probhat_HMSS_ME))

#Nice plot
library(ggplot2)
ggDat_probhat = data.frame(Method = rep(c("0 Nested FP","1 CCP","2 FS with Theta","3 FS Theta Separate","4 FS with t0"),each = S),
                   probhat  = c(pchoice,probhat_CCP,probhat_theta,probhat_HMSSFuture,probhat_HMSS),
                   state    = rep(1:S,5))
ggplot(ggDat_probhat,aes(state,probhat,color=Method)) + theme_bw(15) + geom_line()

#Use valueFS_HMSS to predict theta

#Unweighted:
y = log(pchoice/(1-pchoice)) - valueFS_HMSS[,1]
X = valueFS_HMSS[,-1]
solve(crossprod(X))%*%crossprod(X,y)
fit_CCP$par #CCP estimate
theta       #truth

#Weighted:
y = log(pchoice/(1-pchoice)) - valueFS_HMSS[,1]
X = valueFS_HMSS[,-1]
stateCount = aggregate(id~state,busData,FUN=length)
Wdiag = rep(0,S)
Wdiag[stateCount$state] = stateCount$id 
W = diag(Wdiag)
solve(t(X)%*%W%*%X)%*%(t(X)%*%W%*%y)
fit_CCP$par #CCP estimate
theta       #truth

#Same as weighted:
y = log(pchoice[busData$state]/(1-pchoice[busData$state])) - valueFS_HMSS[busData$state,1]
X = valueFS_HMSS[busData$state,-1]
solve(crossprod(X))%*%crossprod(X,y)

#Now with glm
y = busData$choice
X = valueFS_HMSS[busData$state,-1]
summary(glm(y~X-1,family = "binomial"))
theta

#MEASUREMENT ERROR ---------------------------------------------------------------

source("logistic_with_measurement_error_support_functions.R")
source("logistic_with_measurement_error.R")

#SIMULATIONS TO CONSIDER:
# (A) number of draws in the forward simulation
# (B) number of choice occasions (note this is different from the number of periods of lookahead in the forward simulation)
# (C) The across-choice-occasions variance of the (true) value of the value function
# (D) heavy-tailness of the distribution of the draws -- high skewness and kurtosis. I do not know how you would control this but see if you can
# (E) the straddling of probabilities across the nonlinear region of the logistic function
# 
# Theory from the literature on measurement error and method-of-moments
# indicates that HMSS will do worse than our proposed method when
# (A) (B) and (C) are low and (D) (E) are high.  Because you have
# programmed the HMSS method, you should be able to identify the
# settings of the parameter values and the data sizes that lead to the
# above conditions with respect to (A) through (E) where HMSS does
# relatively poorly and relatively well. This I suggest I work on while
# I complete the coding.

#We need the output from the forward simulation to be an array of size:
#n x p_with_measurement_error x n_replicates

#essentially modified valHMSS,
valueFS_HMSS_ME     = valHMSS_ME(1:S,T=75,R=200,beta,p,pchoice)
valueFS_HMSS_ME_bar = apply(valueFS_HMSS_ME,c(1,2),FUN=mean)

theta

#Original HMSS (assuming known pchoice)
y = log(pchoice[busData$state]/(1-pchoice[busData$state])) - valueFS_HMSS[busData$state,1]
X = valueFS_HMSS[busData$state,-1]
solve(crossprod(X))%*%crossprod(X,y)

#Approach taken in paper
y = busData$choice
X = valueFS_HMSS_ME_bar[busData$state,-1]
summary(glm(y~X-1,family = "binomial"))$coef

X_ME = valueFS_HMSS_ME[busData$state,-1,]

compute_logistic_regr_under_measurement_error(y, # vector of length n
                                              matrix(runif(length(y),-1,1),ncol=1), #matrix of size n x p_without_measurement_error, should include a column of ones if model has an intercept, cannot be a vector even if p_without_measurement_error=1
                                              X_ME, #3D array of size n x p_with_measurement_error x n_replicates, cannot be a matrix even if p_with_measurement_error=1
                                              n_draws_of_mismeasured_x=200,n_cycles_internal_to_arms_gibbs=10,
                                              n_burnin_draws=100,n_reported_draws=3000,method="hbmix", # one of npshrinkage, hbmix, tdist, johnson_cornish_fisher, edgeworth, gram,
                                              scheme="pluginmle") # one of pluginmle, em, mcmc)

#Set up simulations
#structural parameters
theta1  = 2.0   #replacement cost
theta2  = 0.09  #maintenance cost
beta    = 0.9   #discount factor (.9 in paper)
S       = 90    #discretized states (90 in paper)
nbus    = 1     #in this simulation, 1 bus with varying number of choice occasions
nperiod = 5000  #maximum number of choice occasions
theta   = c(theta1,theta2)
T       = 75
#Regardless of method, we need an estimate of p and pchoice, which will be fixed across comparisons and is a function of beta  and theta
p
out_cm     = contraction_mapping(S=S,p=p,MF=lin_cost,params=theta,beta=beta)
pchoice    = out_cm$CP_forward[,2]
busData    = bus_dgp(nbus,nperiod,pchoice,p)

#First, obtain "true" value function estimates (set R to large number)
val_True   = valHMSS_ME_mat(1:S,T=75,R=1000,beta,p,pchoice) #version _mat averages already

#"Truth" to aim for:
y = busData$choice
X = val_True[busData$state,-1]
theta_true = summary(glm(y~X-1,family = "binomial"))$coef
theta_true

#parameters to vary|structural parameters
nobs    = 500   #NUMBER OF CHOICE OCCASIONS TO CONSIDER
R       = 50    #NUMBER OF DRAWS IN THE FORWARD SIMULATION

#busData|number of choice occassions
busData_temp = busData[1:nobs,]
y = busData_temp$choice
X = valHMSS_ME(busData_temp$state,T=75,R=R,beta,p,pchoice)
X_bar = apply(X,c(1,2),FUN=mean)

#Now estimate with three methods (assuming state space is too complex): 
#1) HMSS forward simulation (method of moments)
y_hmss = log(pchoice[busData_temp$state]/(1-pchoice[busData_temp$state])) - X_bar[,1]
X_hmss = X_bar[,-1]
out_hmss = solve(crossprod(X_hmss))%*%crossprod(X_hmss,y_hmss)
out_hmss

#2) Forward simulation with logisitic regression
X_log = X_bar[,-1]
out_log = summary(glm(y~X_log-1,family = "binomial"))
out_log$coef

#3) Forward simulation with logit accounting for missing data
X_me   = X[,-1,]
out_me = compute_logistic_regr_under_measurement_error(y, # vector of length n
                                              matrix(runif(length(y),-1,1),ncol=1), #matrix of size n x p_without_measurement_error, should include a column of ones if model has an intercept, cannot be a vector even if p_without_measurement_error=1
                                              X_me, #3D array of size n x p_with_measurement_error x n_replicates, cannot be a matrix even if p_with_measurement_error=1
                                              n_draws_of_mismeasured_x=200,n_cycles_internal_to_arms_gibbs=10,
                                              n_burnin_draws=100,n_reported_draws=300,method="hbmix", # one of npshrinkage, hbmix, tdist, johnson_cornish_fisher, edgeworth, gram,
                                              scheme="pluginmle") # one of pluginmle, em, mcmc)
out_me$coef











#HMSS 1994: Estimate p(choice|state) and transition probabilities ---------------------

#In reality, we don't know p (state transition probabilities) or pchoice (probability of replace|state)
#So we estimate these

#Equation 3.1 for estimating transition probabilities p
library(dplyr)

busData = busData %>% 
  mutate(rowNumber = 1:n()) %>%
  arrange(id,rowNumber) %>%
  mutate(nextState  = lead(state,1),
         nextId     = lead(id,1),
         deltaState = nextState - state)


phat =  busData %>%
  filter(id == nextId & choice == 0) %>% 
  group_by(deltaState) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n)) %>% .$freq

phat
p

#Equation 3.3 for estimating conditional choice probabilities pchoice (replacement frequencies)
pchoicehat_prep = busData %>%
  group_by(state) %>% 
  summarise(n = n(),choice = sum(choice)) %>%
  mutate(pchoicehat = choice/n)

pchoicehat_rf = rep(0,S)
pchoicehat_rf[pchoicehat_prep$state] = pchoicehat_prep$pchoicehat
#remove 0 and 1 for log(y/(1-y)) estimate
pchoicehat_rf[pchoicehat_rf == 0] = 1e-5
pchoicehat_rf[pchoicehat_rf == 1] = 1 - 1e-5

#VARIATIONS TO ESTIMATE PCHOICEHAT -------------------------------
#use true pchoice
pchoicehat_true = pchoice

#use polynomial estimate
pchoicehat_poly = glm(choice~poly(state,5),busData,family="binomial")
pchoicehat_poly = predict(pchoicehat_poly,newdata = data.frame(state = 1:S),type="response")

#use kernel estimator
kernelFunction = function(state,i,xi) dnorm((state - i)/xi)

getPchoiceKernel = function(xi,tol = 1e-5){
  pchoicehat = rep(NA,S)
  for(s in 1:S){
    pchoicehat[s] = sum((kernelFunction(busData$state,s,xi) * busData$choice)/sum(kernelFunction(busData$state,s,xi)),na.rm=TRUE)
  }
  
  #Make sure no probabilities are 0 or 1
  pchoicehat[pchoicehat == 0] = tol
  pchoicehat[pchoicehat == 1] = 1 - tol
  
  pchoicehat
}

pchoice_kernel.0250 = getPchoiceKernel(.025)
pchoice_kernel.0100 = getPchoiceKernel(.01)
pchoice_kernel.0050 = getPchoiceKernel(.005)
pchoice_kernel.0025 = getPchoiceKernel(.0025)
pchoice_kernel.0010 = getPchoiceKernel(.001)

#COMPARE PCHOICEHAT METHODS -----------------
library(ggplot2)





#SIMULATIONS --------------------------------
stateCount   = aggregate(id~state,busData,FUN = length)
stateCount$n = stateCount$id

getHMSSthetahat = function(T,R,phat,pchoicehat,CoxCorrection = FALSE){
  #same outcome as weighting function
  outHMSS = valHMSS(S,T,R,beta,phat,pchoicehat)
  if(CoxCorrection){
    y = log((pchoicehat[busData$state] + 1/stateCount$n[busData$state])/
              (1 - pchoicehat[busData$state] + 1/stateCount$n[busData$state])) - outHMSS[busData$state,1]
  } else {
    y = log(pchoicehat[busData$state]/(1-pchoicehat[busData$state])) - outHMSS[busData$state,1]  
  }
  X = outHMSS[busData$state,-1]
  solve(crossprod(X))%*%crossprod(X,y)
}


sims     = 5
simArray = array(dim=c(sims,2,9)) #2 variables x 9 variations
T = 50
R = 50
for(sim in 1:sims){
  
  #MLE
  # fit = optim(par=rep(.1,npars),fn=dynamicLogit,method=c("L-BFGS-B"),lower=bounds[1],upper=bounds[2],
  #                 data=busData,S=S,p=p,MF=lin_cost,beta=beta,control=list(fnscale=1),hessian = FALSE)
  # simArray[sim,,2] = fit$par
  
  #Replacement frequencies
  simArray[sim,,2] = getHMSSthetahat(T,R,phat,pchoicehat_rf)
  
  #True pchoice
  simArray[sim,,3] = getHMSSthetahat(T,R,phat,pchoicehat_true)
  
  #Cox Correction
  simArray[sim,,4] = getHMSSthetahat(T,R,phat,pchoicehat_rf,CoxCorrection = TRUE)
  
  #Kernel estimates
  simArray[sim,,5] = getHMSSthetahat(T,R,phat,pchoice_kernel.0250)
  simArray[sim,,6] = getHMSSthetahat(T,R,phat,pchoice_kernel.0100)
  simArray[sim,,7] = getHMSSthetahat(T,R,phat,pchoice_kernel.0050)
  simArray[sim,,8] = getHMSSthetahat(T,R,phat,pchoice_kernel.0025)
  simArray[sim,,9] = getHMSSthetahat(T,R,phat,pchoice_kernel.0010)

  print(sim)
}

t(apply(simArray,c(2,3),FUN=mean))
t(apply(simArray,c(2,3),FUN=sd))

#from CCP
fit_CCP$par
