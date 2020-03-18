library(RcppArmadillo)
library(dplyr)
library(e1071)
#see https://github.com/QuentinAndre/John-Rust-1987-Python/blob/master/Data%20Generation%20and%20Likelihood%20Fit.ipynb
#Goal: replicate Monte Carlo Simulation section of HMSS 1994

#SUPPORT FUNCTIONS -----------------------------------------------------
Rcpp::sourceCpp("HMSSRustwithME.cpp")

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
  while(i <= nobs){
    prev_ind = seq(from=i,by=1,length.out=nbus) #index of the "last" decisions
    next_ind = prev_ind+nbus                    #index of the "next" decisions
    
    prev_states = bus_df$state[prev_ind]  
    new_states  = prev_states
    
    # Generating choices, where 1 = replace, from choice probabilities, conditional on the state s of each bus
    choices = (pchoice[prev_states] > runif(nbus))*1 
    bus_df$choice[prev_ind] = choices
    
    if(tail(next_ind)[1] <= nobs){
      
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

#Set up function to conduct simulations
runSimulations = function(theta1 = 2,theta2 = .09,nsims = 100,nobs = 100, R = 50, print = FALSE) {
  
  #nsims: number of simulations to conduct (except for "truth")
  #nobs: number of data observations used for estimation
  #R: number of forward simulations to average over
  out_parameters = list(theta1 = theta1, theta2 = theta2, nsims = nsims, nobs = nobs, R = R)
  
  #Structural parameters
  p       = c(.349,.639,.012) #transition probabilities
  beta    = 0.9   #discount factor (.9 in paper)
  S       = 90    #discretized states (90 in paper)
  nbus    = 1     #in this simulation, 1 bus with varying number of choice occasions
  nperiod = 2000  #maximum number of choice occasions
  theta   = c(theta1,theta2)
  T       = 75
  #Regardless of method, we need an estimate of p and pchoice, which will be fixed across comparisons and is a function of beta  and theta
  out_cm     = contraction_mapping(S=S,p=p,MF=lin_cost,params=theta,beta=beta)
  pchoice    = out_cm$CP_forward[,2]
  busData    = bus_dgp(nbus,nperiod,pchoice,p)
  
  #True value function and coefficient estimates
  
  #First, obtain "true" value function estimates (set R to large number)
  val_True   = valHMSS_ME_mat(1:S,T=75,R=5000,beta,p,pchoice) #version _mat averages already
  
  #Estimate with either logistic regression or with method of moments
  y = busData$choice
  X = val_True[busData$state,-1]
  true_lr = summary(glm(y~X-1,family = "binomial"))
  theta_lr       = true_lr$coef[,1]
  thetaSE_lr     = true_lr$coef[,2]
  
  y = log(pchoice[busData$state]/(1-pchoice[busData$state])) - val_True[busData$state,1]
  X = val_True[busData$state,-1]
  theta_true_mom = solve(crossprod(X))%*%crossprod(X,y)
  
  #Properties of X
  X = val_True[busData$state,-1]
  X_variance  = apply(X,2,FUN=var)
  X_skewness  = apply(X,2,FUN=skewness)
  X_kurtosis  = apply(X,2,FUN=kurtosis)
  
  out_truth = list(val_True = val_True, 
                  theta_lr   = as.numeric(theta_lr),
                  thetaSE_lr = as.numeric(thetaSE_lr),
                  theta_true_mom  = as.numeric(theta_true_mom),
                  X_variance = X_variance,
                  X_skewness = X_skewness,
                  X_kurtosis = X_kurtosis)
  
  out_simulations = data.frame(sim        = 1:nsims,
                               hmssTheta1   = NA,
                               hmssTheta2   = NA,
                               logRegTheta1 = NA,
                               logRegTheta2 = NA,
                               logRegTheta1SE = NA,
                               logRegTheta2SE = NA,
                               logRegMETheta1 = NA,
                               logRegMETheta2 = NA,
                               logRegMETheta1SE = NA,
                               logRegMETheta2SE = NA)
  
  #Now proceed with nsims simulations
  for(sim in 1:nsims){
    
    #take subsample of busData based on nobs (note does not replace)
    busInd = sample(1:nrow(busData),nobs)
    busData_temp = busData[busInd,]
    
    #Use esimate of pchoice, not true pchoice
    pchoicehat_prep = busData_temp %>% 
      mutate(rowNumber = 1:n()) %>%
      arrange(id,rowNumber) %>%
      mutate(nextState  = lead(state,1),
             nextId     = lead(id,1),
             deltaState = nextState - state) %>%
      group_by(state) %>% 
      summarise(n = n(),choice = sum(choice)) %>%
      mutate(pchoicehat = choice/n)
    
    pchoicehat_rf = rep(0,S)
    pchoicehat_rf[pchoicehat_prep$state] = pchoicehat_prep$pchoicehat
    #remove 0 and 1 for log(y/(1-y)) estimate
    pchoicehat_rf[pchoicehat_rf == 0] = 1e-5
    pchoicehat_rf[pchoicehat_rf == 1] = 1 - 1e-5
    pchoice = pchoicehat_rf
    

    y = busData_temp$choice
    X = valHMSS_ME(busData_temp$state,T=75,R=R,beta,p,pchoice)
    X_bar = apply(X,c(1,2),FUN=mean)
    
    #Now estimate with three methods (assuming state space is too complex): 
    #1) HMSS forward simulation (method of moments)
    y_hmss = log(pchoice[busData_temp$state]/(1-pchoice[busData_temp$state])) - X_bar[,1]
    X_hmss = X_bar[,-1]
    out_hmss = solve(crossprod(X_hmss))%*%crossprod(X_hmss,y_hmss)
    
    #2) Forward simulation with logistic regression
    X_log = X_bar[,-1]
    out_log = summary(glm(y~X_log-1,family = "binomial"))
  
    #3) Forward simulation with logit accounting for missing data
    X_me   = X[,-1,]
    out_me = compute_logistic_regr_under_measurement_error(y, # vector of length n
                                                           matrix(runif(length(y),-1,1),ncol=1), #matrix of size n x p_without_measurement_error, should include a column of ones if model has an intercept, cannot be a vector even if p_without_measurement_error=1
                                                           X_me, #3D array of size n x p_with_measurement_error x n_replicates, cannot be a matrix even if p_with_measurement_error=1
                                                           n_draws_of_mismeasured_x=200,n_cycles_internal_to_arms_gibbs=10,
                                                           n_burnin_draws=100,n_reported_draws=300,method="hbmix", # one of npshrinkage, hbmix, tdist, johnson_cornish_fisher, edgeworth, gram,
                                                           scheme="pluginmle") # one of pluginmle, em, mcmc
    
    #Save output
    out_simulations[sim,2:3] = out_hmss
    
    out_simulations[sim,4:5] = out_log$coef[,1]
    out_simulations[sim,6:7] = out_log$coef[,2]
    
    out_simulations[sim,8:9]   = out_me$coef[-1]
    out_simulations[sim,10:11] = out_me$stderr[-1]
    
    #print
    if(print & sim %% 10 == 0) print(sim)
    
  }
  
  out = list(parameters  = out_parameters,
             truth       = out_truth,
             simulations = out_simulations)
  
  out
  
}

#RUN SIMULATIONS --------------------------------------------------------------------------------

out1  = runSimulations(theta1 = 2.0, theta2 = 0.09, nobs = 10, R = 50);save(out1,file="simulations/out1.Rdata")
out2  = runSimulations(theta1 = 2.0, theta2 = 0.09, nobs = 50, R = 50);save(out2,file="simulations/out2.Rdata")
out3  = runSimulations(theta1 = 2.0, theta2 = 0.09, nobs = 100, R = 50);save(out3,file="simulations/out3.Rdata")
out4  = runSimulations(theta1 = 2.0, theta2 = 0.09, nobs = 1000, R = 50);save(out4,file="simulations/out4.Rdata")
out5  = runSimulations(theta1 = 2.0, theta2 = 0.09, nobs = 10, R = 50);save(out5,file="simulations/out5.Rdata")
out6  = runSimulations(theta1 = 2.0, theta2 = 0.09, nobs = 10, R = 100);save(out6,file="simulations/out6.Rdata")
out7  = runSimulations(theta1 = 2.0, theta2 = 0.09, nobs = 10, R = 200);save(out7,file="simulations/out7.Rdata")
out8  = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 10, R = 50);save(out8,file="simulations/out8.Rdata")
out9  = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 50, R = 50);save(out9,file="simulations/out9.Rdata")
out10 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 100, R = 50);save(out10,file="simulations/out10.Rdata")
out11 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 1000, R = 50);save(out11,file="simulations/out11.Rdata")
out12 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 10, R = 50);save(out12,file="simulations/out12.Rdata")
out13 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 10, R = 100);save(out13,file="simulations/out13.Rdata")
out14 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 10, R = 200);save(out14,file="simulations/out14.Rdata")
out15 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 10, R = 50);save(out15,file="simulations/out15.Rdata")
out16 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 50, R = 50);save(out16,file="simulations/out16.Rdata")
out17 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 100, R = 50);save(out17,file="simulations/out17.Rdata")
out18 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 1000, R = 50);save(out18,file="simulations/out18.Rdata")
out19 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 10, R = 50);save(out19,file="simulations/out19.Rdata")
out20 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 10, R = 100);save(out20,file="simulations/out20.Rdata")
out21 = runSimulations(theta1 = 5.0, theta2 = 0.09, nobs = 10, R = 200);save(out21,file="simulations/out21.Rdata")

#Summarize all output into a table
scenarios = 21

#first create list of all output
all_out = list()
for(s in 1:scenarios){
  file_out = paste("out",s,sep = "")
  load(paste("simulations/",file_out,".rdata",sep=""))
  all_out[[s]] = get(file_out)
  rm(list = ls(pattern=("out."))) #can be done outside the loop, or remove as you go
}

#Set up summary table
summaryTable = data.frame(scenario      = 1:scenarios,
                          theta1        = NA,
                          theta2        = NA,
                          replications  = NA,
                          observations  = NA,
                          forwardsims   = NA,
                          thetahat1_lr   = NA,
                          thetahat2_lr   = NA,
                          thetahat1_odds = NA,
                          thetahat2_odds = NA,
                          X1_variance   = NA,
                          X2_variance   = NA,
                          X1_skewness   = NA,
                          X2_skewness   = NA,
                          X1_kurtosis   = NA,
                          X2_kurtosis   = NA,
                          hmssT1    = NA,
                          hmssT2    = NA,
                          hmssT1empSE   = NA,
                          hmssT2empSE   = NA,
                          logT1 = NA,
                          logT2 = NA,
                          logT1SE = NA,
                          logT2SE = NA,
                          logT1empSE = NA,
                          logT2empSE = NA,
                          logMET1 = NA,
                          logMET2 = NA,
                          logMET1SE = NA,
                          logMET2SE = NA,
                          logT1MEempSE = NA,
                          logT2MEempSE = NA)


#fill in summary table
for(s in 1:length(all_out)){
  
  out = all_out[[s]]
  summaryTable[s,2:6]  = unlist(out$parameters)
  summaryTable[s,7:16] = unlist(out$truth[-c(1,3)])
  
  #mean estimates
  summaryTable[s,c(17:18,21:24,27:30)] = colMeans(out$simulations)[-1]
  
  #empirical standard error
  summaryTable[s,c(19:20,25:26,31:32)] = apply(out$simulations,2,FUN=sd)[c(2,3,4,5,8,9)]
  
}

write.csv(summaryTable,file="simulations/summaryTable.csv",row.names = FALSE)




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
