// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube valHMSS_ME(vec states, int T, int R, double beta, vec phat, vec pchoicehat) {
  //notice we have a vector of states, not just S
  //also: output is an array of draws, not a matrix of averages
  
  int S = states.size();
  int maxS = 90;
  
  //Simulate value function from t = 0 forward
  vec betaVec = zeros<vec>(T + 1);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  mat  outR      = zeros<mat>(S,3); //difference in discounted values for both actions
  cube out       = zeros<cube>(S,3,R); //difference in discounted values for both actions
  vec pStates  = {0,1,2};         //transition states associated with transition probabilities p
  
  double gamma_cnst = .5772;
  
  vec state    = zeros<vec>(T + 1);
  vec choice   = zeros<vec>(T + 1);
  vec pchoice  = zeros<vec>(T + 1); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T + 1,2);
  mat choiceMat     = zeros<mat>(T + 1,2);
  mat errorMat      = zeros<mat>(T + 1,2); //store the expected error terms
  
  double x0,x1,x2;
  //In the paper, states start at i = 1...S (not 0)
  
  //NOTE: R AND S LOOPS HAVE BEEN SWITCHED
  //for R repititions
  for(int r = 0; r < R; r++){
    
    //Loop over each state (obsservation)
    for(int s = 0; s < S; s++){
      
      //for each action: 0 = maintain, 1 = replace
      for(int a = 0; a < 2; a++){
        
        //reset "state" and "maintain" vectors
        // state[0]       = s + 1; //states are from 1...90
        state[0]       = states[s]; //states are from 1...90
        choice[0]      = a;     //when a = 0, maintain, else = 1
        
        //forward simulate T time periods (first time period pre-determined by action a)
        for(int t = 1; t < (T+1); t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          }
          
          //Ensure maximum state is bounded
          if(state[t] > maxS) state[t] = maxS;
          
          //Based on state, make decision for next period
          //Substract 1 for c++ indexing
          pchoice[t]  = pchoicehat[state[t]-1];
          choice[t]   = (pchoice[t] > runif(1)[0])*1;
          
        } //end of T forward simulations
        
        //fill in state, choice, and pchocie matrixes based on action a
        //a = 1 means REPLACE
        stateMat.col(a)    = state;
        choiceMat.col(a)   = choice;
        //If the decision is to replace, use gamma + log(pchoice), otherwise use gamma + log(1 - pchoice)
        errorMat.col(a)    = gamma_cnst - log(choice % pchoice + (1 - choice) % (1 - pchoice));
        
      } //end of action a loop
      
      //sum the discounted utilities
      //t0 = maintain
      //a0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0)%(1-choiceMat.col(0)) + gamma_cnst - log(errorMat.col(0))));
      
      //t0 = replace
      //a1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1)%(1-choiceMat.col(1)) + gamma_cnst - log(errorMat.col(1))));
      
      //Need a1 - a0
      
      //Note: in time = 0 there is no error term
      errorMat.row(0).zeros();
      
      x0 = sum(betaVec % (errorMat.col(1) - errorMat.col(0)));
      
      x1 = sum(betaVec % (-choiceMat.col(1) - -choiceMat.col(0)));
      
      x2 = sum(betaVec % (-stateMat.col(1) % (1-choiceMat.col(1)) - -stateMat.col(0) % (1-choiceMat.col(0)))); 
      
      //place into outEV
      outR(s,0) = x0;
      outR(s,1) = x1;
      outR(s,2) = x2;
      
    } //end of S states loop
    
    out.slice(r) = outR;
    
  } //end of R simulations loop
  
  // //average over simulations
  // out = out/double(R);
  
  return out;
}



// [[Rcpp::export]]
mat valHMSS_ME_mat(vec states, int T, int R, double beta, vec phat, vec pchoicehat) {
  //notice we have a vector of states, not just S
  
  int S = states.size();
  int maxS = 90;
  
  //Simulate value function from t = 0 forward
  vec betaVec = zeros<vec>(T + 1);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  mat out      = zeros<mat>(S,3); //difference in discounted values for both actions
  vec pStates  = {0,1,2};         //transition states associated with transition probabilities p
  
  double gamma_cnst = .5772;
  
  vec state    = zeros<vec>(T + 1);
  vec choice   = zeros<vec>(T + 1);
  vec pchoice  = zeros<vec>(T + 1); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T + 1,2);
  mat choiceMat     = zeros<mat>(T + 1,2);
  mat errorMat      = zeros<mat>(T + 1,2); //store the expected error terms
  
  double x0,x1,x2;
  //In the paper, states start at i = 1...S (not 0)
  
  //Loop over each state (obsservation)
  for(int s = 0; s < S; s++){
    
    //for R repititions
    for(int r = 0; r < R; r++){
      
      //for each action: 0 = maintain, 1 = replace
      for(int a = 0; a < 2; a++){
        
        //reset "state" and "maintain" vectors
        // state[0]       = s + 1; //states are from 1...90
        state[0]       = states[s]; //states are from 1...90
        choice[0]      = a;     //when a = 0, maintain, else = 1
        
        //forward simulate T time periods (first time period pre-determined by action a)
        for(int t = 1; t < (T+1); t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          }
          
          //Ensure maximum state is bounded
          if(state[t] > maxS) state[t] = maxS;
          
          //Based on state, make decision for next period
          //Substract 1 for c++ indexing
          pchoice[t]  = pchoicehat[state[t]-1];
          choice[t]   = (pchoice[t] > runif(1)[0])*1;
          
        } //end of T forward simulations
        
        //fill in state, choice, and pchocie matrixes based on action a
        //a = 1 means REPLACE
        stateMat.col(a)    = state;
        choiceMat.col(a)   = choice;
        //If the decision is to replace, use gamma + log(pchoice), otherwise use gamma + log(1 - pchoice)
        errorMat.col(a)    = gamma_cnst - log(choice % pchoice + (1 - choice) % (1 - pchoice));
        
      } //end of action a loop
      
      //sum the discounted utilities
      //t0 = maintain
      //a0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0)%(1-choiceMat.col(0)) + gamma_cnst - log(errorMat.col(0))));
      
      //t0 = replace
      //a1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1)%(1-choiceMat.col(1)) + gamma_cnst - log(errorMat.col(1))));
      
      //Need a1 - a0
      
      //Note: in time = 0 there is no error term
      errorMat.row(0).zeros();
      
      x0 = sum(betaVec % (errorMat.col(1) - errorMat.col(0)));
      
      x1 = sum(betaVec % (-choiceMat.col(1) - -choiceMat.col(0)));
      
      x2 = sum(betaVec % (-stateMat.col(1) % (1-choiceMat.col(1)) - -stateMat.col(0) % (1-choiceMat.col(0)))); 
      
      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      out(s,2) = out(s,2) + x2;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over simulations
  out = out/double(R);
  
  return out;
}

// [[Rcpp::export]]
mat valHMSS(int S, int T, int R, double beta, vec phat, vec pchoicehat) {
  
  //Simulate value function from t = 0 forward
  vec betaVec = zeros<vec>(T + 1);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  mat out      = zeros<mat>(S,3); //difference in discounted values for both actions
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  
  double gamma_cnst = .5772;
  
  vec state    = zeros<vec>(T + 1);
  vec choice   = zeros<vec>(T + 1);
  vec pchoice  = zeros<vec>(T + 1); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T + 1,2);
  mat choiceMat     = zeros<mat>(T + 1,2);
  mat errorMat      = zeros<mat>(T + 1,2); //store the expected error terms
  
  double x0,x1,x2;
  //In the paper, states start at i = 1...S (not 0)
  
  //Loop over each state
  for(int s = 0; s < S; s++){
    
    //for R repititions
    for(int r = 0; r < R; r++){
      
      //for each action: 0 = maintain, 1 = replace
      for(int a = 0; a < 2; a++){
        
        //reset "state" and "maintain" vectors
        state[0]       = s + 1; //states are from 1...90
        choice[0]      = a;     //when a = 0, maintain, else = 1
        
        //forward simulate T time periods (first time period pre-determined by action a)
        for(int t = 1; t < (T+1); t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          }
          
          //Ensure maximum state is bounded
          if(state[t] > S) state[t] = S;
          
          //Based on state, make decision for next period
          //Substract 1 for c++ indexing
          pchoice[t]  = pchoicehat[state[t]-1];
          choice[t]   = (pchoice[t] > runif(1)[0])*1;
          
        } //end of T forward simulations
        
        //fill in state, choice, and pchocie matrixes based on action a
        //a = 1 means REPLACE
        stateMat.col(a)    = state;
        choiceMat.col(a)   = choice;
        //If the decision is to replace, use gamma + log(pchoice), otherwise use gamma + log(1 - pchoice)
        errorMat.col(a)    = gamma_cnst - log(choice % pchoice + (1 - choice) % (1 - pchoice));
        
      } //end of action a loop
      
      //sum the discounted utilities
      //t0 = maintain
      //a0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0)%(1-choiceMat.col(0)) + gamma_cnst - log(errorMat.col(0))));
      
      //t0 = replace
      //a1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1)%(1-choiceMat.col(1)) + gamma_cnst - log(errorMat.col(1))));
      
      //Need a1 - a0
      
      //Note: in time = 0 there is no error term
      errorMat.row(0).zeros();
      
      x0 = sum(betaVec % (errorMat.col(1) - errorMat.col(0)));
      
      x1 = sum(betaVec % (-choiceMat.col(1) - -choiceMat.col(0)));
      
      x2 = sum(betaVec % (-stateMat.col(1) % (1-choiceMat.col(1)) - -stateMat.col(0) % (1-choiceMat.col(0)))); 
      
      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      out(s,2) = out(s,2) + x2;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over simulations
  out = out/double(R);
  
  return out;
}

// [[Rcpp::export]]
mat futureValHMSS(int S, int T, int R, double beta, vec phat, vec pchoicehat) {
  
  //Simulate value function from t+1 forward
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t + 1);
  
  mat out      = zeros<mat>(S,3); //difference in discounted values for both actions
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  
  double gamma_cnst = .5772;
  
  vec state    = zeros<vec>(T + 1);
  vec choice   = zeros<vec>(T + 1);
  vec pchoice  = zeros<vec>(T + 1); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  mat errorMat      = zeros<mat>(T,2); //store the expected error terms
  
  double x0,x1,x2;
  //In the paper, states start at i = 1...S (not 0)
  
  //Loop over each state
  for(int s = 0; s < S; s++){
    
    //for R repititions
    for(int r = 0; r < R; r++){
      
      //for each action: 0 = maintain, 1 = replace
      for(int a = 0; a < 2; a++){
        
        //reset "state" and "maintain" vectors
        state[0]       = s + 1; //states are from 1...90
        choice[0]      = a;     //when a = 0, maintain, else = 1
        pchoice[0]     = pchoicehat[s];
        
        //forward simulate T time periods (first time period pre-determined by action a)
        for(int t = 1; t < (T+1); t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          }
          
          //Ensure maximum state is bounded
          if(state[t] > S) state[t] = S;
          
          //Based on state, make decision for next period
          //Substract 1 for c++ indexing
          pchoice[t]  = pchoicehat[state[t]-1];
          choice[t]   = (pchoice[t] > runif(1)[0])*1;
          
        } //end of T forward simulations
        
        //fill in state, choice, and pchocie matrixes based on action a
        //a = 1 means REPLACE
        stateMat.col(a)    = state.tail(T);
        choiceMat.col(a)   = choice.tail(T);
        //If the decision is to replace, use gamma + log(pchoice), otherwise use gamma + log(1 - pchoice)
        errorMat.col(a)    = gamma_cnst - log(choice.tail(T) % pchoice.tail(T) + (1 - choice.tail(T)) % (1 - pchoice.tail(T)));
        
      } //end of action a loop
      
      //sum the discounted utilities
      //t0 = maintain
      //a0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0)%(1-choiceMat.col(0)) + errorMat.col(0)));
      
      //t0 = replace
      //a1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1)%(1-choiceMat.col(1)) + errorMat.col(1)));
      
      //Need a1 - a0
      x0 = sum(betaVec % (errorMat.col(1) - errorMat.col(0)));
      
      x1 = sum(betaVec % (-choiceMat.col(1) - -choiceMat.col(0)));
      
      x2 = sum(betaVec % (-stateMat.col(1) % (1-choiceMat.col(1)) - -stateMat.col(0) % (1-choiceMat.col(0)))); 
      
      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      out(s,2) = out(s,2) + x2;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/double(R);
  
  return out;
}

// [[Rcpp::export]]
mat futureValFS(int S, int T, int R, double beta, vec phat, vec pchoicehat, vec thetahat) {
  
  //Simulate value function from t+1 forward
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t + 1);
  
  mat out      = zeros<mat>(S,2); //difference in discounted values for both actions
  vec pStates  = {0,1,2};         //transition states associated with transition probabilities p
  
  double gamma_cnst = .5772;
  
  vec state    = zeros<vec>(T + 1);
  vec choice   = zeros<vec>(T + 1);
  vec pchoice  = zeros<vec>(T + 1); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  mat errorMat      = zeros<mat>(T,2); //store the expected error terms
  
  double x0,x1;
  //In the paper, states start at i = 1...S (not 0)
  
  //Loop over each state
  for(int s = 0; s < S; s++){
    
    //for R repititions
    for(int r = 0; r < R; r++){
      
      //for each action: 0 = maintain, 1 = replace
      for(int a = 0; a < 2; a++){
        
        //reset "state" and "maintain" vectors
        state[0]       = s + 1; //states are from 1...90
        choice[0]      = a;     //when a = 0, maintain, else = 1
        
        //forward simulate T time periods (first time period pre-determined by action a)
        for(int t = 1; t < (T+1); t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          }
          
          //Ensure maximum state is bounded
          if(state[t] > S) state[t] = S;
          
          //Based on state, make decision for next period
          //Substract 1 for c++ indexing
          pchoice[t]  = pchoicehat[state[t]-1];
          choice[t]   = (pchoice[t] > runif(1)[0])*1;
          
        } //end of T forward simulations
        
        //fill in state, choice, and pchocie matrixes based on action a
        //a = 1 means REPLACE
        stateMat.col(a)    = state.tail(T);
        choiceMat.col(a)   = choice.tail(T);
        //If the decision is to replace, use gamma + log(pchoice), otherwise use gamma + log(1 - pchoice)
        errorMat.col(a)    = gamma_cnst - log(choice.tail(T) % pchoice.tail(T) + (1 - choice.tail(T)) % (1 - pchoice.tail(T)));
        
      } //end of action a loop
      
      //sum the discounted utilities
      //t0 = maintain
      x0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0)%(1-choiceMat.col(0)) + errorMat.col(0)));
      
      //t0 = replace
      x1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1)%(1-choiceMat.col(1)) + errorMat.col(1)));
      
      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/double(R);
  
  return out;
}