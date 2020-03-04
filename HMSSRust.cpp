// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat simpleFSgamma(int S, int T, int R, double beta, vec phat, vec pchoicehat, vec thetahat) {
  
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  mat out      = zeros<mat>(S,2); //difference in discounted values for both actions
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  
  double gamma_cnst = .5772;
  
  vec state    = zeros<vec>(T);
  vec choice   = zeros<vec>(T);
  vec pchoice  = zeros<vec>(T); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  mat pchoiceMat     = zeros<mat>(T,2);
  
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
        pchoice[0]     = pchoicehat[s];
        
        //forward simulate T time periods (first time period pre-determined by action a)
        for(int t = 1; t < T; t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
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
        stateMat.col(a)   = state;
        choiceMat.col(a)  = choice;
        pchoiceMat.col(a)  = pchoice;
        
      } //end of action a loop
      
      //sum the discounted utilities
      
      //t0 = maintain
      x0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0)%(1-choiceMat.col(0)) + gamma_cnst - log(pchoiceMat.col(0))));
      
      //t0 = replace
      x1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1)%(1-choiceMat.col(1)) + gamma_cnst - log(pchoiceMat.col(1))));
      
      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/double(R);
  
  return out;
}


// [[Rcpp::export]]
mat simpleHMSS(int S, int T, int R, double beta, vec phat, vec pchoicehat) {
  
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  mat out      = zeros<mat>(S,2); //difference in discounted values for both actions
  vec pStates  = {0,1,2};         //transition states associated with transition probabilities p
  
  vec state    = zeros<vec>(T);
  vec choice   = zeros<vec>(T);
  vec pchoice  = zeros<vec>(T); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  
  double x1,x2;
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
        for(int t = 1; t < T; t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
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
        stateMat.col(a)   = state;
        choiceMat.col(a)  = choice;

      } //end of action a loop
      
      //sum the discounted utilities
      
      //Get a1 - a0 independent of thetahat
      
      //t0 = maintain
      //a0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0)));
      
      //t0 = replace
      //a1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1)));
      
      //more transparent: x1 = sum(betaVec % (-choiceMat.col(1) - -choiceMat.col(0)));
      x1 = sum(betaVec % (choiceMat.col(0) - choiceMat.col(1)));
      
      //more transparent: x2 = sum(betaVec % (-stateMat.col(1) - -stateMat.col(0)));
      x2 = sum(betaVec % (stateMat.col(0) % (1 - choiceMat.col(0)) - stateMat.col(1) % (1 - choiceMat.col(1))));
      
      //place into outEV
      out(s,0) = out(s,0) + x1;
      out(s,1) = out(s,1) + x2;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/double(R);
  
  return out;
}

// [[Rcpp::export]]
mat simpleFS(int S, int T, int R, double beta, vec phat, vec pchoicehat, vec thetahat) {
  
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  mat out      = zeros<mat>(S,2); //difference in discounted values for both actions
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  
  vec state    = zeros<vec>(T);
  vec choice   = zeros<vec>(T);
  vec pchoice  = zeros<vec>(T); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  
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
        pchoice[0]     = pchoicehat[s];
        
        //forward simulate T time periods (first time period pre-determined by action a)
        for(int t = 1; t < T; t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
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
        stateMat.col(a)   = state;
        choiceMat.col(a)  = choice;
        
      } //end of action a loop
      
      //sum the discounted utilities
      
      //t0 = maintain
      x0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0)%(1-choiceMat.col(0))));

      //t0 = replace
      x1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1)%(1-choiceMat.col(1))));

      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/double(R);
  
  return out;
}

// [[Rcpp::export]]
mat simpleHMSS_1994(int S, int T, int R, double beta, vec phat, vec pchoicehat) {
  
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  double gamma_cnst = .5772;
  
  mat out      = zeros<mat>(S,3); //difference in discounted values for both actions
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  
  vec state    = zeros<vec>(T);
  vec choice   = zeros<vec>(T);
  vec pchoice  = zeros<vec>(T); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  mat pchoiceMat    = zeros<mat>(T,2);
  
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
        for(int t = 1; t < T; t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = int(std::min(1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0],S + 0.0));
          } else{
            state[t] = int(std::min(state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0],S + 0.0));
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
        stateMat.col(a)   = state;
        choiceMat.col(a)  = choice;
        pchoiceMat.col(a) = pchoice;
        
      } //end of action a loop
      
      //sum the discounted utilities
      
      //Get a1 - a0 independent of thetahat
      
      //t0 = maintain
      //a0 = -sum(betaVec % (1.0 - pchoiceMat.col(0)) % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0) + gamma_cnst - log(1.0 - pchoiceMat.col(0))));
      
      //t0 = replace
      //a1 = -sum(betaVec % pchoiceMat.col(1) % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1) + gamma_cnst - log(pchoiceMat.col(1))));
      
      
      x0 = -sum(betaVec % (pchoiceMat.col(1)%(gamma_cnst - log(pchoiceMat.col(1))) - 
        (1.0-pchoiceMat.col(0)) % (gamma_cnst - log(1.0-pchoiceMat.col(0)))));
      
      x1 = -sum(betaVec % (pchoiceMat.col(1) % -choiceMat.col(1) - (1.0 - pchoiceMat.col(0)) % -choiceMat.col(0)));
      
      x2 = -sum(betaVec % (pchoiceMat.col(1) % -stateMat.col(1) - (1.0-pchoiceMat.col(0)) % -stateMat.col(0)));
      
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
mat simpleFS_1994(int S, int T, int R, double beta, vec phat, vec pchoicehat, vec thetahat) {
  
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  double gamma_cnst = .5772;
  
  mat out      = zeros<mat>(S,2); //difference in discounted values for both actions
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  
  vec state    = zeros<vec>(T);
  vec choice   = zeros<vec>(T);
  vec pchoice  = zeros<vec>(T); //probability of choosing each state|simulated states
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  mat pchoiceMat    = zeros<mat>(T,2);
  
  double x0,x1;
  //In the paper, states start at i = 1...S (not 0)
  
  vec temp;
  
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
        for(int t = 1; t < T; t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
          } else{
            state[t] = state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0];
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
        stateMat.col(a)   = state;
        choiceMat.col(a)  = choice;
        pchoiceMat.col(a) = pchoice;
        
      } //end of action a loop
      
      //sum the discounted utilities
      
      //t0 = maintain
      //x0 = sum(betaVec % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0) + gamma_cnst - log(pchoiceMat.col(0))));
      x0 = -sum(betaVec % (1.0 - pchoiceMat.col(0)) % (-thetahat[0]*choiceMat.col(0) + -thetahat[1]*stateMat.col(0) + gamma_cnst - log(1.0 - pchoiceMat.col(0))));
      
      //t0 = replace
      //x1 = sum(betaVec % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1) + gamma_cnst - log(pchoiceMat.col(1))));
      x1 = -sum(betaVec % pchoiceMat.col(1) % (-thetahat[0]*choiceMat.col(1) + -thetahat[1]*stateMat.col(1) + gamma_cnst - log(pchoiceMat.col(1))));
      
      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/double(R);
  
  return out;
}


// [[Rcpp::export]]
mat getFS(int S, int T, int R, double beta, vec phat, vec pchoicehat, vec thetahat) {
  
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t + 1.0);
  
  double gamma_cnst = .5772;
  
  mat out      = zeros<mat>(S,2); //difference in discounted values for both actions
  out = out + 0.0;
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  vec onesT    = ones<vec>(T) + 0.0;
  
  vec state    = zeros<vec>(T + 1);
  vec choice   = zeros<vec>(T + 1);
  vec pchoice  = zeros<vec>(T + 1); //probability of choosing each state|simulated states
  
  //store FUTURE simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  mat pchoiceMat    = zeros<mat>(T,2);
  
  double x0,x1;
  //In the paper, states start at i = 1...S (not 0)
  
  vec temp;
  
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
            state[t] = int(std::min(1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0],S + 0.0));
          } else{
            state[t] = int(std::min(state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0],S + 0.0));
          }
          
          //Based on state, make decision for next period
          //Substract 1 for c++ indexing
          pchoice[t]  = pchoicehat[state[t]-1];
          choice[t]   = (pchoice[t] > runif(1)[0])*1;
          
        } //end of T forward simulations
        
        //fill in state, choice, and pchocie matrixes based on action a
        //a = 1 means REPLACE
        stateMat.col(a)   = state.tail(T);
        choiceMat.col(a)  = choice.tail(T);
        pchoiceMat.col(a) = pchoice.tail(T);
        
        stateMat.col(a)   = state.head(T);
        choiceMat.col(a)  = choice.head(T);
        pchoiceMat.col(a) = pchoice.head(T);
        
      } //end of action a loop
      
      //maintain
      x0 = -accu(betaVec % (onesT-pchoiceMat.col(0))%(0 - thetahat[1]*stateMat.col(0) + gamma_cnst - log(onesT-pchoiceMat.col(0))));
      //x0 = -accu(betaVec % choiceMat.col(0)%(0 - thetahat[1]*stateMat.col(0) + gamma_cnst - log(onesT-pchoiceMat.col(0))));
      x0 = -accu(betaVec % (0 - thetahat[1]*stateMat.col(0) + gamma_cnst - log(onesT-pchoiceMat.col(0))));
      
      
      //replace
      x1 = -accu(betaVec % pchoiceMat.col(1)%(-thetahat[0] - thetahat[1]*stateMat.col(1) + gamma_cnst - log(pchoiceMat.col(1))));
      //x1 = -accu(betaVec % choiceMat.col(1)%(-thetahat[0] - thetahat[1]*stateMat.col(1) + gamma_cnst - log(pchoiceMat.col(1))));
      x1 = -accu(betaVec % (-thetahat[0] - thetahat[1]*stateMat.col(1) + gamma_cnst - log(pchoiceMat.col(1))));
      
      
      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/(R + 0.0);
  
  return out;
}



// [[Rcpp::export]]
mat getValDiffhmss(int S, int T, int R, double beta, vec phat, vec pchoicehat) {
  
  vec betaVec = zeros<vec>(T);
  // for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t + 0.0);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t + 1.0);
  
  double gamma_cnst = .5772;
  
  mat out      = zeros<mat>(S,3); //difference in discounted values for both actions
  out = out + 0.0;
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  vec onesT    = ones<vec>(T) + 0.0;
  
  vec state    = zeros<vec>(T + 1);
  vec choice   = zeros<vec>(T + 1);
  vec pchoice  = zeros<vec>(T + 1); //probability of choosing each state|simulated states
  
  //store FUTURE simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  mat pchoiceMat    = zeros<mat>(T,2);
  
  double x0,x1,x2;
  //In the paper, states start at i = 1...S (not 0)
  
  vec temp;
  
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
            state[t] = int(std::min(1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0],S + 0.0));
          } else{
            state[t] = int(std::min(state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0],S + 0.0));
          }
          
          //Based on state, make decision for next period
          //Substract 1 for c++ indexing
          pchoice[t]  = pchoicehat[state[t]-1];
          choice[t]   = (pchoice[t] > runif(1)[0])*1;
          
        } //end of T forward simulations
        
        //fill in state, choice, and pchocie matrixes based on action a
        //a = 1 means REPLACE
        stateMat.col(a)   = state.tail(T);
        choiceMat.col(a)  = choice.tail(T);
        pchoiceMat.col(a) = pchoice.tail(T);
        
      } //end of action a loop
      
      //Equation 5.6
      // x0 = -accu(betaVec % (pchoiceMat.col(1)%(gamma_cnst - log(pchoiceMat.col(1))) +
      //   (onesT-pchoiceMat.col(0))%(gamma_cnst - log(onesT-pchoiceMat.col(0)))));
      
      x0 = -accu(betaVec % (pchoiceMat.col(1)%(-gamma_cnst + log(pchoiceMat.col(1))) +
        (onesT-pchoiceMat.col(0))%(-gamma_cnst + log(onesT-pchoiceMat.col(0)))));
      
      x1 = -accu(betaVec % (onesT + pchoiceMat.col(1) - pchoiceMat.col(0))); //POTENTIAL TYPO
      // x1 = -accu(betaVec % (pchoiceMat.col(1)));
      
      x2 = -accu(betaVec % (pchoiceMat.col(1) % stateMat.col(1) + (onesT-pchoiceMat.col(0)) % stateMat.col(0)));
      
      //place into outEV
      out(s,0) = out(s,0) + x0;
      out(s,1) = out(s,1) + x1;
      out(s,2) = out(s,2) + x2;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/(R + 0.0);
  
  return out;
}


// [[Rcpp::export]]
mat getXHMSSv2(int S, int T, int R, double beta, vec phat, vec pchoicehat,int option = 1) {
  
  vec betaVec = zeros<vec>(T);
  for(int t = 0; t < T; t ++) betaVec[t] = pow(beta,t);
  
  double gamma_cnst = .5772;
  
  mat out      = zeros<mat>(S,3); //difference in discounted values for both actions
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  vec onesT    = ones<vec>(T);
  vec state    = zeros<vec>(T);
  vec choice   = zeros<vec>(T);
  vec pState       = zeros<vec>(T); //probability of choosing each state|simulated states
  vec pStateFix0   = zeros<vec>(T); //probability of choosing each state, except in t = 0 where action is determined
  
  //store simulated states, choices, pchoice|state, and pchoice|state where pchoice[0] is fixed based on a
  mat stateMat      = zeros<mat>(T,2);
  mat choiceMat     = zeros<mat>(T,2);
  mat pStateMat     = zeros<mat>(T,2);
  mat pStateFix0Mat = zeros<mat>(T,2);
  
  vec x0 = zeros<vec>(T); 
  vec x1 = zeros<vec>(T);
  vec x2 = zeros<vec>(T);
  
  //In the paper, states start at i = 1...S (not 0)
  
  //Loop over each state
  for(int s = 0; s < S; s++){
    
    //for R repititions
    for(int r = 0; r < R; r++){
      
      //for each action: 0 = maintain, 1 = replace
      for(int a = 0; a < 2; a++){
        
        //reset "state" and "maintain" vectors
        choice[0] = a;     //when a = 0, maintain, else = 1
        state[0]  = s + 1; //states are from 1...90
        
        pState[0]     = pchoicehat[0];
        pStateFix0[0] = a;
        
        //forward simulate T-1 time periods (first time period pre-determined by action a)
        for(int t = 1; t < T; t++){
          
          //If last period's decision was to replace  (choice = 1), state is reset to 1 plus some transition
          //If last period's decision was to maintain (choice = 0), take the last state plus some transition
          if(choice[t - 1] == 1){
            state[t] = std::min(1 + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0],S + 0.0);
          } else{
            state[t] = std::min(state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,phat)[0],S + 0.0);
          }
          
          //Based on state, make decision for next period
          //Substract 1 for c++ indexing
          pState[t]     = pchoicehat[state[t]-1];
          pStateFix0[t] = pState[t];
          choice[t]     = (pState[t] > runif(1)[0])*1;
          
        } //end of T time periods loop
        
        //fill in state, choice, pState, and pStateFix0 matrixes based on action a
        //a = 1 means REPLACE
        stateMat.col(a)  = state;
        choiceMat.col(a) = choice;
        pStateMat.col(a) = pState;
        pStateFix0Mat.col(a) = pStateFix0;
      
      } //end of action a loop
      
      if(option == 1){ //Standard BBL approach: use realized choice
        x0 = (gamma_cnst - log(pStateMat.col(0))) - (gamma_cnst - log(pStateMat.col(1)));
        x1 = choiceMat.col(0) - choiceMat.col(1);
        x2 = ((onesT - choiceMat.col(0)) % stateMat.col(0)) - ((onesT-choiceMat.col(1)) % stateMat.col(1));
      } else if (option == 2) { //Use probabilities instead of realized choices
        x0 = (gamma_cnst - log(pStateMat.col(0))) - (gamma_cnst - log(pStateMat.col(1)));
        x1 = pStateMat.col(0) - pStateMat.col(1);
        x2 = ((onesT-pStateMat.col(0)) % stateMat.col(0)) - ((onesT-pStateMat.col(1)) % stateMat.col(1));
      } else if (option == 3) { //Use probabilities instead of realized choices BUT fix the first one (except in logs, which is impossible)
        x0 = (gamma_cnst - log(pStateMat.col(0))) - (gamma_cnst - log(pStateMat.col(1)));
        x1 = pStateFix0Mat.col(0) - pStateFix0Mat.col(1);
        x2 = ((onesT-pStateFix0Mat.col(0)) % stateMat.col(0)) - ((onesT-pStateFix0Mat.col(1)) % stateMat.col(1));
      } else if (option == 4) { //HMSS 1994 based on notation
        x0 = pStateMat.col(0) % (gamma_cnst - log(pStateMat.col(0))) + (onesT - pStateMat.col(1)) % (gamma_cnst - log(onesT - pStateMat.col(1)));
        x1 = onesT + pStateMat.col(0) - pStateMat.col(1);
        x2 = pStateMat.col(0) % choiceMat.col(0) + (onesT - pStateMat.col(1)) % choiceMat.col(1);
      } else if (option == 5) { //modified HMSS based on experiments
        x0 = pStateFix0Mat.col(0) % (gamma_cnst - log(pStateMat.col(0))) + (onesT - pStateFix0Mat.col(1)) % (gamma_cnst - log(pStateMat.col(1)));
        x1 = 1 + pStateFix0Mat.col(0) - pStateFix0Mat.col(1);
        x2 = pStateFix0Mat.col(0) % choiceMat.col(0) - (onesT - pStateFix0Mat.col(1)) % choiceMat.col(1);
      }
        
      //place into outEV
      out(s,0) = out(s,0) + sum(betaVec % x0);
      out(s,1) = out(s,1) + sum(betaVec % x1);
      out(s,2) = out(s,2) + sum(betaVec % x2);
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  out = out/R;
  
  return out;
}



// [[Rcpp::export]]
mat getXHMSSt(int S, int T, int R, double beta, vec p, vec pchoicehat) {
  
  vec discFact = zeros<vec>(T);
  for(int t = 0; t < T; t ++) discFact[t] = pow(beta,t);
  
  double gamma_cnst = .5572;
  
  mat outEV    = zeros<mat>(S,3); //discounted expected value for both actions
  mat out      = zeros<mat>(S,2); //difference in discounted values for both actions
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  vec onesT    = ones<vec>(T);
  vec state    = zeros<vec>(T);
  vec maintain = zeros<vec>(T);
  
  vec p_hmss(T);
  vec p1h2;
  vec p1h1;
  vec h1;
  vec h2;
  double x0; 
  double x1;
  double x2;
  
  //In the paper, states start at i = 1...S (not 0)
  
  //Loop over each state
  for(int s = 0; s < S; s++){
    
    //for R repititions
    for(int r = 0; r < R; r++){
      
      //for each action: 0 = maintain, 1 = replace
      for(int a = 0; a < 2; a++){
        
        //reset "state" and "maintain" vectors
        maintain[0] = (1-a); //when a = 0, maintain[0] = 1
        state[0]    = s + 1; //s starts at zero, in the paper states start at 1
        
        //forward simulate T-1 time periods (first time period pre-determined by action a)
        for(int t = 1; t < T; t++){
          
          //If last period maintain = 1, draw a new state based on transition probabilities
          if(maintain[t - 1] == 1){
            state[t] = state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,p)[0];
          } else{
            state[t] = 1 + Rcpp::RcppArmadillo::sample(pStates,1,true,p)[0]; //starting state plus some transition
          }
          
          //Based on state, make decision for next period
          //substract 1 for index to start at 0, not 1
          //phoice is probability of maintaining
          maintain[t] = (pchoicehat[state[t]-1] > runif(1)[0])*1;
          
        } //end of T time periods loop
        
        //fill in p1N
        //subtract 1 so that indexing starts at 0
        for(int t = 0; t < T; t++) p_hmss[t] = pchoicehat[state[t]-1];
        
        //In paper: superscript 1 = replace, 2 = maintain
        if(a == 0){ //if t0 = maintain
          p1h2 = p_hmss;
          h2 = state;
        } else {    //if t0 = replace
          p1h1 = p_hmss;
          h1 = state;
        }
        
      } //end of action a loop
      
      x0 = -sum(discFact % (
        (p1h1 % (gamma_cnst - log(p1h1))) +
        ((1 - p1h2) % (gamma_cnst - log(1 - p1h2)))));
      
      x1 = -sum(discFact % (1 + p1h1 - p1h2));
      
      x2 = -sum(discFact % (p1h1 % h1 + ((1 - p1h2) % h2)));
      
      //place into outEV
      outEV(s,0) = outEV(s,0) + x0;
      outEV(s,1) = outEV(s,1) + x1;
      outEV(s,2) = outEV(s,2) + x2;
      
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  outEV = outEV/R;

  return outEV;
}


// [[Rcpp::export]]
mat getXdiffRust(int S, int T, int R, double beta, vec p, vec pchoicehat) {
  
  vec discFact = zeros<vec>(T);
  for(int t = 0; t < T; t ++) discFact[t] = pow(beta,t);
  
  mat outEV    = zeros<mat>(S,4); //discounted expected value for both actions
  mat out      = zeros<mat>(S,2); //difference in discounted values for both actions
  
  vec pStates  = {0,1,2};      //transition states associated with transition probabilities p
  vec onesT    = ones<vec>(T);
  vec state    = zeros<vec>(T);
  vec maintain = zeros<vec>(T);
  double discMaintain;
  double discReplace;
  
  //Loop over each state
  for(int s = 0; s < S; s++){
    
    //for R repititions
    for(int r = 0; r < R; r++){
      
      //for each action: 0 = maintain, 1 = replace
      for(int a = 0; a < 2; a++){
        
        //reset "state" and "maintain" vectors
        maintain[0] = a;
        state[0]    = s;
        
        //forward simulate T-1 time periods (first time period pre-determined by action a)
        for(int t = 1; t < T; t++){
          
          //If last period maintain = 1, draw a new state based on transition probabilities
          if(maintain[t - 1] == 1){
            state[t] =  state[t-1] + Rcpp::RcppArmadillo::sample(pStates,1,true,p)[0];
          } else{
            state[t] = 0;
          }
          
          //Based on state, make decision for next period
          maintain[t] = (pchoicehat[state[t]] > runif(1)[0])*1;
          
          
        } //end of T time periods loop
        
        //sum the discounted simulated values
        discReplace  = sum(discFact % (onesT - maintain));
        discMaintain = sum(discFact % maintain % state);
        
        //place into outEV depending on time = 0 action a
        if(a == 0){ //maintain
          outEV(s,0) = outEV(s,0) + discReplace;
          outEV(s,1) = outEV(s,1) + discMaintain;
        } else {    //replace
          outEV(s,2) = outEV(s,2) + discReplace;
          outEV(s,3) = outEV(s,3) + discMaintain;
        }
        
      } //end of action a loop
    } //end of R simulations loop
    
  } //end of S states loop
  
  //average over expected values
  outEV = outEV/R;
  
  out.col(0) = outEV.col(2) - outEV.col(0);
  out.col(1) = outEV.col(3) - outEV.col(1); 
  
  return out;
}