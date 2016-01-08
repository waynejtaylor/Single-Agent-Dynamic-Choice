// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec fvdataRcpp(vec const& b1,mat const& RX1,int tbin,int xbin,vec const& Zstate,
               mat const& Xstate, mat const& xtran, int N,int T,vec const& State, bool hetero = false){

  //calculates fv terms given the data 
  cube FV1 = zeros<cube>(tbin,2,T+1);
  mat FVT1nohet = zeros<mat>(N,T);
  cube FVT1het = zeros<cube>(N,T,2);
  vec FVT1out;
  
  int adj, adj0, z2, z3;
  vec tvec(3);
  mat U1;

  for(int t=2;t<=T;t++){    
    for(int s=0;s<=1;s++){
      tvec[0] = 1;
      tvec[1] = t/10.;
      tvec[2] = tvec[1]*tvec[1];
      U1=kron(trans(tvec),join_rows(RX1,s*RX1))*b1;
      FV1(span(),span(s),span(t-1)) = -log(exp(U1)/(1+exp(U1))); //replaced call to 'plogit'
    }
  }
  
  for(int n = 0;n<N;n++){
    adj0=(Zstate[n]-1)*xbin;
    z2=(Zstate[n]-1)*xbin;
    z3=z2+xbin-1;
    
    for(int t=0;t<T;t++){
      adj=z2+Xstate(n,t)-1;
      if(!hetero){
          FVT1nohet(n,t)=as_scalar((xtran(span(adj),span())-xtran(span(adj0),span()))*vectorise(FV1(span(z2,z3),span(State[n]),span(t+1))));  
      } else {
        for(int s=0;s<2;s++){
          FVT1het(n,t,s)=as_scalar((xtran(span(adj),span())-xtran(span(adj0),span()))*vectorise(FV1(span(z2,z3),span(s),span(t+1)))); 
        }
      }
    }
  }
  
  if(!hetero){
    FVT1out = vectorise(FVT1nohet);  
  } else {
    FVT1out = vectorise(FVT1het);
  }
  
  return(FVT1out);
}