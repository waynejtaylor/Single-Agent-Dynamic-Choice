// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <float.h>
using namespace Rcpp;
using namespace arma;

List rwishart_rcpp(int const& nu, mat const& V){

// Wayne Taylor 4/7/2015

// Function to draw from Wishart (nu,V) and IW
 
// W ~ W(nu,V)
// E[W]=nuV

// WI=W^-1
// E[WI]=V^-1/(nu-m-1)
  
  // T has sqrt chisqs on diagonal and normals below diagonal
  int m = V.n_rows;
  mat T = zeros(m,m);
  
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]); //rchisq returns a vectorized object, so using [0] allows for the conversion to double
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
  }}
  
  mat C = trans(T)*chol(V);
  mat CI = solve(trimatu(C),eye(m,m)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  // C is the upper triangular root of Wishart therefore, W=C'C
  // this is the LU decomposition Inv(W) = CICI' Note: this is
  // the UL decomp not LU!
  
  // W is Wishart draw, IW is W^-1
  
  return List::create(
    Named("W") = trans(C) * C,
    Named("IW") = CI * trans(CI),
    Named("C") = C,
    Named("CI") = CI);
}

List rmultireg_rcpp(mat const& Y, mat const& X, mat const& Bbar, mat const& A, int nu, mat const& V) {

// Keunwoo Kim 09/09/2014

// Purpose: draw from posterior for Multivariate Regression Model with natural conjugate prior

// Arguments:
//  Y is n x m matrix
//  X is n x k
//  Bbar is the prior mean of regression coefficients  (k x m)
//  A is prior precision matrix
//  nu, V are parameters for prior on Sigma

// Output: list of B, Sigma draws of matrix of coefficients and Sigma matrix
 
// Model: 
//  Y=XB+U  cov(u_i) = Sigma
//  B is k x m matrix of coefficients

// Prior:  
//  beta|Sigma  ~ N(betabar,Sigma (x) A^-1)
//  betabar=vec(Bbar)
//  beta = vec(B) 
//  Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)

  int n = Y.n_rows;
  int m = Y.n_cols;
  int k = X.n_cols;
  
  //first draw Sigma
  mat RA = chol(A);
  mat W = join_cols(X, RA); //analogous to rbind() in R
  mat Z = join_cols(Y, RA*Bbar);
  // note:  Y,X,A,Bbar must be matrices!
  mat IR = solve(trimatu(chol(trans(W)*W)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  // W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
  mat Btilde = (IR*trans(IR)) * (trans(W)*Z);
  // IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
  mat E = Z-W*Btilde;
  mat S = trans(E)*E;
  // E'E
  
  // compute the inverse of V+S
  mat ucholinv = solve(trimatu(chol(V+S)), eye(m,m));
  mat VSinv = ucholinv*trans(ucholinv);
  
  List rwout = rwishart_rcpp(nu+n, VSinv);
  
  // now draw B given Sigma
  //   note beta ~ N(vec(Btilde),Sigma (x) Covxxa)
  //       Cov=(X'X + A)^-1  = IR t(IR)  
  //       Sigma=CICI'    
  //       therefore, cov(beta)= Omega = CICI' (x) IR IR' = (CI (x) IR) (CI (x) IR)'
  //  so to draw beta we do beta= vec(Btilde) +(CI (x) IR)vec(Z_mk)  
  //    	Z_mk is m x k matrix of N(0,1)
  //	since vec(ABC) = (C' (x) A)vec(B), we have 
  //		B = Btilde + IR Z_mk CI'

  mat CI = rwout["CI"]; //there is no need to use as<mat>(rwout["CI"]) since CI is being initiated as a mat in the same line
  mat draw = mat(rnorm(k*m));
  draw.reshape(k,m);
  mat B = Btilde + IR*draw*trans(CI);
    
  return List::create(
      Named("B") = B, 
      Named("Sigma") = rwout["IW"]);
}

// [[Rcpp::export]]
List bddcMCMCloop(int R,int nsim1,vec const& States,List const&  lgtdata,mat const&  Z,mat const& ST_mat,mat const& R_mat,double beta0,
                  double nu, mat const&  V, mat const&  Deltabar, mat const& ADelta, mat oldVtheta, mat oldVthetai, mat oldDelta, 
                  double stheta, bool fastConv, bool shocks, bool kernelavg, 
                  double scale, mat const& thetas, mat const& Vtheta,
                  int NsAll, int LAll, vec const& hkerns){
    
  int S = States.size();
  int nlgt = lgtdata.size();
  int nvar = V.n_cols;
  int nz = Z.n_cols;
  
  //Storage
  mat emax   = zeros<mat>(S,nlgt);           //Emax values for each state x cross-sectional unit
  mat Deltadraw   = zeros<mat>(R,nvar*nz);
  cube thetadraw  = zeros<cube>(nlgt,nvar,R);
  mat Vthetadraw  = zeros<mat>(R,nvar*nvar);
  
  vec reject = zeros<vec>(R);
  vec llike  = zeros<vec>(R);
  //Initialize variables
  mat oldthetas = zeros<mat>(nlgt,nvar);
  if(fastConv) oldthetas = thetas*scale;  //'thetas' from DGP file
  
  
  int nobs;
  double rej,logl,pold, pnew, logknew, logkold, alpha, u;
  vec y, state, thetao, thetan, val1, val2 ,diff, pval1, prob, alphaminv;
  vec futil_maint,futil_repl,futil_maint_state,futil_repl_state;
  mat sV, root;
  List lgtdatai, multiregout;
  
  //for Kernel estimation
  int Ns, L;
  double kernconstant, theta1_c,theta2_c;
  vec theta1_p,theta2_p,rkern11,rkern22,kernMulti;
  uvec allRowsState(S),indKerni,indTheta,indThetai;
  mat valuei;
  cube value = zeros<cube>(S,nlgt,NsAll);  //Keep the updated emax values from the last Ns MCMC draws (used for kernel weighting)
  
  vec kern = zeros<vec>(NsAll);
  for(int i = 0; i < S; i++) allRowsState[i] = i;
  
  //for shocks
  mat sh     = zeros<mat>(nsim1,2);
  int nhalf = .5*nsim1;
  mat val;
  
  //START OF GIBBS SAMPLER #######################################################################################################
  for(int r = 0; r<R; r++){
    
    if((r+1) % 10 == 0) Rprintf("r = %i \n", r+1);
    
    rej = 0;
    logl = 0;  
    sV = stheta*oldVtheta;
    root = trans(chol(sV));
    if(fastConv){
      oldVtheta  = trans(chol(Vtheta));
      oldVthetai = solve(trimatu(chol(Vtheta)),eye(nvar,nvar));
      root = chol(Vtheta);
    }
    
    //Draw B_i|B-bar, V
    for(int i = 0;i<nlgt;i++){
      lgtdatai = lgtdata[i];
      
      y = as<vec>(lgtdatai["y"]);
      state = as<vec>(lgtdatai["state"]);
      nobs = state.size();
      
      //Draw new thetas
      thetao = vectorise(oldthetas(i,span::all));
      thetan = thetao + root*as<vec>(rnorm(nvar));
      
      //--------------------------------------------------------------
      // Derive expected values
      //--------------------------------------------------------------
      
      futil_maint = vectorise(ST_mat*emax(span::all,i));
      futil_repl  = vectorise(R_mat*emax(span::all,i));
      
      futil_maint_state.set_size(nobs);
      futil_repl_state.set_size(nobs);
      for(int n = 0;n<nobs;n++){
        futil_maint_state[n] = futil_maint[state[n]-1];
        futil_repl_state[n] = futil_repl[state[n]-1];
      }
      //--------------------------------------------------------------
      // Construct the likelihood
      //--------------------------------------------------------------
      
      //USING 'OLD' PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      val1 = thetao[1]*state+beta0*futil_maint_state;
      val2 = thetao[0]+beta0*futil_repl_state;
      diff = -(val1-val2); //Remember, since these are costs, higher values are worse (hence the negative)
      
      pval1 = (exp(diff)/(1+exp(diff)));
      prob = pval1%(1-y) + (1-pval1)%y;
    
      pold = sum(log(prob));
      
      //USING 'NEW' PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      val1 = thetan[1]*state+beta0*futil_maint_state;
      val2 = thetan[0]+beta0*futil_repl_state;
      diff = -(val1-val2);
      
      pval1 = (exp(diff)/(1+exp(diff)));
      prob = pval1%(1-y) + (1-pval1)%y;
    
      pnew = sum(log(prob));
      
      //--------------------------------------------------------------------------
      // Calculate the next iteration parameters of the Metropolis- Hastings algorithm. 
      //--------------------------------------------------------------------------
      
      //heterogeneity
      logknew = as_scalar(-.5*(thetan.t()-Z(i,span::all)*oldDelta) * oldVthetai * (thetan-(Z(i,span::all)*oldDelta).t()));
      logkold = as_scalar(-.5*(thetao.t()-Z(i,span::all)*oldDelta) * oldVthetai * (thetao-(Z(i,span::all)*oldDelta).t()));
      //oldDelta       = nz x nvar
      //Z(i,span::all) = 1 x nz
      //oldVthetai     = nvar x nvar
      //thetan/o       = nvar x 1
      //result is 1xnvar x nvarxnvar x nvarx1 = 1x1
      
      // M-H step
      alpha = exp(pnew + logknew - pold - logkold);
      //if(alpha == "Inf") alpha = -1;
      u = as_scalar(vec(runif(1)));
      if(u < alpha) { 
        oldthetas(i,span::all) = thetan.t();
        logl = logl + pnew;
      } else {
        logl = logl + pold;
        rej = rej+1;
      }
  
    }

    //Draw B-bar and V as a multivariate regression
    multiregout = rmultireg_rcpp(oldthetas,Z,Deltabar,ADelta,nu,V);
    oldDelta    = as<mat>(multiregout["B"]);
    oldVtheta   = as<mat>(multiregout["Sigma"]);
    oldVthetai  = solve(trimatu(chol(oldVtheta)),eye(nvar,nvar));
    
    Deltadraw(r,span::all)  = trans(vectorise(oldDelta));
    Vthetadraw(r,span::all) = trans(vectorise(oldVtheta));
    thetadraw.slice(r)      = oldthetas;
    llike[r]=logl;
    reject[r]=rej/nlgt;
    
    //================================================================
    // This part derives the Emax functions using the simulated data.
    //================================================================
    if(shocks){
      val.set_size(nsim1,2);
      sh.set_size(nsim1,2);
    } else {
      val.set_size(S,2);
    }
                         
    //Loop through each person
    for(int i = 0;i<nlgt;i++){
      
      //----------------------------------------------------------------
      // Derive the transition probabilities.
      //----------------------------------------------------------------
      
      futil_maint = vectorise(ST_mat*emax(span::all,i));
      futil_repl = vectorise(R_mat*emax(span::all,i));
      
      //Are you averaging over shocks?
      if(shocks){
        
        mat rdraws(nhalf, 2, fill::randn);
        sh(span(0,nhalf-1),span::all) = rdraws;
        sh(span(nhalf,nsim1-1),span::all) = -sh(span(0,nhalf-1),span::all);
        
        for(int s = 0; s<S; s++){        
          val(span::all,0) =   thetadraw(i,1,r)*States[s]+beta0*futil_maint[s]+sh(span::all,0); //COST of MAINTAINING
          val(span::all,1) =   thetadraw(i,0,r)+beta0*futil_repl[s]+sh(span::all,1);            //COST of REPLACING
          
          emax(s,i) = as_scalar(sum(min(val,1)))/(double)nsim1;
        }
      } else {
        
        val(span::all,0) =  thetadraw(i,1,r)*States+beta0*futil_maint;  //COST of MAINTAINING
        val(span::all,1) =  thetadraw(i,0,r)+beta0*futil_repl;          //COST of REPLACING
        
        emax(span::all,i) = min(val,1);
      }
    }    
    
    //Note: value(,,0) contains the most recent emax, value(,,1) contains emax from 2 times ago
    if(kernelavg){
      for(int v = 1;v<NsAll;v++) value.slice(v) = value.slice(v-1); //push back the old values
      value.slice(0) = emax;                                        //save the current Emax
    }
    
    if((r>0) & kernelavg){
       
      //---------------------------------------------------------
      // Using the kernel function, update the expected value matrix
      //--------------------------------------------------------
  
      //Note: value(,,1) contains the "last" emax, value(,,2) contains emax from 2 times ago
  
      //We check the last Ns thetadraws and compare the distance with the most recent thetadraw
      Ns = std::min(NsAll,r); //this is the most we can check back over, which might be less than NsAll
      L  = std::min(LAll,r);  //this is the most we can average back over, which might be less than LAll
      
      //These are the indices we will check over      
      indTheta.set_size(Ns);
      for(int n = 0;n<Ns;n++) indTheta[n] = r-n-1;
       
      //----------------------------------
      // Derive the Gaussian Kernel values
      //----------------------------------
      kernconstant = pow(2*M_PI,-L/2.0); //not really needed, since only relative kernel values matter
    
      for(int i = 0;i<nlgt;i++){  
        //current thetas
        theta1_c   =  thetadraw(i,0,r);
        theta2_c   =  thetadraw(i,1,r);
        // MAYBE LATER: theta_c = thetadraw(i,span::all,r);
        
        indThetai = i + (indTheta * (nlgt * nvar));
        
        //past thetas
        theta1_p   =  thetadraw(indThetai);
        theta2_p   =  thetadraw(indThetai + nlgt);
        
        rkern11 = (1/hkerns[0])*exp(-.5*(pow((theta1_c-theta1_p)/hkerns[0],2)));
        rkern22 = (1/hkerns[1])*exp(-.5*(pow((theta2_c-theta2_p)/hkerns[1],2)));
        
        kern = kernconstant*(rkern11%rkern22);
        
        indKerni = sort_index(kern,"descend");  //order N(s) kernel values. Higher values are closer (preferred)
        indKerni = indKerni.head(L);            //keep only the L closet values
        
        //=================================================================
        // Recalculate the Emax function. (see equation 7 on page 1883)
        //=================================================================  
        kernMulti = kern(indKerni);
        //Check to make sure division will work
        if(sum(kernMulti) == 0){
          kernMulti.ones();
        } else {
          kernMulti = kernMulti/sum(kernMulti);
        }
        
        valuei = value(span::all,span(i),span::all);
        valuei = valuei.submat(allRowsState,indKerni);
        
        emax(span::all,i) = valuei * kernMulti; //weight the past values by the kernel density, more weight given to the closer values
      }
    }
    
    //==========================================================
    // Now, generate nsim new data.
    //==========================================================
      
    //(omitted)
  
  }
  
  return(List::create(
    Named("Deltadraw")  = Deltadraw,
    Named("Vthetadraw") = Vthetadraw,
    Named("thetadraw")  = thetadraw,
    Named("llike")      = llike, 
    Named("reject")     = reject,
    Named("emax")       = emax
  ));
  
}