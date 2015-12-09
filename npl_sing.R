library(pracma)

npl_sing=function(inda,indx,zmat,pini,bdisc,fmat,names){
  
  "
  Maximum Likelihood Estimates of structural parameters 
  of a discrete choice single-agent dynamic programming 
  model using the NPL algorithm in Aguirregabiria and Mira (Econometrica, 2002)
   
  Original code in GAUSS by Victor Aguirregabiria
  Converted to R by Wayne Taylor

  Version 12/7/2015
  ---------------------------------------------------------------
  
  INPUTS:
    inda    - (nobs x 1) vector with indexes of discrete decision variable (values of 1,...,J)
  
    indx    - (nobs x 1) vector with indexes of the state vector x (values of 1,..,S)
  
    zmat    - (zmat1,zmat2,...,zmatJ) matrix with the values of the variables z(a=j,x)
                note: each zmat has J columns to represent the utility of choice j given action a
  
    pini    - (numx x J) vector with the initial estimates of the choice probabilities Pr(a=j|x)
  
    bdisc   - Discount factor (between 0 and 1)
  
    fmat    - (fmat1,fmat2,...,fmatJ) matrix with the conditional choice transition probs
  
    names   - (npar x 1) vector with names of parameters
  
   OUTPUTS:
    A list of size K where the k'th entry contains:
    
    tetaest - (npar x 1) matrix with estimates of structural parameters of the k'th stage estimate
  
    varest  - (npar x npar) matrix with asymptotic covariance matrices of estimates for the k'th stage
  
    pest    - (numx x J) matrix with the estimated choice probabilities Pr(d=1|x),...,Pr(d=J|x) for the k'th stage
  ---------------------------------------------------------------"
  
  npar = length(names)
  nobs = length(inda)
  nchoice = max(inda)
  if(ncol(zmat)!=(npar*nchoice)){
    print("Error: The number of columns in 'zmat' does not agree",fill=TRUE)
    print("with the number of 'choices * number of parameters'",fill=TRUE)
  }

  myzero = 1e-12
  eulerc = 0.5772
  numx = nrow(pini)
  convcrit = 1000
  convcons = 1e-6
  tetaest0 = matrix(0,npar,1)
  out = NULL

  #---------------------------------------------------------
  #             ESTIMATION OF STRUCTURAL PARAMETERS
  #---------------------------------------------------------
  ks=1
  while(convcrit>=convcons){
  
    cat("-----------------------------------------------------",fill=TRUE)
    cat("POLICY ITERATION ESTIMATOR: STAGE =",ks,fill=TRUE)
    cat("-----------------------------------------------------",fill=TRUE)

    #1. Obtaining matrices "A=(I-beta*Fu)" and "Bz=sumj{Pj*Zj}" and vector Be=sumj{Pj*ej}
    #-----------------------------------------------------------------------------------
  
    i_fu = matrix(0,numx,numx)
    sumpz = matrix(0,numx,npar)
    sumpe = matrix(0,numx,1)
    j=1
    while (j<=nchoice){
      i_fu = i_fu + pini[,j]*fmat[,(numx*(j-1)+1):(numx*j)] #notice the column references
      sumpz = sumpz + pini[,j]*zmat[,(npar*(j-1)+1):(npar*j)]
      sumpe = sumpe + pini[,j]*(eulerc - log(pini[,j]+myzero)) #NOTE  I ADDED +MYZERO so log() works
      j=j+1 ;
    }
      
    i_fu = diag(numx) - bdisc * i_fu
    
    #2. Solving the linear systems "A*Wz = Bz" and "A*We = Be" using CROUT decomposition
    #-----------------------------------------------------------------------------------
  
    i_fu = lu(i_fu)
    wz = solve(i_fu$L,cbind(sumpz,sumpe))
    wz = solve(i_fu$U,wz)
  
    we = wz[,npar+1]
    wz = wz[,1:npar]
  
    #OR:
    # we=solve(i_fu,sumpe)
    # wz=solve(i_fu,sumpz)

    #3. Computing "ztilda(a,x) = z(a,x) + beta * F(a,x)'*Wz" and "etilda(a,x) = beta * F(a,x)'*We"
    #-----------------------------------------------------------------------------------
  
    ztilda = matrix(0,numx,nchoice*npar)
    etilda = matrix(0,numx,nchoice)
    j=1
    while(j<=nchoice){
      ztilda[,(npar*(j-1)+1):(npar*j)] = zmat[,(npar*(j-1)+1):(npar*j)]+bdisc*fmat[,(numx*(j-1)+1):(numx*j)]%*%wz
      etilda[,j] = bdisc * fmat[,(numx*(j-1)+1):(numx*j)]%*%we  
      j=j+1
    }

    #4. Sample observations of "ztilda" and "etilda"
    #-----------------------------------------------------------------------------------
  
    zobs = ztilda[indx,]
    eobs = etilda[indx,]

    #-----------------------------------------------------------------------------------
    #5. Pseudo Maximum Likelihood Estimation
    
    clogitout=clogit(inda,zobs,eobs,names)
    tetaest1=clogitout$b0
    varest=clogitout$Avarb
    
    #6. Re-Computing probabilities
    #-----------------------------------------------------------------------------------
    
    pini = matrix(0,numx,nchoice)
    j=1
    while(j<=nchoice){
      pini[,j] = ztilda[,(npar*(j-1)+1):(npar*j)]%*%tetaest1 + etilda[,j]
      j=j+1
    }
    pini = pini - apply(pini,1,max)
    pini = exp(pini)
    pini = pini/rowSums(pini)
      
    #7. Convergence Criterion
    #-----------------------------------------------------------------------------------
    convcrit = max(abs(tetaest1-tetaest0))
    tetaest0 = tetaest1
    cat("NPL Criterion =",convcrit,fill=TRUE)

    #8. Save output from current k'th stage
    #------------------------------------------------------------------------------------
    out[[ks]]=list(tetaest=tetaest1,varest=varest,pini=pini)
    
    ks=ks+1
  }
  
  out
}