##########################################################################
#from "intconddata"
#This program recreates Table 1 columns 2,3,5, and 6 from
  #"Conditional Choice Probability Estimation of Dynamic Discrete Choice Models 
  #with Unobserved Heterogeneity" by Arcidiacono and Miller (2011)

#The original code is by the above authors

#Coversion to R was done by Wayne Taylor

#Version 1/5/2016

#Note: The intent of this code was to remain consistent with the original code in terms of structure and naming conventions
  #Therefore, the code has not been optimtized for speed

##########################################################################

library(Rcpp)

source('xgrid.R')
source('wlogitd.R')
source('wlogit.R')
source('likebusML4.R')
sourceCpp('genbus4.cpp')
sourceCpp('fvdataBOTH.cpp')
source('intcond.R')
source('intcondP.R')

set.seed(1)

#true parameter values
alpha=c(2,-.15,1,.9,.4) #Intercept, mileage, heterogeneity, discount factor, Pi
tol=.0000001

MCiter=3       #Monte Carlo iterations
FIML = FALSE   #estimate FIML too? (it takes much longer than CCP)
hetero = FALSE #Is heterogeneity observed? FALSE = cols 1 and 2 TRUE = cols 5 and 6
T=200          #Time periods
if(hetero) T=T/10
N=1000        #Observations per time period

Bccp=NULL #CCP parameter storage
Tccp=NULL #CCP timing
Bfl=NULL  #FIML parameter storage
Tfl=NULL  #FIML timing

if(hetero){
  Lccp=NULL
  LFl=NULL
  Iccp=NULL
  Binit=NULL
}

#Create transition matrices
zval=seq(.25,1.25,.01)
zbin=length(zval)
xval=seq(0,25,.125)
xbin=length(xval)
xtran=matrix(0,zbin*xbin,xbin)
xtranc=array(0,c(xbin,xbin,zbin))
for(z in 1:zbin){
  temp=xgrid(zval[z],xval)
  xtran[(1+(z-1)*xbin):(z*xbin),] = temp$xtran
  xtranc[,,z] = temp$xtranc
}

xtrancRcpp = matrix(xtranc,xbin,xbin*zbin)
tbin=xbin*zbin

#z and x values for each state
zvalr=kronecker(zval,rep(1,xbin))
xvalr=kronecker(rep(1,zbin),xval)/10

#data for reduced form logits
#covers the state space
RX1=cbind(rep(1,zbin*xbin),xvalr,zvalr,xvalr*zvalr,xvalr*xvalr,zvalr*zvalr)

#monte carlos

#starting values for FIML and CCP
alphaf= c(alpha[1:3],log(alpha[4])-log(1-alpha[4]))
alphac= alpha[1:4]

MC=1
while(MC <= MCiter){

  #generating the data 
  genbusout=genbusRcpp(alpha,N,T,xtran,xtrancRcpp,xbin,zbin,xval,zval)
  
  Y=genbusout$Y
  X=genbusout$X
  Z=genbusout$Z
  Xstate=genbusout$Xstate
  Zstate=genbusout$Zstate
  State=genbusout$State
  FVT=genbusout$FVT
  
  y2=as.vector(Y)
  x2=as.vector(X[,1:T])/10
  z2=kronecker(rep(1,T),Z)
  s2=kronecker(rep(1,T),State)
  t2=kronecker(1:T,rep(1,N))/10
  
  if(hetero){
    y2=c(y2,y2)
    x2=c(x2,x2)
    z2=c(z2,z2)
    s2=c(rep(0,N*T),rep(1,N*T)) #restated
    t2=c(t2,t2)
    stemp=c(rep(0,N),rep(1,N))
  }
  
  #estimating FIML----
  if(FIML){
    
    tic = proc.time()[3] #start the timer
    
    if(!hetero){
      bfl=optim(alphaf,likebusML4,Y=Y,State=s2,N=N,T=T,X=X,Zstate=Zstate,Xstate=Xstate,xtran=xtran,tbin=tbin,zbin=zbin,xbin=xbin,xval=xval,Z=Z)
    } else {
      bfl=optim(alphaf,likebusML4,Y=y2,State=stemp,N=N,T=T,X=x2,Zstate=c(Zstate,Zstate),Xstate=c(Xstate,Xstate),xtran=xtran,tbin=tbin,zbin=zbin,xbin=xbin,xval=xval,Z=Z)
    }
    
    toc = proc.time()[3]-tic
    
    Tfl=c(Tfl,toc)
    Bfl=rbind(Bfl,bfl)
  }
  
  #estimating with data ccps----
  
  tic = proc.time()[3] #start the timer
  
  #setting up data for reduced form logit
  xx=cbind(rep(1,N*T),x2,z2,x2*z2,x2*x2,z2*z2,s2,s2*x2,s2*z2,s2*x2*z2,s2*x2*x2,s2*z2*z2)
  xx=cbind(xx,matrix(rep(t2,12),ncol=12)*xx,matrix(rep(t2*t2,12),ncol=12)*xx)
  
  #estimating reduced form logit
  if(!hetero){
    #   b1=rep(0,ncol(xx))
    #   b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=rep(1,N*T),method="BFGS")$par #default fnscale = 1 = min which is what we want
    b1=glm(((y2==0)*1)~xx-1,family='binomial')$coef #MUCH faster, same result
  } else {
    PType=.5*rep(1,2*N*T)
    oPType=rep(0,2*N*T)
    Pi2=c(.5,.5)  
    
    b1 = rep(0,ncol(xx))
    b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=PType,method="BFGS")$par
    #For a binomial GLM prior weights are used to give the number of trials when the response is the proportion of successes
    #So we cannot send them into the "weights" argument
  }
  
  #calculating fv terms
  if(!hetero){
    fvt1 = fvdataRcpp(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,State)  
  } else {
    fvt1 = fvdataRcpp(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,rep(1,N),hetero)  
  }
  
  #estimating the structural parameters
  xccp = cbind(rep(1,N*T),x2*10,s2)
  
  if(!hetero){
    #bccp = alphac  
    #bccp = optim(bccp,wlogit,Y=y2,X=cbind(xccp,fvt1),P=rep(1,N*T),method="BFGS")$par
    bccp = glm(y2~cbind(xccp,fvt1)-1,family='binomial')$coef #much faster than 'optim'
  } else {
    
    #starting the EM algorithm
    
    j=0
    
    bccp = alphac
    
    intcondX=cbind(rep(1,N), X[1:N,1],Z[1:N,1])
    binit=rep(0,3)
    
    cond=0 
    lp=NULL
    while(cond==0){
      
      #updating PType
      ##first getting the type-specific likelihoods
      
      oPType=PType
      
      #replaces call to "likeCPP"
      U1 = cbind(xccp,fvt1)%*%bccp
      Like = (y2*exp(U1)+(1-y2))/(1+exp(U1))
      
      Like2=array(Like,c(N,T,2))
      base=apply(Like2,c(1,3),prod)
      
      #now getting the initial condition parameters 
      intcond_optim=optim(binit,intcond,like=base,X=intcondX,method="BFGS")
      binit = intcond_optim$par
      lp=c(lp,intcond_optim$value)
      
      #and the PType's
      PType=intcondP(binit,base,intcondX)
      PType=kronecker(rep(1,T),PType)
      PType=as.vector(PType)
      
      #estimating reduced form logit
      b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=PType,method="BFGS")$par   
      
      #calculating fv terms
      fvt1=fvdataRcpp(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,rep(1,N),hetero)
      
      bccp = optim(bccp,wlogit,Y=y2,X=cbind(xccp,fvt1),P=PType,method="BFGS")$par 
      
      #CHECKING CONVERGENCE
      if(j>26){
        junk=abs((lp[j]-lp[j-25])/lp[j])<tol
        junk2=abs((lp[j-1]-lp[j-26])/lp[j-1])<tol
        
        cond=junk2*junk
        
        if(j>1000) cond=1
      }
      
      j=j+1
      cat("j: ",j,fill=TRUE)
    }
  }
    
  toc = proc.time()[3]-tic

  Tccp=c(Tccp,toc)
  Bccp=rbind(Bccp,bccp)
  
  if(hetero){
    Iccp=c(Iccp,j)
    Binit=rbind(Binit,binit)  
  }
  
  cat("MC ",MC, " completed",fill=TRUE)
  MC = MC+1
}

Bccp
Tccp