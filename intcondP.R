intcondP=function(b,like,X){
  U1=X%*%b
  p=exp(U1)/(1+exp(U1))
  p=cbind(p,1-p)
  
  PType = (like*p)/(rowSums(p*like)%*%matrix(1,1,2))
  
  PType
}