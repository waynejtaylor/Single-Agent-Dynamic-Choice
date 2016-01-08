intcond=function(b,like,X){
  U1=X%*%b
  p=exp(U1)/(1+exp(U1))
  p=cbind(p,1-p)
  
  Like = -sum(log(rowSums(p*like)))
  
  Like
}