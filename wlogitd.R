wlogitd=function(b,Y,X,P,LikOnly=TRUE){
  U1=X%*%b
  Like=t(P)%*%(log(1+exp(U1))-Y*U1)
  
  if(LikOnly){
    Like
  } else {
    dg=t(P*((1-(1/(1+exp(U1))))-Y))%*%X
    list(Like=Like,dg=dg)  
  }
  
}
