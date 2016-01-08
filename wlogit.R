wlogit=function(b,Y,X,P){
  
  U1=X%*%b
  Like=t(P)%*%(log(1+exp(U1))-Y*U1)
  
  Like
}

