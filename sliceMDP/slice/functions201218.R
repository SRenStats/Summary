#1-combine array into a matirx
combine<-function(array){
  data<-array
  n<-dim(data)[1]
  k<-dim(data)[2]-1
  T<-dim(data)[3]
  x<-matrix(0,nrow=n*T,ncol=k*T)
  y<-c()
  for (t in 1:T){
    x[(t*n-n+1):(t*n),(t*k-k+1):(t*k)]<-data[,-1,t]
    y[(t*n-n+1):(t*n)]<-data[,1,t]
  }
  return(cbind(y,x))
}
#2-thresholding
th<-function(z,lam){
  f<-z
  for(i in 1:length(z)){
    if(z[i]>lam){
      f[i]=z[i]-lam}
    else if(z[i]<(-lam)){
      f[i]=z[i]+lam}
    else{
      f[i]=0}	
  }				
  return(f)}
#3-main function
admm<-function(x,T,k,ep1,max,lam,rho){
  x_train<-x[,-1]
  y_train<-x[,1]
  p<-(T-1)*k
  h0<-matrix(0,p,k)
  F<-cbind(-diag(p),h0)+cbind(h0,diag(p))
  beta_k1<-rep(0,ncol(x_train))
  q_k1<-rep(0,ncol(x_train)-k)
  u_k1<-rep(0,ncol(x_train)-k)
  m<-0
  repeat{
    beta_k<-beta_k1
    q_k<-q_k1
    u_k<-u_k1
    ep0<-0.001
    #########update loss function
    py<-exp(x_train %*% beta_k)
    py_k<-py/(1+py)
    w_k<-rep(1,nrow(x_train))
    for (i in 1:nrow(x_train)) {
      if (py_k[i]>=1-ep0)
      { py_k[i]<-1
      w_k[i]<-ep0 }
      else if (py_k[i]<=ep0)
      { py_k[i]<-0
      w_k[i]<-ep0 }
      else {
        w_k[i]<-py_k[i]*(1-py_k[i]) }
    }
    z_k<-(x_train %*% beta_k)+((y_train-py_k)/w_k)
    b_k<-z_k*sqrt(w_k)
    wm_k<-matrix(sqrt(w_k),nrow(x_train),ncol(x_train))
    A_k<-wm_k*x_train
    ####update coefficients
    beta_k1<-(solve(t(A_k)%*%A_k+rho*(t(F)%*%F)))%*%((rho*t(F)%*%(q_k-u_k))+(t(A_k)%*%b_k))
    q_k1<-th(((F%*%beta_k1)+u_k),(lam/rho))
    u_k1<-u_k+(F%*%beta_k1)-q_k1
    m<-m+1
    #if ((crossprod(beta_k1-beta_k)<ep1)|(m>max))
    if (crossprod(beta_k1-beta_k)<ep1)
      break
  }
  beta_<-c(lam,rho,m,beta_k1)
  return(beta_)  
}
#4-??????lam,??????rho,?????????????????????
admm.error<-function(x,T,k,ep1,max,lams,rho,beta)
{
  n<-length(lams)
  error<-NULL
  for (s in 1:n) 
  {
    beta_h<-admm(x,T,k,ep1,max,lams[s],rho)
    error1<-mean(abs(beta_h[-1:-3]-beta))##????????????????????????
    error2<-crossprod(beta_h[-1:-3]-beta)/length(beta)##????????????????????????
    error3<-max(abs(beta_h[-1:-3]-beta))##????????????????????????
    error4<-c(error1,error2,error3)
    error<-rbind(error,error4)
  }
  return(error) 
}

###?????????lam,rho?????????,??????rho?????????
admm.error<-function(x,T,k,ep1,max,lams,rhos,beta)
{
  m<-length(rhos)
  n<-length(lams)
  error<-NULL
  for (r in 1:m)
  {
    for (s in 1:n) 
    {
      beta_h<-admm(x,T,k,ep1,max,lams[s],rhos[r])
      error1<-mean(abs(beta_h[-1:-3]-beta))##????????????????????????
      error2<-crossprod(beta_h[-1:-3]-beta)/length(beta)##????????????????????????
      error3<-max(abs(beta_h[-1:-3]-beta))##????????????????????????
      error4<-c(beta_h[1:3],error1,error2,error3)
      error<-rbind(error,error4)
    }
  }
  return(error) 
}



admm.msellikebic<-function(x,T,k,ep1,max,lams,rho,beta)
{
  n<-length(lams)
  beta_h<-matrix(nrow=k*T+3,ncol=n)
  mse<-c()
  llike<-c()
  bic<-c()
  for (l in 1:n) 
  {
    beta_h[,l]<-admm(x,T,k,ep1,max,lams[l],rho)
    mse[l]<-crossprod(beta_h[-1:-3,l]-beta)##????????????????????????
    fx_h<-x[,-1]%*% beta_h[-1:-3,l]
    px_h<-exp(fx_h)
    llike[l]<-sum((x[,1]*fx_h)-log(1+px_h))##????????????????????????
    b1<-beta_h[-1:-3,l]
    b2<-b1[-(T*k-k+1):-(T*k)]-b1[-1:-k]
    bic[l]<--2*llike[l]+log(nrow(x))*(T*k-sum(ifelse(b2==0,1,0)))
  }
  return(rbind(mse,llike,bic,beta_h)) 
}