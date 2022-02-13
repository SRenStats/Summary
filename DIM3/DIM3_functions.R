####F1 crtList####
#create a list called dim3 to replace structure array using another function crtList
crtList<-function(datas,numClass){
  dataNum <- dim(datas)[1]
  tTime <- dim(datas)[3]
  lam1<-sum(datas)/(dataNum^2*tTime)
  lam2<-1-lam1
  seLabel<-array(dim=c(dataNum,dataNum,tTime))#label as a sender
  reLabel<-seLabel
  for (t in 1:tTime)
  {
    for (i in 1:dataNum)
    {
      se_value<-rmultinom(dataNum,1, rep(1/numClass,numClass))#a numClass*dataNum dim matrix showing label for all nodes as sender
      re_value<-rmultinom(dataNum,1, rep(1/numClass,numClass))
      for (j in 1:dataNum)
      {
        seLabel[i,j,t]<-which(se_value[,j]==1)#find the label
        reLabel[i,j,t]<-which(re_value[,j]==1)
      }
    }
  }
  m_val <- rep(1,numClass)
  list(datas=datas,lam1=lam1,lam2=lam2,numClass=numClass,dataNum=dataNum,tTime=tTime,reLabel=reLabel,seLabel=seLabel,m_val=m_val)
}

####F2 rdirichlet####
rdirichlet<-function(n,para){
  p <- length(para)
  r <- rgamma(n*p, shape=para, scale = 1)
  r <- r/sum(r)
  c(r)}

####F3 matrix_plus####
mt_pl <- function(A,a){
  A<-as.matrix(A)
  B<-matrix(a,nrow=nrow(A),ncol = ncol(A))
  D<- A+B
  return(D)
}

####F4 count####
count<-function(x,label){
  ct<-rep(0,length(label))
  for (i in 1:length(label)){
    ct[i]<-sum(x==label[i],na.rm=TRUE)
  }
  c(ct)
}

####F5 stirling####
stirling<-function(nn){
  allss<-list()
  allss[[1]]<-1
  if (nn>1){
  for (mm in 2:nn){
    allss[[mm]]<-c(allss[[mm-1]]*(mm-1),0)+c(0,allss[[mm-1]])
    allss[[mm]]<-allss[[mm]]/max(allss[[mm]])
  }}
  return(allss[[nn]])
}

####F6 randconparam####
randconparam<-function(alpha,x,num,numclass,aa,bb,numiter){
  # Modification of Escobar and West.  Works for multiple groups of data.
  # numdata, numclass are row vectors, one element per group.
totalclass <- sum(numclass)

for (ii in 1:numiter){#repeat for accuracy?
  xx <- rbeta(num,(alpha+1),x)
  z1 <- runif(num)
  z2 <- as.vector(x/(alpha+x))
  gammaa <- aa + totalclass - sum(z1<z2)
  gammab <- bb - sum(log(xx))
  alpha <- rgamma(1,gammaa)/gammab
}

return(alpha)
}

####F7 Hyperpara####
Hyperpara <- function(dim3){
  
  alpha_kappa <- randconparam((dim3$alpha+dim3$kappa[2]*2*dim3$dataNum), 2*dim3$dataNum,dim3$dataNum*dim3$tTime,sum(dim3$nohat_m), 1, 1 , 10)
  ratio <- rbeta(1,1+sum(dim3$m_val),1+sum(dim3$nohat_m)-sum(dim3$m_val))
  alphas <- alpha_kappa*ratio
  kappas <- (alpha_kappa - dim3$alpha)/(2*dim3$dataNum)
  gammas <- randconparam(dim3$gamma, sum(dim3$m_val),1,dim3$numClass, 1, 1, 10)
  dim3$alpha <- alphas
  dim3$gamma <- gammas
  dim3$kappa[2] <- kappas
  
  return(dim3)
}

