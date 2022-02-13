gibbs_dev <- function(dim3){
  deviance <- 0
  cu_jointp <- 0
  se_Labels <- dim3$seLabel
  re_Labels <- dim3$reLabel
  numClass <- dim3$numClass
  indexLabel <- dim3$indexLabel
  alpha <- dim3$alpha
  betas <- dim3$betas
  li_jps <- 0
  
  
  for (k_dev in 1:(numClass-1)){
    if (k_dev > 1){
      k_beta <- betas[k_dev]/(1-sum(betas[1:(k_dev-1)]))
    }else{
      k_beta <- betas[1]}
    
    cu_jointp <- cu_jointp+dbeta(k_beta, 1, dim3$gamma,log=TRUE)
  }
  
  pre_Nikt <- matrix(0,dim3$dataNum, numClass)
  
  for (t_dev in 1:dim3$tTime){
    tau_kl <- matrix(0,numClass, numClass)
    tau1_kl <- matrix(0,numClass, numClass)
    
    for (k in 1:numClass){
      for (l in 1:numClass)
      { x_loc<-which((se_Labels[,,t_dev]==indexLabel[k])&(re_Labels[,,t_dev]==indexLabel[l]),arr.ind = TRUE)
      tau1_kl[k,l] <- sum(diag(dim3$datas[x_loc[,1],x_loc[,2],t_dev]))#number of actual connections from group k to l at time t_sli
      tau_kl[k,l]<- nrow(x_loc)#number of possible connections from group k to l
      }}
    
    Nikt <- matrix(0,dim3$dataNum, numClass)
    for (i_dev in 1:dim3$dataNum){
      Nikt[i_dev,]<-count(c(se_Labels[i_dev,,t_dev],re_Labels[,i_dev,t_dev]),indexLabel)#number of nodes in each cluster
    }
    
      for (i_dev in 1:dim3$dataNum){
      for (j_dev in 1:dim3$dataNum){
        seL <- which(indexLabel==se_Labels[i_dev, j_dev, t_dev])
       #why not seL <- se_Labels[i_dev, j_dev, t_dev]
        reL <- which(indexLabel==re_Labels[i_dev, j_dev, t_dev])
        
        tau_kl[seL, reL]<-tau_kl[seL, reL]-1
        tau1_kl[seL, reL]<- tau1_kl[seL, reL]-dim3$datas[i_dev,j_dev,t_dev]
        tau0_kl <- tau_kl-tau1_kl
        
        
        if (dim3$datas[i_dev,j_dev,t_dev]==1){
          like_wei <- mt_pl(tau1_kl,dim3$lam1)/mt_pl(tau_kl,1)
        }else{
          like_wei <- mt_pl(tau0_kl,dim3$lam2)/mt_pl(tau_kl,1)}
        wei_ij <- diag(Nikt[i_dev, ])%*%like_wei%*%diag(Nikt[j_dev, ])/as.vector(4*(dim3$dataNum)^2*dim3$tTime)
        # why 4 rather than  dim3$numClass
        #wei_ij <- diag(Nikt[i_dev, ])%*%like_wei%*%diag(Nikt[j_dev, ])/as.vector(dim3$numClass*(dim3$dataNum)^2*dim3$tTime)
        deviance <- deviance+log(sum(wei_ij))
        
        
        # log joint probability's calculation
        Nikt[i_dev, seL] <- Nikt[i_dev, seL]-1
        Nikt[j_dev, reL] <- Nikt[j_dev, reL]-1
        
        if (dim3$datas[i_dev,j_dev,t_dev]==1)
          like_wei <- (tau1_kl[seL, reL]+dim3$lam1)/(tau_kl[seL, reL]+1)
        else
          like_wei <- (tau0_kl[seL, reL]+dim3$lam2)/(tau_kl[seL, reL]+1)
        
        if (t_dev > 1){
          ps <- (alpha*betas[seL]+dim3$kappa[2]*pre_Nikt[i_dev, seL]+Nikt[i_dev, seL])/(alpha+(dim3$kappa[2]+1)*2*dim3$dataNum-1)
          pr <- (alpha*betas[reL]+dim3$kappa[2]*pre_Nikt[j_dev, reL]+Nikt[j_dev, reL])/(alpha+(dim3$kappa[2]+1)*2*dim3$dataNum-1)
        }else{
          ps <- (alpha*betas[seL]+Nikt[i_dev, seL])/(alpha+2*dim3$dataNum-1)
          pr <- (alpha*betas[reL]+Nikt[j_dev, reL])/(alpha+2*dim3$dataNum-1)}
        
        cu_jointp <- cu_jointp+log(like_wei)+log(ps)+log(pr)
        
        
        li_jps <- li_jps + log(like_wei)
        
        
        Nikt[i_dev, seL] <- Nikt[i_dev, seL]+1
        Nikt[j_dev, reL] <- Nikt[j_dev, reL]+1
        
        tau_kl[seL, reL]<-tau_kl[seL, reL]+1
        tau1_kl[seL, reL]<- tau1_kl[seL, reL]+dim3$datas[i_dev,j_dev,t_dev]
      }
    }
    pre_Nikt <- Nikt
  }
  deviance <- -2*deviance
  
  
  dim3$deviance <- deviance
  dim3$cu_jointp <- cu_jointp 
  dim3$li_jps <- li_jps
  
  return(dim3)
}


