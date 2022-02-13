label_slice<-function(dim3)#slice sampling for label
{
  se_Labels <- dim3$seLabel #sender label
  re_Labels <- dim3$reLabel #reciever label
  dataNum <- dim3$dataNum #data number
  numClass <- dim3$numClass #class number
  indexLabel <- dim3$indexLabel
  alpha <- dim3$alpha
  betas <- dim3$betas  # betas's value need to be re-defined
  
  for (t_sli in 1:dim3$tTime){
  #t_sli <- 1
    tau_kl <- matrix(0,ncol=numClass, nrow=numClass)  #counts for inter-cluster occurancy
    tau1_kl <- matrix(0,ncol=numClass, nrow=numClass) #counts for inter-cluster data occurancy

    for (k in 1:numClass){
      for (l in 1:numClass){
        x_loc<-which((se_Labels[,,t_sli]==indexLabel[k])&(re_Labels[,,t_sli]==indexLabel[l]),arr.ind = TRUE)
        tau1_kl[k,l] <- sum(diag(dim3$datas[x_loc[,1],x_loc[,2],t_sli]))#number of actual connections from group k to l at time t_sli
        tau_kl[k,l]<- nrow(x_loc)#number of possible connections from group k to l
      }
    }
    
    #count each class's data size in each time
    Nikt <- matrix(0,nrow = dataNum, ncol = numClass)
    la_Nikt <- matrix(0,nrow = dataNum, ncol = numClass)
    qpi_value <- matrix(0,nrow = dataNum, ncol = numClass+1)#\omega's value in Maria's paper
    
    # initialize the \pi's value
    for (i_sli in 1:dataNum)
    {
      Nikt[i_sli,]<-count(c(se_Labels[i_sli,,t_sli],re_Labels[,i_sli,t_sli]),indexLabel)#number of nodes in each cluster
      if (t_sli<dim3$tTime)
        la_Nikt[i_sli,]<-count(c(se_Labels[i_sli,,t_sli+1],re_Labels[,i_sli,t_sli+1]),indexLabel)#number of nodes in each cluster at next time point
    }
    
    
    ####step 1 update v,w####

    if (t_sli > 1){
      a_ik <- matrix(alpha*betas[1:numClass], nrow=dataNum,ncol=numClass,byrow = TRUE)+Nikt+dim3$kappa[2]*pre_Nikt # beta distribution's parameters
      b_ik <- matrix(alpha*(1-cumsum(betas[1:numClass])),nrow=dataNum,ncol=numClass,byrow = TRUE)+matrix(2*dataNum,nrow=dataNum,ncol=numClass)-t(apply(Nikt,1,cumsum))+dim3$kappa[2]*(matrix(2*dataNum,nrow=dataNum,ncol=numClass)-t(apply(pre_Nikt,1,cumsum)))
    }else{
      a_ik <- matrix(alpha*betas[1:numClass], nrow=dataNum,ncol=numClass,byrow = TRUE)+Nikt # beta distribution's parameters
      b_ik <- matrix(alpha*(1-cumsum(betas[1:numClass])),nrow=dataNum,ncol=numClass,byrow = TRUE)+matrix(2*dataNum,nrow=dataNum,ncol=numClass)-t(apply(Nikt,1,cumsum))
    }
    pi_value <- matrix(nrow=dataNum,ncol=numClass)
    for (ip in 1:dataNum){
      for (jp in 1:numClass){
        pi_value[ip,jp] <- rbeta(1,a_ik[ip,jp],b_ik[ip,jp]) ###step 1 pi~beta(a,b)
      }
    }
    qpi_value[,1:numClass] <- pi_value*t(apply((matrix(1,nrow=dataNum,ncol=numClass)-cbind(0,pi_value[,-numClass])),1,cumprod))
    qpi_value[,numClass+1] <- 1-rowSums(qpi_value)#extra column for weight
    
    
    ####step 2 update u####
    f2 <- matrix(0,nrow=dataNum, ncol=2*dataNum)
    us_vals <- matrix(0,nrow=dataNum,ncol=2*dataNum)
    for (i_sli in 1:dataNum){
      i_vals <- qpi_value[i_sli,1:numClass]
      f2[i_sli,] <- c(se_Labels[i_sli,,t_sli],re_Labels[,i_sli,t_sli])
      us_vals[i_sli,] <- i_vals[f2[i_sli,]]
    }
    uij<-matrix(runif(2*dataNum^2),nrow = dataNum,ncol=2*dataNum)*us_vals
    rm(us_vals,i_vals)
    
    ####step 3 update no. clusters####
    for (i_sli in 1:dataNum){
      while (qpi_value[i_sli,ncol(qpi_value)]>min(uij[i_sli,],na.rm=TRUE)){#condition of adding new class
        betas[length(betas)]<-betas[length(betas)]*rbeta(1,1,dim3$gamma)#change the old end of betas
        ai <- alpha*betas[length(betas)]#new last 2 of betas
        betas<-c(betas,1-sum(betas))#add new end to betas
        bi <- alpha*betas[length(betas)]#new end of betas
        qpi_value[,ncol(qpi_value)]<-qpi_value[,ncol(qpi_value)]*rbeta(dataNum,ai,bi)#change the end col of qpi
        qpi_value<-cbind(qpi_value,1-rowSums(qpi_value))#add new col to qpi
        Nikt <- cbind(Nikt,0)
        la_Nikt <- cbind(la_Nikt,0)
        numClass <- numClass + 1
        indexLabel <- c(1:numClass)
        tau_kl <- rbind(cbind(tau_kl,0),0)
        tau1_kl <- rbind(cbind(tau1_kl,0),0)
      }
    }
    
    for (i_sli in 1:dataNum){
      for (j_sli in 1:dataNum){
        #i_sli<-1
        #j_sli<-1
        n_tau <- tau_kl
        n_tau1 <- tau1_kl
        # the related label
        seLa <- f2[i_sli, j_sli]#sij send
        reLa <- f2[j_sli, (dataNum+i_sli)]#rij recieve
        
        Nikt[i_sli, seLa] <- Nikt[i_sli, seLa] -1#delete before new allocation
        Nikt[j_sli, reLa] <- Nikt[j_sli, reLa] -1#no. of times that node i in this group-1
        
        ui <- uij[i_sli, j_sli] #send
        uj <- uij[j_sli, (dataNum+i_sli)]#recieve
        
        
        i_label <- which(qpi_value[i_sli,]>ui)
        j_label <- which(qpi_value[j_sli,]>uj)#select the group with bigger weight
        n_tau[seLa, reLa] <- n_tau[seLa, reLa]-1
        n_tau1[seLa, reLa] <- n_tau1[seLa, reLa]-dim3$datas[i_sli,j_sli,t_sli]
        n_tau0 <- n_tau - n_tau1
        #update theta1 is_wei
        if (t_sli<dim3$tTime){
          cu_Ni <- Nikt[i_sli,]
          la_Ni <- la_Nikt[i_sli,]
          la_wei1 <- gamma((alpha*betas[i_label])+la_Ni[i_label]+dim3$kappa[2]*cu_Ni[i_label]+dim3$kappa[2])/gamma(alpha*betas[i_label]+la_Ni[i_label]+(dim3$kappa[2]*cu_Ni[i_label]))
          la_wei2 <- gamma(alpha*betas[i_label]+dim3$kappa[2]*cu_Ni[i_label])/gamma(alpha*betas[i_label]+dim3$kappa[2]*cu_Ni[i_label]+dim3$kappa[2])
          is_wei <- la_wei1*la_wei2
          #update theta2 ir_wei
          cu_Ni <- Nikt[j_sli,]
          la_Ni <- la_Nikt[j_sli,]
          la_wei1 <- gamma(alpha*betas[j_label]+la_Ni[j_label]+dim3$kappa[2]*cu_Ni[j_label]+dim3$kappa[2])/gamma(alpha*betas[j_label]+la_Ni[j_label]+dim3$kappa[2]*cu_Ni[j_label])
          la_wei2 <- gamma(alpha*betas[j_label]+dim3$kappa[2]*cu_Ni[j_label])/gamma(alpha*betas[j_label]+dim3$kappa[2]*cu_Ni[j_label]+dim3$kappa[2])
          ir_wei <- la_wei1*la_wei2
        }
        #update theta3 like_wei        # edge(i,j,t)'s likelihood calculation
        if (dim3$datas[i_sli,j_sli,t_sli]==1){
          like_wei <- (n_tau1[i_label, j_label]+matrix(dim3$lam1,nrow=length(i_label),ncol=length(j_label)))/(n_tau[i_label, j_label]+matrix(1,nrow=length(i_label),ncol=length(j_label)))  # change the denominator with condition lam1+lam2=1 shown in crtList function
          #like_wei <- mt_pl(n_tau1[i_label, j_label],dim3$lam1)/mt_pl(n_tau[i_label, j_label],1)
        }else{
          like_wei <- (n_tau0[i_label, j_label]+matrix(dim3$lam2,nrow=length(i_label),ncol=length(j_label)))/(n_tau[i_label, j_label]+matrix(1,nrow=length(i_label),ncol=length(j_label))) 
          #like_wei <- mt_pl(n_tau0[i_label, j_label],dim3$lam2)/mt_pl(n_tau[i_label, j_label],1)
        }
        
        if (t_sli<dim3$tTime){
          like_wei <- diag(is_wei)%*%like_wei%*%diag(ir_wei)
        }
        
        like_wei<- as.vector(like_wei)/sum(like_wei)#change the matrix to vector by column and normalize
        ath_value <- 1+sum(runif(1)>cumsum(like_wei))#integer range from 1 to length like_wei
        ath_col <- ceiling(ath_value/length(i_label))
        ath_row <- ath_value-(ath_col-1)*length(i_label)
        ath_col-ath_row
        
        se_Labels[i_sli,j_sli,t_sli]<-i_label[ath_row]#new allocation
        re_Labels[i_sli,j_sli,t_sli]<-j_label[ath_col]#new allocation
        
        Nikt[i_sli, i_label[ath_row]] <- Nikt[i_sli, i_label[ath_row]]+1#add back after new allocation
        Nikt[j_sli, j_label[ath_col]] <- Nikt[j_sli, j_label[ath_col]]+1
      }
    }
    pre_Nikt <- Nikt
  }
  dim3$seLabel <- se_Labels
  dim3$reLabel <- re_Labels
  dim3$numClass <- numClass
  dim3$betas <- betas
  dim3$indexLabel <- c(1:numClass)
  return(dim3)
}