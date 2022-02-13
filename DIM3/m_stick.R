m_stick<-function(dim3){
  #library(gmp)
  #sampling m's value
  #dim3 is the structure used
  indexLabel <- dim3$indexLabel
  betas <- dim3$betas
  alpha <- dim3$alpha
  
  #method 1 to remove empty clusters to adjust the right class number
  #empty_clu <- NULL
  #for (i_m  in 1:dim3$numClass){
  #if ((sum(dim3$seLabel == indexLabel[i_m])==0)&&(sum(dim3$reLabel == indexLabel[i_m],na.rm=TRUE)==0))
  #empty_clu <- c(empty_clu, i_m)
  #}
  #if (length(empty_clu)!=0){
  #  indexLabel <- indexLabel[-empty_clu] 
  #  dim3$numClass <- dim3$numClass-length(empty_clu)
  #  betas <- betas[-empty_clu] 
  #}
  
  #method 2 to remove empty label
  actual1_clu <- unique(c(as.vector(dim3$reLabel),as.vector(dim3$seLabel)))
  actual2_clu <- sort(actual1_clu[!is.na(actual1_clu)])#ordered nonempty label
  indexLabel <- indexLabel[actual2_clu] 
  dim3$numClass <- length(actual2_clu)
  betas <- betas[actual2_clu] 
  
  
  
  #get the m_val value
  m_val <- rep(0, dim3$numClass)
  nohat_m <- rep(0, dim3$numClass)
  pre_Nikt <- matrix(0,nrow=dim3$dataNum, ncol=dim3$numClass)
  Nikt <- matrix(0,nrow=dim3$dataNum, ncol=dim3$numClass)
  for (t_m in 1:dim3$tTime){
    #t_m<-1
    for (i_m in 1:dim3$dataNum){
      Nikt[i_m,] <- count(c(dim3$seLabel[i_m,,t_m],dim3$reLabel[,i_m,t_m]),indexLabel)
    }
    for (i_m in 1:dim3$dataNum){
      #i_m <-1
      t_Table <- rep(0, dim3$numClass)
      nohat_t <- rep(0, dim3$numClass)
      for (k_m in 1:dim3$numClass){
        #k_m<-1
        i_max <- Nikt[i_m,k_m]
        if (i_max > 0)
        {
          i_stir <- stirling(i_max)*((alpha*betas[k_m]+dim3$kappa[2]*pre_Nikt[i_m, k_m])^(1:i_max))#prb
          i_stir <- i_stir/sum(i_stir)
          i_nt <- 1+sum(runif(1) > cumsum(i_stir))
          i_p <- (dim3$kappa[2]*pre_Nikt[i_m, k_m])/(alpha*betas[k_m]+dim3$kappa[2]*pre_Nikt[i_m, k_m])
          #i_p <- (alpha*betas[k_m])/(alpha*betas[k_m]+dim3$kappa[2]*pre_Nikt[i_m, k_m]/(2*dim3$dataNum))
          i_val <- rbinom(1,i_nt, i_p)
          t_Table[k_m] <- i_nt-i_val
          nohat_t[k_m] <- i_nt
        }
      }
      m_val <- m_val + t_Table
      nohat_m <- nohat_m + nohat_t
      
    }
    pre_Nikt <- Nikt
  }
  
  
  # disp(m_val);
  
  dim3$nohat_m <- nohat_m
  dim3$indexLabel <- indexLabel
  dim3$m_val <- m_val
  dim3$betas <- betas
  
  return(dim3)
}

