rm(list=ls())
#setwd("Z:/DIM3")
setwd("C:/Users/DELL/Desktop/DIM3")
library(R.matlab)
source("DIM3_functions.R")
source("label_slice.R")
source("m_stick.R")
source("gibbs_dev.R")


###initialization####
datas<-readMat("enrondata.mat")$datas#adjacency matrix
numClass <- 4
dim3 <- crtList(datas,numClass)#initialization list with 9 dim
dim3$kappa <- c(0,0.01)
dim3$indexLabel <- c(1:numClass)
rm(numClass)
dim3$gamma <- 0.3
dim3$alpha <- 1


Niteration <- 5
###Gibbs Sampling loop####
ite_numc <- rep(0, Niteration)
deviance_numc <- rep(0, Niteration)
st_like <- rep(0, Niteration)
like_seL <- dim3$seLabel
like_reL <- dim3$reLabel
selec_like <- rep(0, Niteration)
cu_like <- rep(0, Niteration)
st_dims <- NULL


#for (n_ite in 1:Niteration)
#{
  n_ite <- 1
  dim3$betas <- rdirichlet(1,c(dim3$m_val,dim3$gamma))#sample vector beta from dirichlet distribution 
  dim3 <- label_slice(dim3)#update label,numClass
  dim3 <- m_stick(dim3)#update m
  dim3 <- Hyperpara(dim3)#update hyper parameters
  dim3 <- gibbs_dev(dim3)
  cu_jp <- dim3$cu_jointp
  cu_likes <- dim3$li_jps
  
  cu_like[n_ite] <- cu_likes
  
  spe_k <- ceiling(n_ite/Niteration*5)
  
  if (cu_likes > st_like[spe_k]){
    st_like[spe_k] <- cu_likes
    like_seL <- dim3$seLabel
    like_reL <- dim3$reLabel
    selec_like[spe_k] <- n_ite
    st_dims[[spe_k]] <- dim3#save list
  }
  
  deviance_numc[n_ite] <- dim3$deviance
  ite_numc[n_ite] <- dim3$numClass
  
}

