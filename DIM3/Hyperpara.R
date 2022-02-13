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