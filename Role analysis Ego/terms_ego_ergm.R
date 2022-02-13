ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
obs.S<-list()
print("Calculating all network statistics...")
for (i in 1:N)
  {
  #obs.S[[i]]<-ergmMPLE(as.formula(paste("x[[i]]",ergmformula)),output="array") # pseudolikelihood change stats
  obs.S[[i]]<-ergmMPLE(as.formula(paste("x[[i]]",ergmformula))) # pseudolikelihood change stats
  obs.S[[i]]$offset<-ergm.offset[i]
  }

