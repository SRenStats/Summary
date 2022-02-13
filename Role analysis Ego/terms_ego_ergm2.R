ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
obs.S<-list()#list of every ego-network with 4 objects
print("Calculating all network statistics...")
for (i in 1:N)
  {
  #obs.S[[i]]<-ergmMPLE(as.formula(paste("x[[i]]",ergmformula)),output="array") # pseudolikelihood stats
  obs.S[[i]]<-ergmMPLE(as.formula(paste("sim.x[[i]]",ergmformula))) # pseudolikelihood statistics
  obs.S[[i]]$offset<-init$ergm.offset[i]
  }

