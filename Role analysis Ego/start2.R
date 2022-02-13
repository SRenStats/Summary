require(boot) 
require(ergm) 
require(sna) 
LOWESTLL=-1e8 # lowest allowed loglikelihood to avoid numerical instabilities
tol=1e-6 # for convergence checking

if (START=="READ") # read in a network and extract ego-networks as x
  {
  N=net$gal$n
  directed=is.directed(net)
  Y<-as.matrix.network(net)
  ### find each ego-network; use K steps out from each node
  # use Y+t(Y) to jump backwards across links at the second step out
  #Y2<-Y+t(Y)
  #net2<-network(Y2)
  #x_in<-gapply(net2,c(1,2),1:N,"*",1,distance=K) 
  x_in<-gapply(net,1,1:N,"*",1,distance=K) 
    x<-x_in
  for (i in 1:N)
    {
    x_in[[i]]<-c(i,x_in[[i]])
    x[[i]]<-as.matrix(Y[x_in[[i]],x_in[[i]]])
    }
  x<-lapply(x,network,directed=directed)
  if (length(Covariate)!=0)
  {
    for (i in 1:N)
    {
      set.vertex.attribute(x[[i]],"gdp",gdp_value[x_in[[i]]])
    }
  }
  rm(x_in)
  # only look at egos bigger than MINSIZE
  if (MINSIZE>1)
    {
    keep<-lapply(x,network.size)>=MINSIZE
    all_net<-net
    N=all_net$gal$n
    net<-network(as.sociomatrix(all_net)[(1:N)[keep],(1:N)[keep]],directed=is.directed(all_net))
    for (att in list.vertex.attributes(all_net))
      set.vertex.attribute(net,att,get.vertex.attribute(all_net,att)[keep])
    rm(all_net)
    x<-x[keep]
    }
  N=length(x)
  }
if (START=="SIM")
  {
  #### simulate new data
  source("sim_networks2.R")
  #directed=TRUE
  x<-sim.x
  rm(sim.x)
  }
if (START=="RSIM")
  {
  #### read simulated data
  START="SIM"
  x<-sim.x
  N=length(x)
  }
