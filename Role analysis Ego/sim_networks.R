# code to (approximately) simulate from the EgoERGM model
# we want ERGMS simulated from G groups with distinct ERGM parameters. We then want the these ego
# networks to come together as a single network. One way is to have all egonets as disconnected
# components...

require(statnet)
### create the underlying parameters
sim.theta<-rbind(c(-3,1,0), c(-1,-2,-1), c(-2,0,2)) # easy version
#sim.theta<-rbind(c(-0.05,0.05,0), c(-0.05,-0.05,-0.05), c(-0.05,0,0.05)) # hard version
sim.ego.terms<-c("edges", "gwesp(0.8,fixed=TRUE)", "gwdegree(decay=0.8,fixed=TRUE)")
G=nrow(sim.theta)
M=50 # M egonets per group
N=G*M 
### assign the egos to groups
sim.K<-rep(NaN,N)
for (g in 1:G)
  sim.K[((g-1)*M+1):(g*M)]<-rep(g,M)
### simulate the egos
NN=20 # number of nodes per ego
ergmformula <- paste("~", paste(sim.ego.terms,collapse="+"),sep="")
sim.x<-list()
for (i in 1:N)
  sim.x[[i]]<-simulate(as.formula(paste("network(NN,directed=FALSE)",ergmformula)), coef=c(sim.theta[(1:G)[sim.K[i]],]))

save(sim.ego.terms,sim.theta,sim.K,sim.x, file="simulated_networks.rdata")

