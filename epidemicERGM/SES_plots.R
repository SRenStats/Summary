######## SES ##########

# load packages
library(statnet)
library(EpiModel)
library(ndtv)

# 1) Network with SES
# initialize empty network with desired attributes
N <- 200
net <- network.initialize(n = 200, directed = FALSE)
net <- set.vertex.attribute(net, "SES", rep(0:1, each = 100))

# Setting model parameters
# NOTE: this is where we can add attributes that will define structure in the way we like based on an ERGM formula
formation <- ~edges + nodematch("SES", diff=TRUE) + nodefactor("SES") 
target.stats <- c(150, 30, 15, 20) # correspond to formation formula above

# because networks generated are dynamic, we also set parameters to determine the rate of dissolution of ties
coef.diss <- dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("SES", diff=TRUE)), 
                               duration = c(1, 10, 10))

# Fitting model with desired parameters
mod1 <- netest(net, formation, target.stats, coef.diss, edapprox = T)


# Simulate epidemic
param <- param.net(inf.prob = 0.3, act.rate = 1, rec.rate=.05)
status.vector <- rbinom(N, 1, 0.03)   # seed initial infected individuals
status.vector <- ifelse(status.vector == 1, "i", "s")
init <- init.net(status.vector = status.vector)

# Structural model features
control <- control.net(type = "SIR", nsteps = 150, nsims = 10, epi.by = "SES", ncores=2,use.pids = FALSE)

# simulate!
sim1 <- netsim(mod1, param, init, control)



# The rest is just plotting


# Plots of how the epidemic spread
par(mfrow = c(1,1), mar = c(5,5,3,3))#Down,Left,Up,Right
plot(sim1, mean.line = FALSE, qnts = FALSE, sim.lines = TRUE)

# Plots of how the simulated dynamic networks looked at different timepoints
par(mfrow = c(1,2), mar = c(0,0,1,0))
plot(sim1, type = "network", at = 1, col.status = TRUE,main = "Prevalence at t1")
plot(sim1, type = "network", at = 30, col.status = TRUE,main = "Prevalence at t30")

nw <- get_network(sim1, sim = 1)
out <- network.extract(nw, at=30)
plot(out, vertex.col="SES")


# 2. transmisson tree
tm <- get_transmat(sim1,sim=2)
par(mfrow = c(1,1), mar = c(5,5,3,3))#Down,Left,Up,Right

#par(ps = 12, cex = 1, cex.main = 1)
#par(mfrow = c(1,1), mar = c(0,0,1,0))
plot(tm,style="transmissionTimeline")



# transition network
nw1 <- as.network(tm,directed=TRUE)
nw1
names <- get.vertex.attribute(nw1,"vertex.names")
#get.vertex.attribute(nw1,"at")
colors <- ifelse(names>N/2,"tomato","steelblue")

#as.network.matrix(nw1)
ggnet2(nw1, label = "vertex.names",size = "outdegree", max_size = 10,label.size = 5,arrow.size = 10, arrow.gap = 0.02,edge.label = "at",edge.label.color = "grey10",edge.label.size = 3,color=colors)

# Phylo
#tmPhylo <- as.phylo.transmat(tm)
#plot(tmPhylo, show.node.label = TRUE,root.edge = TRUE,cex = 1)

################################3

# 1. extract networks at specific points showing race
nw <- get_network(sim1, sim = 1)#a series of networks of sim1
#sim1$network$sim1          #alternative 
#sim1[[6]]
out <- network.extract(nw, at=51)#with 300 steps
out
plot(out, vertex.col="SES")

# 2. network statistics of steps=300 networks
sim1$stats$nwstats$sim1 #300*4
sim1[[5]]

# 3. netwrok transmission matrix with information like at, sus, inf, infDur
sim1$stats$transmat$sim1 #7 columns with all the transmission information of the sim1
tm<- get_transmat(sim1, sim = 1)#the same
plot(tm,style="transmissionTimeline",show.node.label = FALSE,root.edge = TRUE,cex = 0.5)
plot(as.network(tm),show.node.label = TRUE,root.edge = TRUE,cex = 0.5)
plot(as.phylo.transmat(tm),show.node.label = TRUE,root.edge = TRUE,cex = 0.5)

mod1 <- netsim(est1, param, init, control)
tm <- get_transmat(sim1, sim = 1)
tmPhylo <- as.phylo.transmat(tm)
plot(tmPhylo, show.node.label = TRUE,
     root.edge = TRUE,
     cex = 0.5)

# 4. extract s/i/r numbers
sim1[[4]][[2]]
s1 <- sim1$epi$s.num
s2 <- sim1$epi$s.num.SES0
s3 <- sim1$epi$s.num.SES1


# Make an animated plot of the networks over a specific duration
slice.par<-list(start=1, end=30, interval=1, aggregate.dur=1,rule="latest")
compute.animation(nw,slice.par=slice.par)

render.d3movie(nw,
               filename="SES_network.html",
               edge.col="darkgray",displaylabels=TRUE,
               label.cex=.6,label.col="blue",
               output.mode = 'HTML',
               vertex.col="SES")
