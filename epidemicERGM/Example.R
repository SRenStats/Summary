library(statnet)
library(EpiModel)
library(ndtv)
set.seed(10)
nw <- network.initialize(n = 100, directed = FALSE)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
param <- param.net(inf.prob = 0.5)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 40, nsims = 10, verbose = FALSE,
                       use.pids = FALSE)
sim1 <- netsim(est1, param, init, control)

###########PLOTS
par(mfrow = c(1,1), mar = c(5,5,3,3))#Down,Left,Up,Right
plot(sim1, mean.line = FALSE, qnts = FALSE, sim.lines = TRUE)

tm <- get_transmat(mod1)
tmPhylo <- as.phylo.transmat(tm)
plot(tmPhylo, show.node.label = TRUE,
     root.edge = TRUE,
     cex = 0.5)


nw <- get_network(mod1, sim = 1)#a series of networks of sim1
out <- network.extract(nw, at=31)#with 300 steps
out
# 0. plot netsim
# mod <- mod1
# plot(mod)
# plot(mod, type = "epi", grid = TRUE)
# plot(mod, type = "epi", popfrac = TRUE)
# plot(mod, type = "epi", y = "si.flow", qnts = 1, ylim = c(0, 4))
# 1. compartmental plots 
comp_plot(mod1,at=25)

# 2. transmisson tree
tm <- get_transmat(mod1)
plot(tm,style="transmissionTimeline")



# transition network
nw1 <- as.network(tm,directed=TRUE)
ggnet2(nw1, label = TRUE,size = "outdegree", max_size = 10,label.size = 5,arrow.size = 10, arrow.gap = 0.025)

# Phylo
tmPhylo <- as.phylo.transmat(tm)
plot(tmPhylo, show.node.label = TRUE,root.edge = TRUE,cex = 1)
