library("epinet")
library(network)
library(ggplot2)
library(sna)
#library(ggnet)
library(GGally)

set.seed(1)
N <- 25
mycov <- data.frame(id = 1:N, xpos = runif(N), ypos = runif(N))
dyadCov <- BuildX(mycov, binaryCol = list(c(2, 3)),binaryFunc = "euclidean")
eta <- c(0, -7)

net <- SimulateDyadicLinearERGM(N = N, dyadiccovmat = dyadCov, eta = eta)

epi <- SEIR.simulator(M = net, N = N, beta = 1, ki = 3, thetai = 7,ke = 3, latencydist = "gamma")

plot(epi, e.col = "slategrey", i.col = "red")

net1 <- network(net)
colors <- rep(c("steelblue"),N)
colors[11]<-"red"
colors[ c(12,21)]<-"tomato"
colors[ c(13,19)]<-"pink"
colors[ c(16,25)]<-"yellow"

ggnet2(net1, size = 10, label = TRUE, label.size = 5,color = colors)

#mcmcinput <- MCMCcontrol(nsamp = 1000000, thinning = 100,etapropsd = c(1, 1))
#priors <- priorcontrol(bprior = c(0, 4), tiprior = c(1, 15),teprior = c(1, 15), etaprior = c(0, 10, 0, 10), kiprior = c(1, 7), keprior = c(1, 7), priordists = "uniform")
#out <- epinet(~ xpos.ypos.L2Dist, epidata = epi, dyadiccovmat = dyadCov,mcmcinput = mcmcinput, priors = priors)

library("epinet")
data("Hagelloch", package = "epinet")
HagellochTimes
HagellochDyadCov

