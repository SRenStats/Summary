# code to perform regular equivalence partitioning of a network
require(network)
require(blockmodeling)
if (!exists("R2")) R2=2
D<-REGE.for(M=as.sociomatrix(net))$E
hmethods<-c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
if (!exists("h")) h=1
clu<-cutree(hclust(d=as.dist(1-D),method=hmethods[h]),k=G)
par(mar=rep(2,4))
plot(net,vertex.col=clu,vertex.cex=R2,main=paste("using method", hmethods[h]))

