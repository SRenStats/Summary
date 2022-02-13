set.seed(0)
library(hergm)
data("kapferer", package = "hergm")#39 undirected
#kapferer
#gplot(kapferer, gmode = "graph", mode = "kamadakawai",displaylabels = FALSE)
object.k <- hergm(kapferer ~ edges_ij, parametric = TRUE, max_number = 2, sample_size = 1e+4,method = "bayes") #dirichlet distribution prior with 1K itrations

object.k.ml <- hergm(kapferer ~ edges_ij, max_number = 2, sample_size = 1e+4,method = "ml") #dirichlet distribution prior with 1K itrations


data("bali", package = "hergm")#17 undirected
mat <- as.matrix(bali,directed=is.directed(bali))
set.seed(0)
#get.vertex.attribute(bali,"group")
object.b5 <- hergm(bali ~ edges_ij + triangle_ijk, max_number = 5, variational=T, method="bayes",sample_size = 1e+4)#default is ml
object.b2 <- hergm(bali ~ edges_ij + triangle_ijk, max_number = 2, variational=T, method="bayes",sample_size = 1e+4)#default is ml

#try to plot the traceplot to se more details
object <- object.b2
View(object$hergm_theta)
par(mfrow=c(1,1));plot(object$hergm_theta[,5],type = "l")
par(mfrow=c(3,3))
for (i in setdiff(1:12,c(4,8,12))){
  plot(object$hergm_theta[,i],type = "l")
}
object$mcmc.diagnostics


#object.b5$mcmc.diagnostics
#object.b$parameterization
gf2 <- gof(object.b2)
table(object.b$indicator[8000,],get.vertex.attribute(bali,"group"))

object.b.ml <- hergm(bali ~ edges_ij + triangle_ijk, max_number = 2, method="ml",sample_size = 1e+4)#default is ml

library(hergm)
data("sampson") #32 undirected
object.s <- hergm(samplike ~ edges_ij + mutual_ij + ttriple_ijk, max_number = 3,method="bayes",sample_size = 1e+4)
object.s 
table(object.s$indicator[8000,],get.vertex.attribute(samplike,"group"))

object <- object.s
View(object$hergm_theta)
par(mfrow=c(1,1));plot(object$hergm_theta[,5],type = "l")
par(mfrow=c(3,3))
for (i in setdiff(1:12,c(4,8,12))){
  plot(object$hergm_theta[,i],type = "l")
}
object$mcmc.diagnostics



object.s.ml <- hergm(samplike ~ edges_ij + mutual_ij + ttriple_ijk, max_number = 3,method="ml",sample_size = 1e+4)

setwd("C:/Users/DELL/Desktop/Node_ERGM")
mat <- as.matrix(read.table("ELadv.dat"))
net<-network(mat)
atts<-read.csv("ELattr.csv")
atts[,5] <- atts[,5]/10
atts[,6] <- atts[,6]/100
for (i in 1:ncol(atts))
  set.vertex.attribute(net,names(atts)[i],atts[,i])
lawyer <- mat
object.l.ml <- hergm(lawyer~ edges_ij + mutual_ij + ttriple_ijk, max_number = 2,method="ml",sample_size = 1e+4)

#####Assessment
mat <- as.matrix(bali);n<-nrow(mat)
terms<-c('edges','triangle');Nterms <- length(terms);Nterms_nd <- Nterms-1
ergmformula <- paste("~", paste(terms,collapse="+"),sep="")
z <- rep(2,n);z[c(2,7,10:14,16)]<-1;z;table(z)
theta_hat <- matrix(nrow=2,ncol=2)
theta_hat[,1]<-c(-1.295,0.682)
theta_hat[,2]<-c(-0.671,1.501)
theta_hat

i<-2
net_sub <- as.network(mat[which(z==i),which(z==i)],directed = F)
simulate(as.formula(paste("net_sub",ergmformula)),nsim=20,coef=theta_hat[,i],output='stats')
summary(as.formula(paste("net_sub",ergmformula)))

#
mat <- as.matrix(bali);n<-nrow(mat)
terms<-c('edges','triangle');Nterms <- length(terms);Nterms_nd <- Nterms-1
ergmformula <- paste("~", paste(terms,collapse="+"),sep="")
z <- rep(2,n);z[c(2,7,10:14,16)]<-1;z;table(z)
theta_hat <- matrix(nrow=2,ncol=2)
theta_hat[,1]<-c(-1.295,0.682)
theta_hat[,2]<-c(-0.671,1.501)
theta_hat

i<-2
net_sub <- as.network(mat[which(z==i),which(z==i)],directed = F)
simulate(as.formula(paste("net_sub",ergmformula)),nsim=20,coef=theta_hat[,i],output='stats')
summary(as.formula(paste("net_sub",ergmformula)))

