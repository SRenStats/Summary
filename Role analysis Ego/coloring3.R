# 3-colouring algorithm
#   using this: 
#   A<-matrix(0,10,10)
#   A[1,1:4]<-A[2,1:4]<-A[3,c(1,2,5)]<-A[4,c(1,2,5)]<-A[5,3:6]<-A[6,5:8]<-A[7,c(6,9,10)]<-A[8,c(6,9,10)]<-A[9,7:10]<-A[10,7:10]<-1;diag(A)<-1
#   gives the same partitioning as the paper
A<-as.sociomatrix(net)
n<-network.size(net)
dd=sum(A)/(2*n)
eigA<-eigen(A)
sel<-c(which.min(abs(eigA$values-(2*dd))),which.min(abs(eigA$values-(dd))),which.min(abs(eigA$values+(dd))))
P<-t(eigA$vectors[,sel])
S<-t(P)%*%P

R3<-cutree(hclust(dist(S),method="ward"),k=3)

#dists<-as.matrix(dist(S))
#R<-1:n # each it's own class at first
#diag(dists)<-max(dists)
#while (length(unique(R))>3)
#  {
#  tmp<-which(dists==min(dists),arr.ind=1)[1,] # find the closest pair
#  ord<-c(which.min(R[tmp]),which.max(R[tmp]))
#  R[tmp[ord[2]]]<-R[tmp[ord[1]]] # join the two closest
#  dists[tmp[1],tmp[2]]<-max(dists);dists[tmp[2],tmp[1]]<-max(dists);
#  }
#R[which(R==min(R))]<-1
#R[which(R==unique(R)[2])]<-2
#R[which(R==max(R))]<-3
