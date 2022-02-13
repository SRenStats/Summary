
log_lp<-function(s,p,theta_initial,epsilon,eps,lam,ep=1e-4,iter){ 
            #s is the input sample covariance matrix and theta_initial is the initial estimates of the precision matrix
            #p is the value between 0 and 1			
			if(!is.matrix(s)) 
                stop("s has to be a matrix")
			if(nrow(s)!=ncol(s))
		        stop("s must be a square matrix")
			if(!is.matrix(theta_initial)) 
                stop("theta_initial has to be a matrix")
			# if(dim(s)!=dim(theta_initial))
                # stop			
			if(!is.vector(epsilon))
		        stop("epsilon must be a vector")
		    if(length(epsilon)!=ncol(s))
		        stop
			if(!(p<=1&p>0))
			    stop
		    if(lam<0)
			     stop
			theta<-theta_initial
			W<-matrix(0,nrow(s),ncol(s))#the weight of each reweighted l1 regularization
			k<-1
			repeat{
			     theta_old<-theta
			     for(l in 1:(nrow(W)-1)){
		         for(r in (l+1):ncol(W)){
	 W[l,r]<-(1/(sum((abs(theta[,l]))^p)+epsilon[l])+1/(sum((abs(theta[,r]))^p)+epsilon[r]))*lam*p/((abs(theta_initial[l,r]))^(1-p)+eps)
			     W[r,l]<-W[l,r]}}
                 rho<-W
				  
				 theta<-glasso(s,rho,penalize.diagonal=FALSE)$wi
				 
				 
				 
				 
				 if(norm(theta_old-theta,"F")<=ep|k>=iter)
			     break
				 k<-k+1
					}
			 list(theta=theta,iter=k)
			 }





data<-read.csv("C:/Users/Zzz/Desktop/all2015.csv",header=TRUE)
#1:31  2:31 2:32 3:32 3:33 4:33 4:34 5:34 5:35 6:35 6:36 7:36 7:37 8:37 8:38 9:38 9:39 10:39 10:40 11:40 11:41 12:41 12:42 13:42 13:43 14:43 14:44 15:44 15:45 16:45 16:46 
x<-as.matrix(data)

   # for(j in 1:31){
         # x[,j]<-x[,j]-mean(x[,j])}
    # for(j in 1:31){
        # x[,j]<-x[,j]/sd(x[,j])}
n<-nrow(x)
	  s<-(1/n)*t(x)%*%x

lambda<-seq(2000,4000,1)
#lambda<-seq(0.01,1,0.01)
		  m<-length(lambda)
		  bic<-rep(0,m)
		  edge<-matrix(0,nrow(s),ncol(s))
		  for(j in 1:m){
		      rho<-lambda[j]
		      theta<-glasso(s,rho,penalize.diagonal=FALSE)$wi
			  for(k in 1:(nrow(theta)-1)){
			     for(r in (k+1):ncol(theta)){
				  if(theta[k,r]==0){
				  edge[k,r]<-0}
				  else{
				  edge[k,r]<-1}
				  }}
			     
			 bic[j]<-sum(edge)}
			   
			   index<-min(which(bic==120))
			   rho<-lambda[index]
			   
			   fit<-glasso(s,rho,penalize.diagonal=FALSE)$wi
p<-0.5
theta_initial<-fit
epsilon<-rep(0.1,ncol(x))
eps<-0.1
iter<-2

lambda<-seq(1,500,0.2)#1
		 #lambda<-seq(0.01,1,0.01)#0.5
		  
		  
		  m<-length(lambda)
		  bic<-rep(0,m)
		  edge<-matrix(0,nrow(s),ncol(s))
		  for(j in 1:m){
		      lam<-lambda[j]
		      theta<-log_lp(s,p,theta_initial,epsilon,eps,lam,ep=1e-3,iter)$theta
			  for(k in 1:(nrow(theta)-1)){
			     for(r in (k+1):ncol(theta)){
				  if(theta[k,r]==0){
				  edge[k,r]<-0}
				  else{
				  edge[k,r]<-1}
				  }}
			    
				bic[j]<-sum(edge) }
			
			   index<-min(which(bic==50))
			   lam1<-lambda[index]
			  
			   fit<-log_lp(s,p,theta_initial,epsilon,eps,lam=lam1,ep=1e-3,iter)


theta_esti_0<-fit$theta
th0<-A(theta_esti_0)
adjm<-th0
data<- read.csv("C:/Users/Zzz/Desktop/API(1月).csv",header=TRUE)
b<-colnames(data)

 a<-NULL
 for(i in 1:nrow(adjm)){
 if(any(adjm[i,])!=0)
a<-cbind(a,i)
 }
 adjm<-adjm[a,a]
 rownames(adjm)<-b[a]
 colnames(adjm)<-b[a]

g1 <- graph.adjacency( adjm,mode="undirected")
E(g1)$color<-"black"
V(g1)$color<-"green"
#par(mfcol=c(2,1))
degree<-degree(g1)
V(g1)[14]$color<-"darkorchid1"
 V(g1)[c(12,15,16,17)]$color<-"orangered2"
#vertex.lable<-a[-1]
 plot.igraph(g1,layout=layout.fruchterman.reingold,vertex.size=22,vertex.label.cex=.8,main="2-32 log_l1/2")

library(showtext)
 showtext.auto(enable = TRUE)
 font.add('kaishu','simkai.ttf')
 showtext.begin()
 pdf('hub2014.pdf')
  plot.igraph(g1,layout=layout.fruchterman.reingold,vertex.size=22,vertex.label.cex=.9,main="",vertex.label.family="kaishu")

showtext.end()
 dev.off()









