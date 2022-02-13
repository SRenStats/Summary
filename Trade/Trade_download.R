##对数据重新排序
transform<-function(data1){
  data2<-matrix(nrow=n+1,ncol=T+1)
  data2[1,-1]<-seq(2001,2017,1)
  data2[-1,1]<-country2
  for (t in 2:ncol(data2))
  {
    for (j in 2:nrow(data2))
    {
      c<-which(data1[,1]==data2[j,1])
      d<-which(data1[1,]==data2[1,t])
      if (length(c)!=0&length(d)!=0)
        data2[j,t]<-data1[c,d]
    }
  }
  return(data2)
}
country1<-read.csv("C:/Users/Thinkpad/Desktop/country1.csv",header=FALSE)
country1<-as.matrix(country1)
country2<-read.csv("C:/Users/Thinkpad/Desktop/country2.csv",header=FALSE)
country2<-as.matrix(country2)

n<-240
T<-17
matrix<-array(dim=c(n,n,T))

setwd("~/tradetxt")
getwd()
name1<-"Trade_Map_-_List_of_supplying_markets_for_a_product_imported_by_"
name3<-".txt"


i<-i+1
name2<-country1[i]
name2
name<-paste(name1,name2,name3,sep = "")
x<-read.table(name,header=FALSE,skip=1)
x<-as.matrix(x)
lable0<-read.table(name,header=FALSE,encoding = "UTF-8",nrows=1)
lable<-lable0[seq(1,length(lable0),4)]
lable[1]<-"Exporters"
lable<-as.matrix(lable)
data1<-rbind(lable,x)
data2<-transform(data1)
for (col in 1:T)
{
  matrix[,i,col]<-data2[-1,col+1]
}


i<-i-1
country1[i]
matrix[,i,]<-matrix[,1,]






m<-rbind(matrix1[,,1],matrix1[,,2],matrix1[,,3],matrix[,,4],matrix1[,,5],
         matrix1[,,6],matrix1[,,7],matrix1[,,8],matrix[,,9],matrix1[,,10],
         matrix1[,,11],matrix1[,,12],matrix1[,,13],matrix1[,,14],matrix1[,,15],
         matrix1[,,16])
write.csv(m,file="C:/Users/Thinkpad/Desktop/trade.csv")

bkx<-read.csv("C:/Users/Thinkpad/Desktop/bkx.csv",header=FALSE)
bkxa<-intersect(bkx[,1],bkx[,2])
names1<-country1[bkxa]

Ts<-seq(1:16)
T1<-length(Ts)
n1<-length(bkxa)
data3<-array(dim=c(n1,T1,n1))
matrix1<-matrix[bkxa,bkxa,Ts]
for (i in 1:n1)
{
  matrix1[i,i,]<-0
}

for (k in 1:n1)
{
  for (j in Ts)
  {
    data3[,j,k]<-matrix1[,k,j]
  }
}

qs<-c()
for (i in 1:n1)
{
  qs[i]<-sum(is.na(data3[,,i]))
}

look<-cbind(order(qs,decreasing = T),names1[order(qs,decreasing = T)],qs[order(qs,decreasing = T)])
bkxa2<-order(qs,decreasing = T)[58:118]
bkxa2<-bkxa2[-17]
names2<-names1[bkxa2]

n2<-length(bkxa2)
data4<-array(dim=c(n2,T1,n2))
matrix2<-matrix1[bkxa2,bkxa2,]

for (k in 1:n2)
{
  for (j in Ts)
  {
    data4[,j,k]<-matrix2[,k,j]
  }
}

qs2<-c()
for (i in 1:n2)
{
  qs2[i]<-sum(is.na(data4[,,i]))
}


look2<-cbind(order(qs2,decreasing = T),names2[order(qs2,decreasing = T)],qs2[order(qs2,decreasing = T)])
m2<-rbind(matrix2[,,1],matrix2[,,2],matrix2[,,3],matrix2[,,4],matrix2[,,5],
          matrix2[,,6],matrix2[,,7],matrix2[,,8],matrix2[,,9],matrix2[,,10],
          matrix2[,,11],matrix2[,,12],matrix2[,,13],matrix2[,,14],matrix2[,,15],
          matrix2[,,16])
colnames(m2)<-names2
write.csv(m2,file="C:/Users/Thinkpad/Desktop/trade61.csv")








