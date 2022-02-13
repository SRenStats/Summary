setwd("Z:/Trade_preprocessing")

name1<-"Z:/Trade_preprocessing/Trade_Raw_Data/comtrade ("
#name2<-1
name3<-").csv"

trade_pair <- NULL
for ( i in 1:40){
  name2 <- i
  name<-paste(name1,name2,name3,sep = "")
  a <- read.csv(name)[,c(2,11,14,32)]
  trade_pair <- rbind(trade_pair,a)
}

trade_pair[,2] <- as.character( trade_pair[,2])
trade_pair[,3] <- as.character( trade_pair[,3])

country<-read.csv("gdp.csv",header = F)[,2]
year<-c(c(1992:2000),c(2017:2019))

trade_array<-array(dim=c(60,60,12))

for (k in 1:12){
  a <- which(trade_pair[,1]==year[k])
  for (i in 1:60){
    c <- which(trade_pair[,2]==country[i])
    for (j in setdiff(1:60,i)){
      r <- which(trade_pair[,3]==country[j])
      p <- intersect(intersect(a,c),r)
      if (length(p)==1)
      trade_array[j,i,k] <- trade_pair[p,4]
    }
    trade_array[i,i,k] <- 0
  }
}
write.csv(trade_array[,,10],"trade2017.csv")
