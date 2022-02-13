#ergm(net~edges+gwidegree(decay=0.8,fixed=TRUE)+gwodegree(decay=0.8,fixed=TRUE))
ergmMPLE(net~edges+mutual+nodeicov("gdp")+nodeocov("gdp"),output="fit")$coef


#ego.terms<-c("edges", "mutual", "nodeicov("gdp")","nodeocov("gdp")")


ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")

ergmformula<-c(as.character(edges+mutual+nodeicov("gdp")+nodeocov("gdp")))


ego.terms<-c('edges', 'mutual', 'nodeicov("gdp")','nodeocov("gdp")')

ergmformula <- c("~edges+mutual+nodeicov("gdp")+nodeocov("gdp")")

ego.terms<-c("edges", "mutual")
ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
form<-nonsimp_update.formula(as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)
Nterms<-length(ego.terms)
theta<-matrix(0,N,Nterms)
for (i in 1:N)
  theta[i,]<-ergmMPLE(form,output="fit")$coef 


ego.terms<-c('edges', 'mutual', 'nodeicov("gdp")','nodeocov("gdp")')

#ego.terms<-c('edges', 'mutual')
ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep='')
form<-nonsimp_update.formula(as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)
Nterms<-length(ego.terms)
theta<-matrix(0,N,Nterms)
for (i in 1:2)
  theta[i,]<-ergmMPLE(form,output="fit")$coef 



