library("deSolve")
library("mvtnorm")
#load("gentable.Rdata")
y1 <- as.matrix(phetable[,1:11]/100)
y2 <- as.matrix(phetable[,12:22])
y3 <- as.matrix(phetable[,23:33])

LE<-function(par,light)
{ 
  K<-par[1];
  a<-par[2];
  b<-par[3];
  c<-par[4];
  d<-par[5];
  return(K/(1+a*exp(-b*light))+c*light^d)  
}
ME<-function(par,light)
{
  K<-par[1];
  a<-par[2];
  b<-par[3];
  c<-par[4];
  d<-par[5];
  return(K*(1-a*exp(-b*light))+c*light^d)
}

KE<-function(par,light)
{
  K<-par[1];
  a<-par[2];
  b<-par[3];
  c<-par[4];
  d<-par[5];
  return(K*exp(-a/(light^b))+c*light^d)
}

s.mle1 <- function(s.par,s.y,s.l){
  A <- sum((s.y - LE(s.par,s.l))^2 )
  A
}
s.mle2 <- function(s.par,s.y,s.l){
  A <- sum((s.y - ME(s.par,s.l))^2 )
  A
}

s.mle4 <- function(s.par,s.y,s.l){
  A <- sum((s.y - KE(s.par,s.l))^2 )
  A
}


eret<-logi.etr.est(genetr,y1,c(1,length(genetr[,1])))
pret<-mits.qp.est(genqp,y2,c(1,length(genqp[,1])))
nret<-korf.qn.est(genqn,y3,c(1,length(genqn[,1])))
