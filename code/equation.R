
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




LEind<-function(par,light)
{ 
  K<-par[1];
  a<-par[2];
  b<-par[3];
  return(K/(1+a*exp(-b*light)))  
}
MEind<-function(par,light)
{
  K<-par[1];
  a<-par[2];
  b<-par[3];
  return(K*(1-a*exp(-b*light)))
}
REind<-function(par,light)
{ 
  K<-par[1];
  a<-par[2];
  b<-par[3];
  m<-par[4];
  return(K*(1-a*exp(-b*light)^(1/1-m)))  
}

KEind<-function(par,light)
{
  K<-par[1];
  a<-par[2];
  b<-par[3];
  return(K*exp(-a/(light^b)))
}
GEind<-function(par,light)
{ 
  K<-par[1];
  a<-par[2];
  b<-par[3];
  return(K*exp(-a*exp(-b*light)))
}
BEind<-function(par,light)
{
  K<-par[1];
  a<-par[2];
  b<-par[3];
  return(K*(1-a*exp(-b*light))^3)
}
WEind<-function(par,light)
{
  K<-par[1];
  a<-par[2];
  b<-par[3];
  return(K*(1-exp(-(light/a)^b)))
}



s.mleind1 <- function(s.par,s.y,s.l){
  A <- sum((s.y - LEind(s.par,s.l))^2 )
  A
}
s.mleind2 <- function(s.par,s.y,s.l){
  A <- sum((s.y - MEind(s.par,s.l))^2 )
  A
}
s.mleind3 <- function(s.par,s.y,s.l){
  A <- sum((s.y - REind(s.par,s.l))^2 )
  A
}
s.mleind4 <- function(s.par,s.y,s.l){
  A <- sum((s.y - KEind(s.par,s.l))^2 )
  A
}
s.mleind5 <- function(s.par,s.y,s.l){
  A <- sum((s.y - GEind(s.par,s.l))^2 )
  A
}
s.mleind6 <- function(s.par,s.y,s.l){
  A <- sum((s.y - BEind(s.par,s.l))^2 )
  A
}
s.mleind7 <- function(s.par,s.y,s.l){
  A <- sum((s.y - WEind(s.par,s.l))^2 )
  A
}