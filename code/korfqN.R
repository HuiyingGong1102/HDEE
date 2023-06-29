re_korf.qn.est <- function(gen,phe,loci){
  
  snp <- gen
  colnames(snp)<-NULL
  nn <- dim(snp)[1]
  y <- phe
  p.res <- rep(NA,50)
  
  single <- as.numeric(snp[i,])
  na.nu <- which(is.na(single))
  if(length(na.nu)>0){
    ny <- y[-na.nu,]
    nsnp <- single[-na.nu]
  }else{
    ny <- y
    nsnp <- single
  }
  st <- table(nsnp)
  h0 <- try(korf.qn.h0(ny),TRUE)
  if (class(h0)=="try-error"){
    h0<-NA
  }
  par1 <- h0[-1]
  h1 <- try(korf.qn.h1(par1,ny,nsnp),TRUE)
  if (class(h1)=="try-error"){
    h1<-NA
  }
  if(is.na(h0)||is.na(h1)){
    LR<-NA
    allpar<-c(LR,rep(NA,11))
  }else{
    LR <- 2*(h0[1]-h1[1])
    allpar<-c(LR,h0[1],h1)
  }
  cat("snp", loci, "=", allpar, "\n");
  p.res[1:length(allpar)] <- allpar
  
  return(p.res)
}




korf.qn.est <- function(gen,phe,interval=c(1,10)){
  
  snp <- gen
  colnames(snp)<-NULL
  nn <- dim(snp)[1]
  y <- phe
  n1 <- interval[1]
  n2 <- interval[2]
  if(n2 >nn)
    n2 <- nn
  p.res <- matrix(NA,nrow=length(c(n1:n2)),ncol=50)
  for(i in n1:n2){
    single <- as.numeric(snp[i,])
    na.nu <- which(is.na(single))
    if(length(na.nu)>0){
      ny <- y[-na.nu,]
      nsnp <- single[-na.nu]
    }else{
      ny <- y
      nsnp <- single
    }
    st <- table(nsnp)
    h0 <- try(korf.qn.h0(ny),TRUE)
    if (class(h0)=="try-error"){
      h0<-NA
    }
    par1 <- h0[-1]
    h1 <- try(korf.qn.h1(par1,ny,nsnp),TRUE)
    if (class(h1)=="try-error"){
      h1<-NA
    }
    if(is.na(h0)||is.na(h1)){
      LR<-NA
      allpar<-c(LR,rep(NA,11))
    }else{
      LR <- 2*(h0[1]-h1[1])
      allpar<-c(LR,h0[1],h1)
    }
    cat("snp", i, "=", allpar, "\n");
    p.res[(i-(n1-1)),1:length(allpar)] <- allpar
  }
  return(p.res)
}

korf.qn.h0<-function(ny){
  
  init.par<-c(0.8744,167.749,3.405,0,0)
  covar.par <- c(0.5,0.1)
  parin <- c(covar.par,init.par)
  
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  while(loop_k<max_iter && max_err>epsi){
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[3:7])
      AA <- curve.mle(nnpar,ny,light)
      AA
    }
    r1.covar <- optim(parin[1:2],mle.covar1,method = "BFGS",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mle(nnpar,ny,light)
      AA
    }
    r1 <- optim(c(parin[3:7]),mle.1,method = "BFGS",control=list(maxit=32000))    
    new1 <- r1$par
    
    nparin <- c(new.covar1,new1)
    newpar <- c(nparin)
    max_err <- max(abs( oldpar - newpar) );
    parin <- nparin
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  LL <- curve.mle(parin,ny,light)
  return(c(LL,parin))
}

curve.mle<-function( par,ny,light)
{
  len.cov <- 2
  par.covar <- par[1:len.cov]
  n  <- length(ny[,1])
  sigma <- AR1.get_mat(par.covar,light,1)
  
  curve.par <- par[(len.cov+1):(len.cov+ 5)]
  mu<-KE(curve.par,light)
  
  yy <- ny
  fy <- dmvnorm(yy,mu,sigma)
  A <- -sum(log(fy))
  return(A)
}

korf.qn.h1<- function(par1,ny,nsnp){
  
  snp1 <- nsnp
  par.covar <- par1[1:2]
  
  index <- table(snp1)
  snp.type <- as.numeric(names(index))
  
  g.par <- c()
  SNP.index <- list()
  for(i in 1:length(snp.type)){
    SNP.n <- which(snp1==snp.type[i])
    if(length(SNP.n)==1){
      SNP.p<-ny[SNP.n]
    } else{
      SNP.p <- c(colMeans(ny[SNP.n,],na.rm=T))
    }
    s.m <- optim(par1[3:7],s.mle4,s.y=SNP.p,s.l=light,
                 method="BFGS",control=list(maxit=32000))
    g.par <- c(g.par,s.m$par)
    SNP.index[[i]] <- SNP.n
  }
  parin <- c(par.covar,g.par)
  n.par <- length(parin)
  
  loop_k <- 1;
  max_iter <- 100;
  epsi <- 10^-4;
  max_err <- 1;
  
  
  while(loop_k<max_iter && max_err>epsi){
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      nnpar <- c(npar,parin[3:n.par])
      AA <- double.mle(nnpar,ny,light=light,SNP.index,snp.type)
      AA
    }
    r1.covar <- optim(parin[1:2],mle.covar1,method = "BFGS",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    mle1.g <- function(npar){
      nnpar <- c(new.covar1,npar)
      AA <- double.mle(nnpar,ny,light=light,SNP.index,snp.type)
      AA
    }
    # r1.g <- optim(c(parin[3:n.par]),mle1.g,method = "BFGS",control=list(maxit=32000))
    r1.g <- optim(c(parin[3:n.par]),mle1.g,method = "Nelder-Mead",control=list(maxit=32000))
    newpar <- c(new.covar1,r1.g$par)
    max_err <- max(abs( oldpar - newpar) );
    parin <- newpar
    loop_k <- loop_k+1; 
  }
  
  LL <- double.mle(parin,ny,light=light,SNP.index,snp.type)
  return(c(LL,newpar))
}

double.mle<- function(par,ny,light,SNP.index=SNP.index,snp.type=snp.type)
{
  n  <- length(ny[,1])
  
  len.cov <- 2
  par.covar <- par[1:len.cov]
  if(par.covar[1]>1||par.covar<0)
    return(NaN)
  sigma <- AR1.get_mat(par.covar,light,1)
  len.gen <- 5
  len <- 0
  A1 <- c()
  for(i in 1:length(snp.type)){
    mu.g <- par[(len.cov+len+1):(len.cov+len+len.gen)]
    mu <- KE(mu.g, light)
    yy1 <- ny[SNP.index[[i]],]
    nyy1 <- yy1
    fy1 <- dmvnorm( nyy1, mu, sigma)
    A1 <- c(A1,-sum(log(fy1)))
    len <- len + len.gen
  }
  A <- sum(A1)
  #cat("LL=",A,"\n")
  return (A);
}