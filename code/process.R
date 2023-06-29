getlogp1 <- function(ret,gen,npar){
  SNP <- gen
  SNP<-as.matrix(SNP)
  LR <- ret[,1]
  LR[which(is.na(LR))] <- 1
  snpnames <- 1:dim(SNP)[1]
  testcross <- c()
  intercross <- c()
  for(i in 1:length(LR)){
    index <- which(is.na(SNP[i,]))
    if(length(index)>0){
      nsnp <- SNP[i,-index]
    }else{
      nsnp <- SNP[i,]
    }
    tt <- table(nsnp)
    if(length(tt)==2){
      testcross<- c(testcross,i)
    }
    if(length(tt)==3){
      intercross<- c(intercross,i)
    }
  }
  df <- rep(NA,length(LR))
  df[testcross] <- npar
  df[intercross] <- 2*npar
  p.value <- pchisq(LR,df,lower.tail = F)
  logp <- -log10(p.value)
  # fdrp <- p.adjust(p.value,method="fdr",n=length(p.value))
  # logp <- -log10(fdrp)
  return (logp)
}
getlogp2 <- function(ret,gen){
  SNP <- gen
  SNP<-as.matrix(SNP)
  LR <- ret[,1]
  LR[which(is.na(LR))] <- 1
  snpnames <- 1:dim(SNP)[1]
  testcross <- c()
  intercross <- c()
  for(i in 1:length(LR)){
    index <- which(is.na(SNP[i,]))
    if(length(index)>0){
      nsnp <- SNP[i,-index]
    }else{
      nsnp <- SNP[i,]
    }
    tt <- table(nsnp)
    if(length(tt)==2){
      testcross<- c(testcross,i)
    }
    if(length(tt)==3){
      intercross<- c(intercross,i)
    }
  }
  df <- rep(NA,length(LR))
  df[testcross] <- 6
  df[intercross] <- 12
  p.value <- pchisq(LR,df,lower.tail = F)
  logp <- -log10(p.value)
  return (logp)
}



getp <- function(ret,gen,npar){
  SNP <- gen
  SNP<-as.matrix(SNP)
  LR <- ret[,1]
  LR[which(is.na(LR))] <- 1
  snpnames <- 1:dim(SNP)[1]
  testcross <- c()
  intercross <- c()
  for(i in 1:length(LR)){
    index <- which(is.na(SNP[i,]))
    if(length(index)>0){
      nsnp <- SNP[i,-index]
    }else{
      nsnp <- SNP[i,]
    }
    tt <- table(nsnp)
    if(length(tt)==2){
      testcross<- c(testcross,i)
    }
    if(length(tt)==3){
      intercross<- c(intercross,i)
    }
  }
  df <- rep(NA,length(LR))
  df[testcross] <- npar
  df[intercross] <- 2*npar
  p.value <- pchisq(LR,df,lower.tail = F)
  return (p.value)
}
