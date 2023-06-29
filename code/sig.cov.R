
cov.test <- function(pheno,snp){
  y <- pheno
  snps <- snp
  nsnp <- dim(snps)[1]
  time <- dim(y)[2]
  p.matrix <- matrix(rep(0,nsnp*time),nrow=time)
  fdr.matrix <- matrix(rep(0,nsnp*time),nrow=time)
  for(i in 1:time){
    for(j in 1:nsnp){
      tmpsnp <- snps[j,]
      miss <- which(is.na(tmpsnp))
      if(length(miss)>0){
        newsnp <- tmpsnp[-miss]
        yy <- y[-miss,i]
      }else{
        newsnp <- tmpsnp
        yy <- y[,i]
      }
      symbol <- names(table(as.character(unlist(c(newsnp)))))
      index1 <- which(newsnp==symbol[1])
      y11 <- yy[index1]
      index0 <- which(newsnp==symbol[2])
      y0 <- yy[index0]
      if(length(index1)==1||length(index0)==1){
        p.value <- 1
      }else{
        var.t <- var.test(y0,y11)
        if(var.t$p.value>0.05)
          var.i <- TRUE
        else
          var.i <- FALSE
        p.value <- t.test(y0,y11,var.equal = var.i)$p.value#summary(snp.aov)[[1]][[1,"Pr(>F)"]] 
      }
      p.matrix[i,j] <- p.value
    }
    fdr.matrix[i,] <- p.adjust(p.matrix[i,],method="fdr")      
  }
  
  colnames(p.matrix)<- rownames(snps)
  colnames(fdr.matrix)<- rownames(snps)
  
  logp <- -log10(p.matrix)
  thre1 <- -log10(0.05/nsnp)
  thre2 <- -log10(0.01/nsnp)
  return(list(p.value=p.matrix,fdr=fdr.matrix,logp=logp,thre1=thre1,thre2=thre2))
}

etr.sig.time <- cov.test(pheno=y1,snp=genetr)
nsig1 <- c()
for (i in 1:11) {
  nsig1 <- c(nsig1,length(which(etr.sig.time$logp[i,]>etr.sig.time$thre1)))
}

for (ii in 4:11) {
  ETRdat<-cbind(as.numeric(rownames(genetr)),snpinfoetr$CHR,snpinfoetr$POS,snpinfoetr$REF,snpinfoetr$ALT,etr.sig.time$logp[ii,])
  ETRdat<-as.data.frame(ETRdat)
  colnames(ETRdat)<-c("seq","CHR","POS","REF","ALT","logp")
  ETRdat$seq<-as.numeric(ETRdat$seq)
  ETRdat$POS<-as.numeric(ETRdat$POS)
  ETRdat$logp<-as.numeric(ETRdat$logp)
  ETRdat<-ETRdat[order(ETRdat[,1]),]
  rownames(ETRdat)<-NULL
  scaffoldsum <- length(table(ETRdat[,2]))
  chrnum<-c(rep(1,1395),rep(2,440),rep(3,503),rep(4,603),rep(5,696),rep(6,713),rep(7,437),rep(8,454),rep(9,198),rep(10,433),rep(11,408),
            rep(12,383),rep(13,298),rep(14,644),rep(15,423),rep(16,664),rep(17,394),rep(18,330),rep(19,375),rep(20,19),rep(21,4),rep(22,7),
            rep(23,9),rep(24,3),rep(25,6),rep(26,2),rep(27,4),rep(29,9),rep(30,1),rep(31,4),rep(32,4),rep(33,7),rep(34,2),rep(35,8),rep(36,6),
            rep(37,4),rep(38,7),rep(39,5),rep(41,2),rep(42,7),rep(45,1),rep(47,1),rep(49,4),rep(50,3),rep(51,1),rep(52,26),rep(56,1),rep(61,3),
            rep(63,3), rep(64,1),rep(65,2),rep(66,1),rep(69,2),rep(76,1),rep(78,1),rep(79,1),rep(80,1),rep(82,8),rep(86,3),rep(87,1),rep(92,1),
            rep(93,2),rep(102,1),rep(104,1),rep(107,3),rep(109,2),rep(111,5),rep(112,2),rep(116,1),rep(119,3),rep(121,1),rep(123,2))
  ETRdat$CHR<-chrnum
  scaffoldc <- as.numeric(names(table(ETRdat[,2])))
  for(i in 1:scaffoldsum){
    ETRdat[which(ETRdat[,2]==scaffoldc[i]),2] <- i
  }
  for (i in 1:scaffoldsum){ 
    ndx <- which(ETRdat[, 2]==i)
    lstMrk <-max(ETRdat[ndx, 3])
    if (i < scaffoldsum) ndx2 <-which(ETRdat[, 2]==i+1)
    if (i < scaffoldsum) ETRdat[ndx2, 3]<-ETRdat[ndx2, 3]+lstMrk
  }
  
  bpMidVec <- vector(length=scaffoldsum)
  
  for (i in 1:scaffoldsum){
    ndx <- which(ETRdat[, 2]==i)
    posSub <- ETRdat[ndx, 3]
    bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
  }
  
  library(ggplot2)
  p1<- ggplot()
  p1<-p1+ geom_point(data=ETRdat,aes(x=POS, y=logp,colour=as.factor(as.numeric(CHR))),size=2,alpha=0.8)
  p1<-p1+scale_color_manual(values=rep(c('#FED439FF', '#D2AF81FF'), 36))
  # p1<-p1+geom_hline(yintercept=20,linetype=1, col='black', lwd=1)
  p1<-p1+geom_hline(yintercept=etr.sig.time$thre1,linetype=1, col='red', lwd=1)
  p1<-p1+scale_x_continuous(breaks=c(bpMidVec[1:19]),labels=c(1:19))
  # p1<-p1+xlab('Chromosome') + ylab('-log10 (p-value)')+theme_bw()+theme(legend.position = "none")
  p1<-p1+xlab('') + ylab('')+theme_bw()+
    theme(legend.position = "none",
          axis.text.x=element_text(hjust=0.5,size=8,color='black'),
          axis.text.y=element_text(vjust=0.2,size=10,color='black'),
          axis.ticks.length = unit(.2, "cm"))
  p1<-p1+scale_y_continuous(expand=c(0,0),breaks=seq(range(ETRdat$logp)[1],range(ETRdat$logp)[2],5),
                            labels=seq(range(ETRdat$logp)[1],range(ETRdat$logp)[2],5),
                            limits=c(0,range(ETRdat$logp)[2]+0.5))
  
  filename <- paste("./cov/ETR-",ii,".tiff",sep = "")
  tiff(filename,width=10,height=6,units = "cm",res=300)
  p1
  dev.off()
}


###################################


qp.sig.time <- cov.test(pheno=y2,snp=genqp)
nsig2 <- c()
for (i in 1:11) {
  nsig2 <- c(nsig2,length(which(qp.sig.time$logp[i,]>qp.sig.time$thre1)))
}


for (ii in 1:11) {
  qPdat<-cbind(as.numeric(rownames(genqp)),snpinfoqp$CHR,snpinfoqp$POS,snpinfoqp$REF,snpinfoqp$ALT,qp.sig.time$logp[ii,])
  qPdat<-as.data.frame(qPdat)
  colnames(qPdat)<-c("seq","CHR","POS","REF","ALT","logp")
  qPdat$seq<-as.numeric(qPdat$seq)
  qPdat$POS<-as.numeric(qPdat$POS)
  qPdat$logp<-as.numeric(qPdat$logp)
  qPdat<-qPdat[order(qPdat[,1]),]
  rownames(qPdat)<-NULL
  scaffoldsum <- length(table(qPdat[,2]))
  chrnum<-c(rep(1,1417),rep(2,456),rep(3,451),rep(4,584),rep(5,410),rep(6,531),rep(7,325),rep(8,457),rep(9,231),rep(10,433),rep(11,407),
            rep(12,331),rep(13,368),rep(14,1048),rep(15,391),rep(16,314),rep(17,352),rep(18,381),rep(19,717),rep(20,31),rep(21,7),rep(22,7),
            rep(23,1),rep(24,4),rep(25,25),rep(26,2),rep(27,4),rep(29,6),rep(30,1),rep(31,8),rep(32,10),rep(33,7),rep(34,4),rep(35,5),rep(36,10),
            rep(37,13),rep(38,10),rep(39,4),rep(41,5),rep(42,44),rep(45,9),rep(47,5),rep(48,2),rep(49,14),rep(50,3),rep(51,2),rep(52,30),rep(53,1),
            rep(56,2),rep(57,2),rep(61,1),rep(63,6),rep(64,6),rep(65,2),rep(66,2),rep(67,1),rep(69,4),rep(71,1),rep(74,1),rep(75,7),rep(76,3),
            rep(77,4),rep(79,2),rep(80,3),rep(82,31),rep(85,1),rep(86,1),rep(87,3),rep(90,1),rep(92,2),rep(93,1),rep(95,6),rep(99,3),
            rep(102,2),rep(104,1),rep(107,2),rep(111,7),rep(112,1),rep(119,5),rep(123,8))
  qPdat$CHR<-chrnum
  scaffoldc <- as.numeric(names(table(qPdat[,2])))
  for(i in 1:scaffoldsum){
    qPdat[which(qPdat[,2]==scaffoldc[i]),2] <- i
  }
  for (i in 1:scaffoldsum){ 
    ndx <- which(qPdat[, 2]==i)
    lstMrk <-max(qPdat[ndx, 3])
    if (i < scaffoldsum) ndx2 <-which(qPdat[, 2]==i+1)
    if (i < scaffoldsum) qPdat[ndx2, 3]<-qPdat[ndx2, 3]+lstMrk
  }
  
  bpMidVec <- vector(length=scaffoldsum)
  
  for (i in 1:scaffoldsum){
    ndx <- which(qPdat[, 2]==i)
    posSub <- qPdat[ndx, 3]
    bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
  }
  
  library(ggplot2)
  p4<- ggplot()
  p4<-p4+ geom_point(data=qPdat,aes(x=POS, y=logp,colour=factor(as.numeric(CHR),levels=as.character(seq(1,80)))),size=2,alpha=0.8)
  p4<-p4+scale_color_manual(values=rep(c('#46732EFF', '#D5E4A2FF'), 40))
  # p4<-p4+geom_hline(yintercept=16,linetype=1, col='black', lwd=1)
  p4<-p4+geom_hline(yintercept=qp.sig.time$thre1,linetype=1, col='red', lwd=1)
  p4<-p4+scale_x_continuous(breaks=c(bpMidVec[1:19]),labels=c(1:19))
  p4<-p4+xlab('') + ylab('')+theme_bw()+
    theme(legend.position = "none",
          axis.text.x=element_text(hjust=0.5,size=8,color='black'),
          axis.text.y=element_text(vjust=0.2,size=10,color='black'),
          axis.ticks.length = unit(.2, "cm"))
  p4<-p4+scale_y_continuous(expand=c(0,0),breaks=seq(range(qPdat$logp)[1],range(qPdat$logp)[2],5),
                            labels=seq(range(qPdat$logp)[1],range(qPdat$logp)[2],5),
                            limits=c(0,range(qPdat$logp)[2]+0.5))
  
  filename <- paste("./cov/qP-",ii,".tiff",sep = "")
  tiff(filename,width=10,height=6,units = "cm",res=300)
  p4
  dev.off()
}

#######################

qn.sig.time <- cov.test(pheno=y3,snp=genqn)
nsig3 <- c()
for (i in 1:11) {
  nsig3 <- c(nsig3,length(which(qn.sig.time$logp[i,]>qn.sig.time$thre1)))
}



for (ii in 1:11) {
  qNdat<-cbind(as.numeric(rownames(genqn)),snpinfoqn$CHR,snpinfoqn$POS,snpinfoqn$REF,snpinfoqn$ALT,qn.sig.time$logp[ii,])
  qNdat<-as.data.frame(qNdat)
  colnames(qNdat)<-c("seq","CHR","POS","REF","ALT","logp")
  qNdat$seq<-as.numeric(qNdat$seq)
  qNdat$POS<-as.numeric(qNdat$POS)
  qNdat$logp<-as.numeric(qNdat$logp)
  qNdat<-qNdat[order(qNdat[,1]),]
  rownames(qNdat)<-NULL
  scaffoldsum <- length(table(qNdat[,2]))
  chrnum<-c(rep(1,1357),rep(2,938),rep(3,501),rep(4,556),rep(5,592),rep(6,707),rep(7,442),rep(8,524),rep(9,199),rep(10,435),rep(11,464),
            rep(12,337),rep(13,375),rep(14,686),rep(15,334),rep(16,403),rep(17,319),rep(18,285),rep(19,347),rep(20,47),rep(21,4),rep(22,4),
            rep(23,1),rep(24,1),rep(25,18),rep(26,4),rep(27,1),rep(28,3),rep(29,5),rep(30,2),rep(31,6),rep(32,7),rep(33,2),rep(34,8),rep(35,2),
            rep(37,6),rep(38,4),rep(41,1),rep(42,8),rep(45,10),rep(47,5),rep(48,3),rep(49,2),rep(50,5),rep(51,1),rep(52,10),rep(53,7),rep(56,1),
            rep(61,1),rep(66,3),rep(73,1),rep(75,1),rep(76,1),rep(81,1),rep(86,1),rep(87,2),rep(92,1),rep(94,2),rep(107,3),rep(111,2),rep(122,1),rep(123,1))
  qNdat$CHR<-chrnum
  scaffoldc <- as.numeric(names(table(qNdat[,2])))
  for(i in 1:scaffoldsum){
    qNdat[which(qNdat[,2]==scaffoldc[i]),2] <- i
  }
  for (i in 1:scaffoldsum){ 
    ndx <- which(qNdat[, 2]==i)
    lstMrk <-max(qNdat[ndx, 3])
    if (i < scaffoldsum) ndx2 <-which(qNdat[, 2]==i+1)
    if (i < scaffoldsum) qNdat[ndx2, 3]<-qNdat[ndx2, 3]+lstMrk
  }
  
  bpMidVec <- vector(length=scaffoldsum)
  
  for (i in 1:scaffoldsum){
    ndx <- which(qNdat[, 2]==i)
    posSub <- qNdat[ndx, 3]
    bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
  }
  
  library(ggplot2)
  p4<- ggplot()
  p4<-p4+ geom_point(data=qNdat,aes(x=POS, y=logp,colour=factor(as.numeric(CHR),levels=as.character(seq(1,80)))),size=2,alpha=0.8)
  p4<-p4+scale_color_manual(values=rep(c('#F05C3BFF', 'sienna1'), 36))
  # p4<-p4+geom_hline(yintercept=16,linetype=1, col='black', lwd=1)
  p4<-p4+geom_hline(yintercept=qp.sig.time$thre1,linetype=1, col='red', lwd=1)
  p4<-p4+scale_x_continuous(breaks=c(bpMidVec[1:19]),labels=c(1:19))
  p4<-p4+xlab('') + ylab('')+theme_bw()+
    theme(legend.position = "none",
          axis.text.x=element_text(hjust=0.5,size=8,color='black'),
          axis.text.y=element_text(vjust=0.2,size=10,color='black'),
          axis.ticks.length = unit(.2, "cm"))
  p4<-p4+scale_y_continuous(expand=c(0,0),breaks=seq(range(qNdat$logp)[1],range(qNdat$logp)[2],5),
                            labels=seq(range(qNdat$logp)[1],range(qNdat$logp)[2],5),
                            limits=c(0,range(qNdat$logp)[2]+0.5))
  
  filename <- paste("./cov/qN-",ii,".tiff",sep = "")
  tiff(filename,width=10,height=6,units = "cm",res=300)
  p4
  dev.off()
}
