setwd("/home/path/")

Exp<-read.csv("Gene_Expression.csv", header = T, row.names = 1, sep = ",") # column are samples and row are genes
View(Exp)
sampletype<-read.table("Sample_Annoation.txt", sep = ",")[2] # Sample annoation files
View(sampletype)
sampletype <- as.vector(t(sampletype))
names(sampletype) <- colnames(Exp)
Normal_s <- names(sampletype[which(sampletype=="control")])
CRC_s <- names(sampletype[which(sampletype=="Tumor")])

expNormal <- Exp[,Normal_s]
rownames(expNormal) <- rownames(Exp)
expCRC <- Exp[,CRC_s]
rownames(expCRC) <- rownames(Exp)

IFnetwork<-read.csv("GeneRegulatoryNetwork.csv", header = T)
head(IFnetwork)

ntwgenes <- c(t(IFnetwork[,1]),t(IFnetwork[,2]))
ntwgenes <- ntwgenes[!(duplicated(ntwgenes))]
expgenes <- intersect(rownames(Exp),ntwgenes)
notexpgs <- setdiff(ntwgenes,expgenes)


dr <- NULL
for(i in 1:nrow(IFnetwork)){
  ifrom <- intersect(IFnetwork[i,1],notexpgs)
  ito <- intersect(IFnetwork[i,2],notexpgs)
  if(length(ifrom)!=0||length(ito)!=0){
    dr <- c(dr,i)
  }
}
In <- IFnetwork[-dr,]

caledge <- function(exp,network){
  edgefrome <- exp[as.vector(t(network[,1])),]
  edgetoe <- exp[as.vector(t(network[,2])),]
  edgev <- edgefrome-edgetoe
  rownames(edgev) <- NULL
  return(edgev)
}


eNormal <- caledge(expNormal,In)
eCRC<- caledge(expCRC,In)

library(limma) # to identify the dysregulated gene pair
depc <- function(evn,evc,In){
  esetm <- cbind(evn,evc)
  strain <- c(rep(1,ncol(evn)),rep(2,ncol(evc)))
  design <- model.matrix(~-1+factor(strain))
  colnames(design) <- c("control","Tumor")
  fit <- lmFit(esetm, design)
  contrast.matrix <- makeContrasts(Tumor-control, levels=design) 
  fit <- eBayes(fit)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2) 
  options(digits=2)
  top<-topTable(fit2, coef=1, n=1000*1000, adjust="BH")
  avgev <- apply(evc,1,mean)
  In <- cbind(In,avgev)
  In <- In[rownames(top),]
  evc <- evc[rownames(top),]
  top <- cbind(top,In,evc)
  return(top)
}

topDEGP <- depc(eNormal,eCRC,In)

topDEGP<-na.omit(topDEGP)
DEGP<-topDEGP
DEGP_regulate<-data.frame(Gene1=DEGP$Gene1,Gene2=DEGP$Gene2,avgev=DEGP$avgev,logFC=DEGP$logFC, PValu=DEGP$P.Val,DEGP$adj.P.Val)

#write.csv(DEGP_regulate, file="DysregNet_File.txt")

DEGP_regulate<-subset(DEGP_regulate, DEGP_regulate$PValu <0.05)

DEGP_gene<-union(DEGP_regulate$Gene1,DEGP_regulate$Gene2)

regulateDS<-function(DEGP_regulate,d){
  mean_c<-c()
  for(j in 1:length(d)){
    w<-which(DEGP_regulate$Gene1==d[j]|DEGP_regulate$Gene2==d[j])
    if(length(w)!=0){mean_c<-c(mean_c,sum(abs(DEGP_regulate$logFC[w])))}#abs logFC
    else{mean_c<-c(mean_c,0)}
  }
  DS<-data.frame(DS=mean_c[order(-mean_c)])
  rownames(DS)<-d[order(-mean_c)]
  return(DS)
}
coverDS<-function(DEGP_regulate,ds){
  d<-rownames(ds)
  cover_gene<-c()
  cover_gene_length<-c()
  cover_gene_pre<-c()
  for(i in 1:length(d)){
    w<-which(DEGP_regulate$Gene1==d[i])
    if(length(w)!=0){
      cover_gene<-union(cover_gene,DEGP_regulate$Gene2[w])
      cover_gene<-union(cover_gene,d[i])
      cover_gene_length<-c(cover_gene_length,length(cover_gene))
      cover_gene_pre<-c(cover_gene_pre,length(cover_gene)/length(d))
    }else{
      cover_gene<-union(cover_gene,cover_gene)
      cover_gene<-union(cover_gene,d[i])
      cover_gene_length<-c(cover_gene_length,length(cover_gene))
      cover_gene_pre<-c(cover_gene_pre,length(cover_gene)/length(d))
    }
  }
  ds$cover=cover_gene_length
  ds$percentage=cover_gene_pre
  return(ds)
}

DEGP_DS<-coverDS(DEGP_regulate,regulateDS(DEGP_regulate,DEGP_gene))
