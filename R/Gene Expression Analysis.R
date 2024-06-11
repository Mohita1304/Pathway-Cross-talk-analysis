setwd("/home/folder Name")  # set the path where all the Raw data (.cel) are present

#  Normalization of RAW microarray Data 
library(affy)

# Read Raw data (.CEL file) 
mydata<-ReadAffy() 
eset <- rma(mydata)  # RMA normalization
## Identification of Differentially expressed genes using limma 

#matrix construction
Sample<-factor(rep(c("control","Tumor","control","Tumor","control","Tumor","control","Tumor")))
Sample
design.mat<-model.matrix(~0+Sample)   #design matrix 
colnames(design.mat)<-levels(Sample)
design.mat
contrast.mat<-makeContrasts(Diff=Tumor-Control,levels=design.mat) #create Contrast matrix
contrast.mat
# identification of DEGs using limma 
library(limma)
fit<-lmFit(eset,design.mat)
fit2<-contrasts.fit(fit,contrast.mat)
fit3<-eBayes(fit2)
DEGs<-topTable(fit3,coef="Diff",p.value=0.05,adjust.method="fdr",lfc=log2(1),number=nrow(eset)) # default are top 10



##Calculation of correlation between the genes

library(Hmisc)

# Read the expression File  (Gene as row and Samples in column)

Data<-read.csv("Normalised_ExpData_Normal.csv", row.names = 1) # read the files Gene in row and samples in columns

# Calculation of Pearson correlation 

CorMatrix_pValue <- rcorr(t(as.matrix(Data)))

# Flatten Correlation Matrix Function and Create data frame to store flattened correlation matrix

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)                   # Extract upper triangle of correlation matrix
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],   # Extract correlations
    p = pmat[ut]          # Extract corresponding p-values
  ) 
}

# Flatten correlation matrix

matrix_file<-flattenCorrMatrix(CorMatrix_pValue$r, CorMatrix_pValue$P)

#Subset significant correlations based on criteria
sig_CorMatrix<- subset(matrix_file,cor>=0.7|cor<=-0.7 & p<=0.05)
head(sig_CorMatrix)

#Rename column names

colnames(sig_CorMatrix)[1]="Gene1"
colnames(sig_CorMatrix)[2]="Gene2"
head(sig_CorMatrix)

# save the correlation in a new file
write.table(sig_CorMatrix, file="Corr_File.txt")

# Gene pairs with similar correlation pattern  in all the three datasets 

cor_file1 <- read.table("Corr_File1.txt", header = TRUE)
cor_file2 <- read.table("Corr_File2.txt", header = TRUE)
merged_cor <- merge(cor_file1, cor_file2, by = c("Gene1", "Gene2")) # merge if gene pairs is present in both the files 
same_sign <- subset(merged_cor, sign(cor.x) == sign(cor.y)) # to check if the corr pattern is same i.e positive or negative correlation

# rename the column name
colnames(same_sign)[3]="Cor_GSE24514"
colnames(same_sign)[4]="P_GSE24514"
colnames(same_sign)[6]="P_GSE37364"
colnames(same_sign)[5]="Cor_GSE37364"
