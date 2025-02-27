#!/usr/bin/env Rscript

#### This program has to provide the mutual information thresholds per each cell type based on a collection of mutual information values calculated
#### among lists of the same type deriving from the bootstrap

args = commandArgs(trailingOnly=TRUE)

FILE <- args[1]      #input file
data <- read.table(FILE,sep="\t",header=T,check.names=FALSE)
thresholds <- data.frame(ct =colnames(data) ,quantile = rep(NA,ncol(data)),min = rep(NA,ncol(data)),max=rep(NA,ncol(data)),mean=rep(NA,ncol(data)))
row.names(thresholds) <- colnames(data)

for (c in colnames(data)){
  qt <- quantile(data[,c],0.05,na.rm=T)
  MIN <- min(data[,c],na.rm = T)
  MAX <- max(data[,c],na.rm = T)
  MEAN <-mean(data[,c],na.rm = T)
  SD <- sd(data[,c],na.rm = T)
  VAR <- var(data[,c],na.rm = T)
  thresholds[c,"quantile"]<-qt   #first quantile
  thresholds[c,"min"]<-MIN       #minimum entropy
  thresholds[c,"max"]<-MAX       #maximum entropy
  thresholds[c,"mean"]<-MEAN     #mean entropy
}

write.table(thresholds,"thresholdsMI.tsv",sep = "\t",quote = F,row.names = F)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
}
