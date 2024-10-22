#!/usr/bin/env Rscript

# The program RDS_RData_formats_converter.R is necessary to convert an ".RData" or ".RDS" file to a ".tsv" format. 

#### Libaris required ####
options(warn=-1)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
args = commandArgs(trailingOnly=TRUE)

start.time <- Sys.time()

#### The function RDS2TSV converts the files from ".RDS"  format to ".tsv" format ####
 
RDS2TSV <- function(rds){
	RDS <- readRDS(file = rds)       			#read the rds
  	table <- data.frame(RDS@assays$RNA@counts)		#convert the rds to a tsv
  	name <- paste(strsplit(rds,".rds")[[1]],".tsv",sep = "")   
  	#file.remove(rds)
  	return(list(table,name,"RDS to TSV conversion: DONE!"))
}

#### The function RData2TSV converts the files from ".RData"  format to ".tsv" format ####

RData2TSV <- function(rdata){
	name <- paste(strsplit(rdata,".RData")[[1]],".tsv",sep = "")    #read the RData
	load(rdata)
  	table <- data.frame(sm)           #convert the RData to a tsv
  	#file.remove(rdata)
	return(list(table,name,"RData to TSV conversion: DONE!"))
}

FILE <- args[1]      #input file
extensionRDS <- grepl("rds",FILE,fixed = TRUE)
extensionRData <- grepl("RData",FILE,fixed = TRUE)   #see the extension of the file
if (extensionRDS==TRUE && extensionRData==FALSE){
	conversion <- RDS2TSV(FILE)
  	write.table(conversion[[1]],conversion[[2]],sep = "\t",quote=F,col.names=T,row.names=T)  #write the comverted tsv
        cat(conversion[[3]])
}else if (extensionRDS==FALSE && extensionRData==TRUE){
	conversion <- RData2TSV(FILE)
	write.table(conversion[[1]],conversion[[2]],sep = "\t",quote=F,col.names=T,row.names=T)  #write the comverted tsv
  	cat(conversion[[3]])
}else{
	cat("File presenting incorrect extension")
}

end.time <- Sys.time()
print(warnings())

time.taken <- difftime(end.time,start.time,units="hours")
cat(paste("The total time required to convert the file in TSV is:",as.character(time.taken),sep=" "))

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
}


