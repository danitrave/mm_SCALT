#!/usr/bin/env Rscript

# This program is responsible of doing the naive classification of each cell exploiting the hypergeometric distribution and the 
# list of genes of each cell type previously downloaded. 
# The program takes three inputs:
# -the read counts;
# -a number indicating a threshold below which the cell is considered not good since it expressed a small number of genes;
# -the gene notation.
# The program considers each cell indipendently. The genes expressed in each cell are collected and tested with each list of genes 
# in each cell type using an hypoergeometric distribution testing how likely is to find such intersection between the two list of genes
# just by chance. In each test i.e. per each cell type test, a FDR if obtained. The lower, the more evident is that the such event
# would not have happend just by chance suggesting the most probable annotation for the cell. 
# Hence, each cell will have one FDR per each cell type tested and the lowest will suggest the naive annotation for the cell.
# Cells having a number of genes below the user specified threshold are annotated as "Unknown" and their FDR is 1 for all cell 
# types tested. Moreover, if the cell passes the threshold, it will be subjected to the naive classification but it must pass a second
# test to be finally annotated. If the FDR difference between the lowset FDR and the second lowest one is smaller than 10e-3, 
# the cell is somehow 'umbiguous' so it will be annotated as "Unknown" too and th FDR are set to 1 for each list of genes tested.
# Finally, if the cell passes all thresholds, it will be annotated. 
# The program generate two outputs:
# -one table having a number of row equal to the number of cells and two columns. The first is the most probable annotation, 
#  the second is its FDR of the hypergeometric test;
# -another table having a number of rows equal to the number of cells and N columns one per each cell type tested. Each cell of the table
#  is filled with the FDR that the cell gained in the hypergeometric tests in that cell type.

args = commandArgs(trailingOnly=TRUE)

start.time <- Sys.time()

# The function hypergeometric_calculation <- function(k,K,N,n) is resposible to do the hypergeometric test for the cell and cell type
# under analysis. In the hypergeometric test, the parameters indicate the following numbers:
# -q is the number of genes resulting from the intersect between the genes expressed by the cell and the genes specific for a given cell type;
# -m is the number of genes resulting from the intersect between the genes expressed by your cell and the collection of genes 
#  from the UNION (without replacement) of all lists of genes from all cell types tested;
# -n is the number of genes in the union minus m;
# -k is the number of genes specific for that cell type.

hypergeometric_calculation <- function(k,K,N,n){      #function that performs the hypergeomtric test
	p <- phyper(k,K,N,n,lower.tail = F)                  #p-value. The test is left sided
	return(p)
}

# The function list_of_cell_types <- function(V) has the role of generating an R list were each element is the collection of genes
# specific for a cell type downloaded from the Human Protein Cell Atlas and saved in the directory called cells_tables_hpca.
# Each element in the list is flagged wth the name of the cell type to which the list of genes refers to. Each list of genes is a 
# vector of characters.

list_of_cell_types <- function(V,N){
	L <- vector(mode='list', length=length(V))       #generate and empty list having as much elements as the number of cell types considered
	NAMES <- rep("",length(V))                       #vector of characters that will be filled with the name of the cell types
	for (x in 1:length(V)){
		name <- chartr("_"," ",strsplit(V[x], split = ".tsv")[1])    #name of the cell type
		read <- paste("./user_lists/",V[x],sep="")                #colled the ensembl id of the genes of this cell type in a vector
		C <-  read.table(read,sep="\t",header=T,row.names=1)
		L[[x]] <- as.vector(C[,N])       #put the vector of characters in the list respecting the gene notation required
		NAMES[x]<-name
	}
	names(L)<-NAMES              #update the names of each vector of characters with the corresponding cell type name
	return(L)
}

counts <- args[1]                                          #counts table
name_counts = strsplit(counts, split = ".tsv")[1]
table <- read.table(counts,sep="\t",header=T,row.names=1)                #read counts table
threshold <- strtoi(args[2])                                          #number of genes to keep a cell 
notation <- args[3]							#notation. Either ensembl_id or gene_symb
custom_cells <- list.files("./user_lists")
numberCellTypes_analyzed <- length(custom_cells)
CELLS_TYPES <- list_of_cell_types(custom_cells,notation)               #R-list having the N list of genes from HPCA per cell type as elements
UNION <- c()
for (k in names(CELLS_TYPES)){                            #do the union of all genes in CELLS_TYPES
	genes_to_union <- as.vector(unlist(CELLS_TYPES[k]))
	UNION <- union(UNION,genes_to_union)
}

sample_cells <- names(table)           #names of the cells to be annotated 
genes <- row.names(table)                  #complete list of all genes 
df <- data.frame(matrix(NA,nrow=length(sample_cells),ncol=length(CELLS_TYPES)))     #empty dataframe that will contain all p-values from the tests
naive_annotation <- data.frame(CELL_ANNOTATION = rep(NA,length(sample_cells)),FDR = rep(NA,length(sample_cells)))   #empty datafranme that will contain the final naive annotation and the correspoding p-value for each cell
row.names(df) <- sample_cells
row.names(naive_annotation) <- sample_cells   #put the same names of the cells in the output tables 
for (i in sample_cells){                   #consider each cell in the read couts table indipendently
	TF <- as.vector(table[i]>0)            #genes expressed by the current cell
	cell_genes <- genes[TF]
	if (length(cell_genes) <= threshold){           #if the cell expresses less genes than the threshold, annotate is as "unknown"
		df[i,] <- rep(1.0,length(CELLS_TYPES))
		naive_annotation[i,"CELL_ANNOTATION"] <- "unknown"
		next
	}
	minimum_pvalue <- 1              #base line p-value and annotation
	annotation <- "unknown"
	cell_p_values <- c()              #empty vector that will be filled with the FDR of each hypergeometric test
	for (j in names(CELLS_TYPES)){             #test each cell type specific gene list on the current cell
		genes_test <- as.vector(unlist(CELLS_TYPES[j]))    #vector of genes specific for the current cell type J
		common <- intersect(cell_genes,genes_test)        #intersect between the genes expressed by the current cell and the current cell type
		q <- length(common)                         #q value
		m <- length(intersect(cell_genes,UNION))    #m value
		n <- length(UNION)-m                        #n value
		K <- length(genes_test)                     #k value
		p_value <- hypergeometric_calculation(q,m,n,K)*numberCellTypes_analyzed        #FDR from the hypergeometric test
		cell_p_values <- c(cell_p_values,p_value)
		if (p_value < minimum_pvalue){             #survival condition: the lowest FDR will be saved and its corresponding naive annotation
			minimum_pvalue <- p_value
			annotation <- j
		}	
	}

	cell_p_values <- sort(cell_p_values,decreasing=F)           #sort the FDRs of this cell from the lowest to the highest
	cell_p_values <- replace(cell_p_values,cell_p_values>1,1.0)   #replace FDR greater than 1 with 1 since FDR cannot be hugher than 1
	if (cell_p_values[1]==0 || cell_p_values[2]==0){      #CASE ONE: low quality cell --> "unknown"
		df[i,] <- cell_p_values
		naive_annotation[i,"CELL_ANNOTATION"]<-"unknown"
		naive_annotation[i,"FDR"] <- 1.0
	}else if ((-log10(cell_p_values[1]))-(-log10(cell_p_values[2])) >= 3.0){   #CASE TWO: successful cell --> collect annotation
		df[i,] <- cell_p_values
		naive_annotation[i,"CELL_ANNOTATION"]<-annotation
		naive_annotation[i,"FDR"] <- minimum_pvalue
	}else{												  #CASE THREE: unsucessful hypergeometric test --> "unknown"
		df[i,] <- cell_p_values
		naive_annotation[i,"CELL_ANNOTATION"]<-"unknown"
		naive_annotation[i,"FDR"] <- 1.0
	}
}

print(warnings())
annotation2write <- paste(name_counts,"_naive_annotation.tsv",sep="")       
write.table(naive_annotation,annotation2write,sep = "\t",quote=F,col.names=T,row.names=F)   #save the FDR table
pvalue_table2write <- paste(name_counts,"_FDR_table.tsv",sep="")
write.table(df,pvalue_table2write,sep="\t",quote=F,col.names=T,row.names=T)                 #save the annotation table

end.time <- Sys.time()
time.taken <- difftime(end.time,start.time,units="hours")
print(paste("The total time required was:",as.character(time.taken),sep=" "))

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
}

