# This program simply takes the probabilities ratios and order them per each cell type providing the index of the gene (specifically in the MI test).
# Logically, genes are sorted from the one having the highest ratio (more specific) to the one having the lowest ratio (less specific)
# As output, the program generates a table with the ranking of genes (indexes) per each cell type.

args = commandArgs(trailingOnly=TRUE)

FILE <- args[1]
table <- read.table(FILE,sep = "\t",header = T,row.names = 1,quote = "")  #import the table with the ranking data
ranked_Ordered_Index_expression <- data.frame(apply(table, 2, order,decreasing=TRUE))-1
write.table(ranked_Ordered_Index_expression,"ranking_genes.tsv",              #and retrieve the index. Substract -1 because the indexes will be processed in python
            sep = " ",quote=F,col.names=T,row.names=T)     #save the matrix having the ordered indexes of the most expressed genes 
#per each cell type from the most specific to the least specific

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
}

