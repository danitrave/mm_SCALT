#### PROBABILTY RATIO --> (EXPRESSED IN ONE CELL TYPE)/(EXPRESSED IN ALL CELLS) ####

# This program simply takes the probabilities ratios and order them per each cell type providing the index of the gene.
# Logically, genes are sorted from the one having the highest ratio (more specific) to the one having the lowest ratio (less specific)
# As output, the program generates a table with the ranking of genes (indexes) per each cell type.

table <- read.table("genesProbabilities_ratios.tsv",sep = "\t",header = T,row.names = 1,quote = "",check.names=FALSE)  #import the table with the ranking data
ranked_Ordered_Index_expression <- data.frame(apply(table, 2, order,decreasing=TRUE))-1
colnames(ranked_Ordered_Index_expression) <- colnames(table)
write.table(ranked_Ordered_Index_expression,"genesRanking.tsv",              #and retrieve the index. Substract -1 because the indexes will be processed in python
            sep = " ",quote=F,col.names=T,row.names=T)     #save the matrix having the ordered indexes of the most expressed genes 
#per each cell type from the most specific to the least specific
