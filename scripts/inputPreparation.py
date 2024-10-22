#!/usr/bin/python

""" Libraries required  """

import sys
import os
import pandas as pd
import numpy as np
from numpy.linalg import norm
from datetime import datetime
import concurrent.futures
import json
import warnings
from operator import add
from operator import itemgetter
import random
warnings.filterwarnings("ignore")

start_time = datetime.now()

''' The function table_indexing(counts) is used to convert the counts table in an "indexed" format i.e. a dictionary having the gene as keys
	and the corresponding counts for all cells as values.'''

def table_indexing(counts):
    data = open(counts,"r")
    tsv = data.readlines()
    line_zero = tsv[0].strip("\n").split("\t")      #check is the first row and the last one have the same number of elements
    last_line = tsv[-1].strip("\n").split("\t")
    if len(line_zero) < len(last_line):
        header ="\t".join(["genes/cells"]+line_zero)
    elif len(line_zero) == len(last_line):
        header = "\t".join(["genes/cells"]+line_zero[1:])
    else:
        data.close()
        return print("ERROR: the counts table is not properly formatted.")   #not properly formatted file
    data.close()

    d = {}
    for line in tsv[1:]:
        L = line.strip("\n").split("\t")     #index the counts table
        d[L[0]]=L[1:]

    return d,header
    
''' The function countsRefinement(readCounts,nota,mode) is responible of refine the counts table depending on the mode used and using the desired 
	gene notation. If the user is using the "likelihood" mode, the counts will be refiened reporting only counts of genes listed in the file 
	"GRCh38.109_ensembID_geneSymbol.tsv" located in the ./config directory. If the user is using the "anno" mode or "naive" mode, this passage is not required.'''

def countsRefinement(readCounts,nota,mode,kct):
    nameCounts = readCounts.split(".tsv")[0]
    table = table_indexing(readCounts)        #index tranformation of counts
    tsv = table[0]
    header = table[1]
    pruned = open(nameCounts+"_adj.tsv","w")    #new file
    pruned.write(header+"\n")
    if mode == "likelihood":            #case one: you are doing the likelihood test using the lists from DISCO and HPA. Adj counts inserting missing genes
        if kct == "cell_types" or kct == "low_granularity_cell_types":
            configuration = pd.read_csv("config/GRCh38.109_ensembID_geneSymbol.tsv",sep="\t",header=0,index_col=0)  #gene configuration file
            genes2validate = list(configuration[nota])        #genes to preserve in the counts 
            for gene in genes2validate:
                if gene in tsv:
                    line = "\t".join([gene]+tsv[gene])    #if the original counts present the gene, write the counts
                    pruned.write(line+"\n")
                else:
                    zeros = ["0"]*len(header.split("\t")[1:])   #else, write a collection of zeros 
                    line = "\t".join([gene]+zeros)
                    pruned.write(line+"\n")
            pruned.close()
            return readCounts
        elif kct == "custom":       #case two: you are doing the likelihood test using the custom lists. Adj the counts inserting missing genes if you are testing the custom lists using a different database which does not have the same genes as the "trainig" one
            genesFromLists = []
            customLists = os.listdir("./custom")      #collect the genes from the custom lists without repetitions
            for l in customLists:
                df = pd.read_csv("custom/"+l,sep="\t",header=0,index_col=False)
                genes2introduce = list(df[nota])
                for k in genes2introduce:              #check if the gene was already considered to be introduced in the counts in the next step
                    if k not in genesFromLists:
                        genesFromLists += [k,]
            
            for x in genesFromLists:            #add the missing genes to the counts with a row of zeros
                if x not in tsv:
                    tsv[x]=["0"]*len(header.split("\t")[1:])

            for gene in tsv:               #finally generate the new adjusted counts
                line = "\t".join([gene]+tsv[gene])  
                pruned.write(line+"\n")
            pruned.close()
            return readCounts
    elif mode == "anno":          
        ###configuration = pd.read_csv("config/GRCh38.109_ensembID_geneSymbol.tsv",sep="\t",header=0,index_col=0)  #gene configuration file
        ###genes2validate = list(configuration[nota])        #genes to preserve in the counts 
        ###for gene in genes2validate:
            ###if gene in tsv:
                ###line = "\t".join([gene]+tsv[gene])    #if the original counts present the gene, write the counts
                ###pruned.write(line+"\n")
            ###else:
                ###zeros = ["0"]*len(header.split("\t")[1:])   #else, write a collection of zeros 
                ###line = "\t".join([gene]+zeros)
                ###pruned.write(line+"\n")
        ###pruned.close()
        ####return readCounts
        for gene in tsv:
            line = "\t".join([gene]+tsv[gene])  
            pruned.write(line+"\n")
        pruned.close()
        return readCounts
    elif mode == "naive":                      #"naive" mode
        userLists = os.listdir("./user_lists")          #we need to check if the genes in the lists are either present or not in the tsv
        genes2check = []
        for l in userLists:
            df=pd.read_csv("user_lists/"+l,sep="\t",header=0,index_col=0)
            genes = list(df.iloc[:,0])
            for g in genes:
                if g not in genes2check:        #collect all the genes present in the user-defined cell type lists
                    genes2check += [g,]
        for gene in tsv:
            line = "\t".join([gene]+tsv[gene])   #in first place, add all the genes from the counts table
            pruned.write(line+"\n")

        for c in genes2check:           #check if the genes from the user-defiend lists are present. If not, add the gene with zeros
            if c not in tsv:
                zeros = ["0"]*len(header.split("\t")[1:])   #else, write a collection of zeros
                line = "\t".join([c]+zeros)
                pruned.write(line+"\n")
        pruned.close()
        return readCounts

''' The function annotationRefinement(annotation) is used only when the "anno" mode is on. The program simply format the annotation table 
	in the desidered table organization. '''
	
def annotationRefinement(annotation):
    nameAnnotation = annotation.split(".tsv")[0]
    df = pd.read_csv(annotation,sep="\t",header=0,index_col=None)   #read original file
    anno = list(df.iloc[:,0])
    new_d = {"CELL_ANNOTATION":anno} 
    new_df = pd.DataFrame.from_dict(new_d)      #new file
    new_df.columns = ["CELL_ANNOTATION"]
    new_df.to_csv(nameAnnotation+"_adj.tsv",sep="\t",header=True,index=False)   #save the new tsv
    return annotation

if __name__ == "__main__":
    counts_table = sys.argv[1]      #counts table 
    anno_table = sys.argv[2]        #annotation table
    mode = sys.argv[3]              #kind of mode. Either "likelihood", "anno" or "naive"
    notation = sys.argv[4]          #kind of notation of genes. Either "ensembl_id" or "gene_symbol"
    kind_cell_types = sys.argv[5]   #kind of cell types to use. Either "cell_types", "custom" or "naive"
    if mode == "likelihood" or mode == "naive":
        refined_counts = countsRefinement(counts_table,notation,mode,kind_cell_types)    #likelihood or naive mode
        try:
            os.system("zip -r originalTables_zipped.zip "+refined_counts)   #zip original file
        except:
            print("Error: cannot zip the original counts table!")

        try:
            os.system("rm -r "+refined_counts)								#remove useless files 
        except:
            print("Error: cannot eliminate the original read counts matrix after zipping it")
    elif mode == "anno":													#anno mode 
        refined_counts = countsRefinement(counts_table,notation,mode,kind_cell_types)		#counts
        refined_annotation = annotationRefinement(anno_table)				#anno
        try:
            os.system("zip -r originalTables_zipped.zip "+refined_counts+" "+refined_annotation)   #zip original files
        except:
            print("Error: cannot zip the original counts table!")

        try:
            os.system("rm -r "+refined_counts+" "+refined_annotation)		#remove useless files 
        except:
            print("Error: cannot eliminate the original read counts matrix after zipping it")
    else:
        print("Mode "+mode+' not available. Please select one of the following: "likelihood", "anno" or "naive"')

end_time = datetime.now()
print('Duration input preparation: {}'.format(end_time - start_time))
