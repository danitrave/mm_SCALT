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
import math
warnings.filterwarnings("ignore")

start_time = datetime.now()

''' The function likelihood_ratio_test(counts) is responsible of performing the likelihood-ratio based test. Each cell is tested indipendently from
    the others. Each cell is tested on over 471 cell-type specific list of genes. Upon likelihood calculation, the AIC is calculated. This compares
    the likelihood for a cell of being of that cell type with respect to the likelihood of being the "general cell type". The higher is the difference
    the lower will be the p-value so being more significant the classification.'''

def likelihood_ratio_test(counts):
    global threshold
    global cell_types_statistics     #global variables to be provided as input from the user
    global cell_types_names
    global notation
    results = []                  #final p-values table and deltas
    differences = []
    #### Counts loading and filtering by threshold ###
    genesExpressed_perCell = list(np.count_nonzero(counts, axis=0))
    trueFalse = ["PASS" if x >= threshold else "EXCLUDE" for x in genesExpressed_perCell]  #remove cells having less than "Threshold" genes expressed 
    counts[counts>0]=1                           #trasform all counts to a 0-1 table
    #### Evaluate each cell ####
    cells = list(counts.columns)
    for i in cells:
        sums = []
        current_diff = []              #evaluate each cell type list indipendently from the others on the same cell
        for ct in cell_types_names:
            grep_counts = counts.loc[cell_types_statistics[ct][notation],i]      #collect the counts of the genes specific for the current cell type on the current cell
            genes = list(grep_counts.index)
            zipped = list(zip(grep_counts,genes))
            expressed_genes = []
            not_expressed_genes = []
            for e in zipped:            #separate the genes specific for the current cell type and expressed by the cell from those specific for the cell type but not expressed by the cell
                if e[0]==1:
                    expressed_genes += [e[1],]
                else:
                    not_expressed_genes += [e[1],]
            EXPRESSED = cell_types_statistics[ct].loc[expressed_genes,"cell_type_prob"]      #collect probabilites P(G|CT)
            GLOBAL_EXPRESSED = cell_types_statistics[ct].loc[expressed_genes,"general_prob"]   	#collect probabilites P(G)
            NOT_EXPRESSED = cell_types_statistics[ct].loc[not_expressed_genes,"minus_cell_type_prob"]    #1-P(G|CT)
            GLOBAL_NOTEXPRESSED = cell_types_statistics[ct].loc[not_expressed_genes,"minus_general"]       #P(P)
            score_ct = sum(np.log10(EXPRESSED).tolist()) + sum(np.log10(NOT_EXPRESSED).tolist())           #Likelihood of being of that cell type
            mean_type = sum(np.log10(GLOBAL_EXPRESSED).tolist()) + sum(np.log10(GLOBAL_NOTEXPRESSED).tolist())   #Likelihood of being the mean type
            diff_scores = score_ct - mean_type     #likelihood comparisons
            if diff_scores <= 0.0:                 #if the likelihood ratio is negative, there is no difference between the two models
                score = 1.0                        #so the p-value is 1 and the delta is 0
                d = 0.0
            else:
                score = math.exp((mean_type-(score_ct))/2)    #else, calculate the AIC
                d = mean_type-(score_ct)
            sums+=[score,]
            current_diff+=[d,]

        results += [sums,]
        differences += [current_diff,]
    return results,trueFalse,differences
    
''' The function parallelizer(counts2paraller,cpus) is responsible of parallizing the computation if multiple processors are give in input.'''

def parallelizer(counts2paraller,cpus):
    if cpus == 1:
        df = pd.read_csv(counts2paraller,sep="\t",header=0,index_col=0)   #one processor
        return likelihood_ratio_test(df)
    else:
        df = pd.read_csv(counts2paraller,sep="\t",header=0,index_col=0)   #multiple processors
        n = df.shape[1]
        s = int(n/cpus)
        k = 0
        slicing = []
        while k < n:
            if k+s >= n:
                slicing += [[k,n],]     #dividie counts table slicing columns in "sub-tables"
                k+=s
            else:
                slicing += [[k,k+s],]
                k+=s
        slices = []
        for SLICE in slicing:
            slices += [df.iloc[:,SLICE[0]:SLICE[1]],]   #collect counts referred to the sub-columns
        final_outcome = []
        final_trueFalse = []
        final_deltas = []
        with concurrent.futures.ProcessPoolExecutor() as executor:    #run multiple processors
            output = executor.map(likelihood_ratio_test,slices)
            for res in output:
                final_outcome.extend(res[0])
                final_trueFalse.extend(res[1])
                final_deltas.extend(res[2])
        return final_outcome,final_trueFalse,final_deltas

''' The function cellTypes(DIR,nota) collectes the cell-type specific lists of genes in a dictionary. The cell types lists are collected from the
	directory "DIR" which can be either "cell_types", "custom" or "naive" depeding on the mode used. Genes are collected depending on the gene 
	notation insert by the user.'''

def cellTypes(DIR,nota):
    ls = os.listdir(DIR+"/")     #collect the lists of genes from the designed directory
    dict_cell_types = {}
    for c in ls:
        name = c.split(".tsv")[0]
        df = pd.read_csv(DIR+"/"+c,sep="\t",header=0,index_col=False)   #open the cell type list
        df.index = list(df[nota])           #collect genes with the proper notation and save the list in a dictionary
        dict_cell_types[name]=df
    return dict_cell_types

if __name__ == "__main__":
    counts = sys.argv[1]            #Sample
    threshold = int(sys.argv[2])    #Minimum number of genes that a cell must have to be validated. The default is 250.
    notation = sys.argv[3]          #Kind of notation to consider
    dirCellTyeps = sys.argv[4]      #name of the directory to use as cell types
    CPUS = int(sys.argv[5])         #Number of CPUs. The default is 1.
    cell_types_statistics = cellTypes(dirCellTyeps,notation)
    cell_types_names = sorted(list(cell_types_statistics.keys()))
    
    lrt_run = parallelizer(counts,CPUS)
    df_lrt = pd.DataFrame(lrt_run[0],columns=cell_types_names)
    filteredCells = pd.DataFrame.from_dict({"geneExpr_Filter":lrt_run[1]})
    filteredCells.to_csv(counts.split(".tsv")[0]+"_genesExpressed_filter.tsv",sep="\t",header=True, index=True)
    df_lrt.dropna(axis=0,inplace=True)
    #### Only for the benchmark ####
    #anno = pd.read_csv("benchmark_labels.tsv",sep="\t",header=0,index_col=None)
    #ref = list(anno["CELL_ANNOTATION"])
    #df_lrt.index = ref
    df_lrt.to_csv(counts.split(".tsv")[0]+"_p_values.tsv",sep="\t",header=True,index=True)
    df_deltas = pd.DataFrame(lrt_run[2],columns=cell_types_names)
    df_deltas.dropna(axis=0,inplace=True)
    #df_deltas.index = ref
    df_deltas.to_csv(counts.split(".tsv")[0]+"_deltas.tsv",sep="\t",header=True,index=True)

end_time = datetime.now()
print('Duration likelihood test: {}'.format(end_time - start_time))
