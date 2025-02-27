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
from operator import itemgetter
import itertools
import seaborn as sns

''' The function entropy_hausekeeping(entropies) is able to calculate an entropy thresholds (upper bound and lower bound) based on the value of entropy'''

def entropy_threshold(d,low_expressed):
    entrops = []
    #### collect entropy values of genes which are expressed in at least 10% of the cells in one cell type ####
    for k in d:
        if k in low_expressed:
            continue
        else:
            entrops+=[d[k],]

    #### calculate upper bound and lower bound of entropy #####
    lower_bound = np.quantile(entrops,0.05)
    upper_bound = np.quantile(entrops,0.95)
    return lower_bound,upper_bound

''' The function entropy(t,n) is used to calculate the entropy of each gene based on the probability of a gene to be expressed in a cell type.
	Moreover, and more importantly, the function collectes the list of housekeeing genes in the config/ directory and validates each gene 
	present in the current counts table. If the genes is among those reporting in the "config/housekeepingGenes.tsv" file, the gene is not considered
	as a potential cell type specific gene. The genes in the "config/housekeepingGenes.tsv" file were empirically retrieved from over 5 milion cells. 
	Additionally, genes having a max probability of being expressed in any cell type lower than 10%, are not considered.'''

def entropy(t,n):
    table = pd.read_csv(t,sep="\t",header=0,index_col=0)
    housekingTable = pd.read_csv("config/housekeepingGenes.tsv",sep="\t",index_col=0,header=0)   #table with the pre-computed housekeeping genes from 5 milion cells 
    if n == "ensembl_id":
        houkeeping = list(housekingTable["ensembl_id"])     #take the genes with the desidered notation
    elif n == "gene_symbol":
        houkeeping = list(housekingTable["gene_symbol"])
    else:
        print("Error: wrong notation provided!")
    indexes = list(table.index)
    d = {}
    black_list_genes = []    #list that will contain the genes not to be considered as potential cell type specific
    c1 = 0
    c2 = 0
    for i in indexes:
        curr_entropy = 0			#final entropy
        L = table.loc[i,:]
        M_PROB = max(L)
        if M_PROB < 0.10:              #genes for which the max probability is less then 10% must be excluded
            black_list_genes += [i,]
        for p in L:                     #entropy calculation
            if p == 0.0:
                e = 0.0
            else:
                e = -1*((p)*(np.log(p)))   #single cell type entropy
            curr_entropy += e
        d[i]=curr_entropy
    boundary_threshold = entropy_threshold(d,black_list_genes)   #estimate thresholds
    lower_bound = boundary_threshold[0]   
    upper_bound = boundary_threshold[1]
    #### Add to the black list genes those that do not respect the thresholds ####
    for y in d:
        if d[y] < lower_bound:
            if y not in black_list_genes:
                black_list_genes+=[y,]
        elif d[y] > upper_bound:
            if y not in black_list_genes:
                black_list_genes+=[y,]
    genes = []
    entropies = []
    for k in d:
        genes += [k,]
        entropies += [d[k],]

    df = pd.DataFrame.from_dict({"entropies":entropies})    #entropy table
    df.index = genes
    df.to_csv("genes_entropy.tsv",sep="\t",header=True, index=True)
    new = pd.DataFrame.from_dict({"genes":black_list_genes})
    new.to_csv("genes2remove.tsv",sep="\t",header=True, index=True)   #genes to exclude table which include
    return black_list_genes

if __name__ == "__main__":
    pt = sys.argv[1]       #table with the probabilities per cell type
    notation = sys.argv[2]
    genes2remove = entropy(pt,notation)

