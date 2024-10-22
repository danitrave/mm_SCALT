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

''' The function entropy_hausekeeping(entropies) is able to calculate an entropy threshold based on the value of entropy of the known housekeeing genes in Homo sapiens.
	NOTE: currently, this function is not used.'''

def entropy_hausekeeping(entropies):
    global notation
    tableHKG = pd.read_csv("config/Housekeeping_GenesHuman.csv",sep=";",header=0,index_col=False)
    hkglist = list(tableHKG["Gene.name"])
    genes_config = pd.read_csv("config/GRCh38.109_ensembID_geneSymbol.tsv",sep="\t",header=0,index_col=0)
    indexes = list(genes_config["gene_symbol"])
    genes_config.index = indexes
    d = {}
    symb = list(genes_config["gene_symbol"])
    ens = list(genes_config["ensembl_id"])
    genes2consider = []
    for i in range(len(symb)):
        d[symb[i]]=ens[i]
    if notation == "gene_symbol":
        for x in hkglist:
            if x in d:
                genes2consider += [x,]
            else:
                continue
    elif notation == "ensembl_id":
        for x in hkglist:
            if x in d:
                genes2consider += [d[x],]
            else:
                continue
    entropies_values = []
    for y in genes2consider:
        if y in entropies:
            entropies_values += [entropies[y],]
    quantile = np.quantile(entropies_values,0.05)
    return quantile

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
        ###if L.shape[0]>1:
            ###black_list_genes += [i,]
            ###d[i]=0
            ###continue
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
    #quantile_threshold = entropy_hausekeeping(d)
    for g in d:
        if g in houkeeping:             #if the gene is among the housekeeping, exclude it
            if g not in black_list_genes:
                black_list_genes += [g,]
            else:
                continue

    genes = []
    entropies = []
    for k in d:
        genes += [k,]
        entropies += [d[k],]

    df = pd.DataFrame.from_dict({"entropies":entropies})    #entropy table
    df.index = genes
    df.to_csv("genes_entropy.tsv",sep="\t",header=True, index=True)
    new = pd.DataFrame.from_dict({"genes":black_list_genes})
    new.to_csv("genes2remove.tsv",sep="\t",header=True, index=True)   #genes to exclude table
    return black_list_genes

if __name__ == "__main__":
    pt = sys.argv[1]       #table with the probabilities per cell type
    notation = sys.argv[2]  #notation required to calculate the minimum threshold for the entropy
    genes2remove = entropy(pt,notation)

