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
import statistics
from operator import itemgetter
import itertools
import seaborn as sns
warnings.filterwarnings("ignore")

start_time = datetime.now()

''' The function probability_calculation(L) calculates the mean of the probabilite of a gene to be expressed in a cell type starting from the collection
	of probabilites tables calculated from the boostrap pipeline.'''

def probability_calculation(L,A):
    if A == None:
        m = sum(L)/len(L)                #mean
        m = m.replace(1.0, 0.9999)          #adjust if the probability is 1 because it is an estimate of the probability
        m.to_csv("genesCellTypes_probabilities.tsv",sep="\t")   #save the probability table
        return m
    else:
        m = sum(L)/len(L)                #mean
        m = m.replace(1.0, 0.9999)          #adjust if the probability is 1 because it is an estimate of the probability
        genesAvailable = list(set(list(m.index)) & set(A))
        filtered_m = m.loc[genesAvailable,:]
        return filtered_m,genesAvailable
    
''' The function general_probabilities(R) calculates the mean of the probabilite of a gene to be expressed in general starting from the collection
	of probabilites tables calculated from the boostrap pipeline.'''

def general_probabilities(R):
    g = sum(R)/len(R)                 #mean
    g.to_csv("genesGeneral_probabilities.tsv",sep="\t")   #save the probability table
    return g
    
''' The function probabilities_ratio(ctp,gp) calculates the mean of the ratio of proabilities previously described starting from the collection
	of probabilites tables calculated from the boostrap pipeline.'''

def probabilities_ratio(ctp,gp):
    cells = list(ctp.columns)
    for ct in cells:
        DIV = ctp.loc[:,ct].div(gp["probs"],axis=0,fill_value=float(0))  #P(G|CT)/G(G)
        ctp[ct]=DIV
    ctp = ctp.fillna(float(0))       #change NA to 0
    ctp = ctp.replace([np.inf, -np.inf], 0)
    ctp.to_csv("genesProbabilities_ratios.tsv",sep="\t",index=True,header=True)   #save the probabilites ratio table

if __name__ == "__main__":
    N = sys.argv[1]    #ensembl_id or gene_symbol
    samples = os.listdir("./boostraps_samples")
    dfs = []
    general = []
    for i in samples:     #collect the tables from each boostrap sample
        p = pd.read_csv("./boostraps_samples/"+i+"/probabilities_tables/cell_type_probabilities.tsv",sep="\t",header=0,index_col=0)  #cell type probability table
        dfs += [p,]
        g = pd.read_csv("./boostraps_samples/"+i+"/probabilities_tables/global_probabilities.tsv",sep="\t",header=0,index_col=0)	#general probability table
        general += [g,]
    ###cell_type_probabilites=probability_calculation(dfs)
    backgroundFileCheck = os.path.exists("config/background_probs.tsv")
    if backgroundFileCheck == False:
        general_probabilities=general_probabilities(general)
        cell_type_probabilites=probability_calculation(dfs,None)
        probabilities_ratio(cell_type_probabilites,general_probabilities)

        #### Extract the total list of genes needed for further analysis ####
        genes = list(general_probabilities.index)
        dfg = pd.DataFrame.from_dict({"genes":genes})
        dfg.to_csv("TABLE_OF_GENES.tsv",sep="\t")
    else:
        read_general = pd.read_csv("config/background_probs.tsv",sep="\t",header=0,index_col=0)
        general_probabilities = pd.DataFrame.from_dict({"probs":list(read_general["probs"])},dtype=float)
        general_probabilities.index = list(read_general[N])
        check_duplicates = {}
        count_only_one = []
        for g in general_probabilities.index:
            if g not in check_duplicates:
                check_duplicates[g]=1
            else:
                check_duplicates[g]+=1
                count_only_one += [g,]
        keep = []
        genees = list(general_probabilities.index)
        already_considered = []
        for e in range(len(genees)):
            number = check_duplicates[genees[e]]
            if number == 1:
                keep += [e,]
            else:
                if genees[e] not in already_considered:
                    keep += [e,]
                    already_considered += [genees[e],]
                else:
                    continue
        filtered_general_probabilities = general_probabilities.iloc[keep,]
        cell_type_probabilites=probability_calculation(dfs,list(filtered_general_probabilities.index))
        filtered2_general_probabilities = filtered_general_probabilities.loc[cell_type_probabilites[1]]
        filtered2_general_probabilities.to_csv("genesGeneral_probabilities.tsv",sep="\t",index=True,header=True)
        cell_type_probabilites[0].to_csv("genesCellTypes_probabilities.tsv",sep="\t")   #save the probability table
        probabilities_ratio(cell_type_probabilites[0],filtered2_general_probabilities)

        #### Extract the total list of genes needed for further analysis ####
        genes = list(filtered2_general_probabilities.index)
        dfg = pd.DataFrame.from_dict({"genes":genes})
        dfg.to_csv("TABLE_OF_GENES.tsv",sep="\t")

end_time = datetime.now()
print('Duration probaility calculation: {}'.format(end_time - start_time))
