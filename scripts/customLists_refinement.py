#!/usr/bin/python

""" Libraries required  """

import sys
import os
import pandas as pd
import numpy as np
from numpy.linalg import norm
from datetime import datetime
import time
import concurrent.futures
import json
import warnings
import statistics
from operator import itemgetter
import itertools
import argparse

warnings.filterwarnings("ignore")

start_time = datetime.now()

''' The function finalLists_generator(CUSTOM,notation) is responsible of refining the cell-type specific lists of genes generated either using the 
	annotation pipeline or the naive one. The program "CTs_lists_generator.py" generates lists having two columns i.e. the genes and the ratio of
	probailities. However the lists require the single probabilites values to be used in the likelihood ratio test. This is the role of the function.'''

def finalLists_generator(CUSTOM,notation,DIR):
    PG = pd.read_csv("genesGeneral_probabilities.tsv",sep="\t",header=0,index_col=0)    #probability of a gene to be expressed im general
    PCTG = pd.read_csv("genesCellTypes_probabilities.tsv",sep="\t",header=0,index_col=0)   #probability of a gene to be expressed in  each cell type
    for c in CUSTOM:
        ct = c.split("_mostExpressedGenes.txt")[0].split(DIR+"/")[1]     #old list of genes
        df = pd.read_csv(c,sep="\t",header=None,index_col=None)
        df.columns = [notation,"ratio"]
        genes = list(df[notation])
        cell_type_prob = np.array(PCTG.loc[genes,ct])    #probability of specific genes in the current cell type
        general_prob = np.array(PG.loc[genes,"probs"])         #probability of specific genes in general
        minus_cell_type_prob = 1 - cell_type_prob  #probability of specific genes not to be expressed in the current cell type
        minus_general = 1 - np.array(general_prob)           #probilities of specific genes not to be expressed in general
        new_df = pd.DataFrame.from_dict({notation:genes,"prob_ratio":list(df["ratio"]),"cell_type_prob":list(cell_type_prob),"general_prob":list(general_prob),"minus_cell_type_prob":list(minus_cell_type_prob),"minus_general":list(minus_general)})
        name = ct.replace("_",".").lower().capitalize()
        try:
            os.system("rm "+c)   #remove the old list of genes
        except:
            pass
        new_df.to_csv(DIR+"/"+name+".tsv",sep="\t",header=True, index=False)    #save the new list of genes in the proper directory

if __name__ == "__main__":
    nota = sys.argv[1]       #notation. Either "ensembl_id" or "gene_symbol"
    mode = sys.argv[2]       #mode. Either "anno" or "naive"
    if mode == "anno":
        custom = os.listdir("./custom")
        custom = ["custom/"+i for i in custom]
        finalLists_generator(custom,nota,"custom")
    elif mode == "naive":
        custom = os.listdir("./naive")
        custom = ["naive/"+i for i in custom]
        finalLists_generator(custom,nota,"naive")
    else:
        print("Error: wrong made inserted")

end_time = datetime.now()
print('Duration lists refinement: {}'.format(end_time - start_time))
