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
warnings.filterwarnings("ignore")

start_time = datetime.now()

""" The function RATIO_PROBABILITY_MATRIX_CREATION(single,general) is responsible of calculation the ratio of probilities defined as the 
	probaility of a gene to be expressed in one cell type over the probability of a gene to be expressed in general. """
	
def RATIO_PROBABILITY_MATRIX_CREATION(single,general):   #cell type probabilities and general probabilites
    cells = list(single.columns)
    for ct in cells:
        DIV = single.loc[:,ct].div(general["probs"],axis=0,fill_value=float(0))  #P(G|CT)/P(G) 
        single[ct]=DIV
    single = single.fillna(float(0))
    single = single.replace([np.inf, -np.inf], 0)
    single.to_csv("probabilities_ratio_matrix.tsv",sep="\t",index=True,header=True)

""" The function cellTypeSpecificExpressionProbabilities(counts,annotation,genes) is responsible of calcating the probability of a gene
	to be expressed in one cell type. Starting from the table with all counts from the cells of the equilibrate cell type table and
	the corresponding annotation, first the function separates the counts based on the cell type annotation, then it calculates the 
	probability of a gene to be expressed in that cell type based on the number of cells of that cell type. Again, counts are negletted
	in the sense that is considers only if a gene is expressed or not in a cell. """

def cellTypeSpecificExpressionProbabilities(counts,annotation,genes):
    genes_table = pd.read_csv(genes,sep="\t",header=0,index_col=0)    #read the counts table 
    g = list(genes_table.loc[:,"genes"])
    I = [e for e in range(len(list(counts.columns)))]
    counts[counts>1]=1                             #transform counts as 1 if the gene is expressed or 0 if not
    CELLS = list(annotation["CELL_ANNOTATION"])
    IDS = list(annotation.index)
    IDS_CELLS_ZIP = list(zip(IDS,CELLS))           #zip indexes of the cell and their correspinding annotation (example (1,T-cell))
    D = {}
    CELLS_COUNTS = {}    #number of cells per each cell type
    for i in IDS_CELLS_ZIP:
        if i[1] not in CELLS_COUNTS:      #count the number of cells of that cell type
            CELLS_COUNTS[i[1]]=1
        else:
            CELLS_COUNTS[i[1]]+=1
        currCellCount = list(counts.iloc[:,i[0]])
        if i[1] not in D:                            #sum the transformed count per each gene of cells coming from the same cell type
            D[i[1]]=currCellCount
        else:
            update = list(map(add,D[i[1]],currCellCount))
            D[i[1]]=update
    final = {}
    for x in D:
        for y in CELLS_COUNTS:
            if x==y:
                DIVIDED = [z/CELLS_COUNTS[y] for z in D[x]]    #probabilty calculation
                final[y.replace(" ",".").replace("–","-").replace("/","_").replace("(","_").replace(")","_").upper()]=DIVIDED
    TABLE = pd.DataFrame.from_dict(final,dtype=float)
    TABLE.index=g
    C = list(TABLE.columns)
    TABLE.to_csv("cell_type_probabilities.tsv",sep="\t",index=True,header=True)
    return TABLE

""" The function generalExpressionProbabilies(f1) simply reads the table with the cells (equilibrate cell number) and calculates the 
	probibility of a gene to be expressed in this collection of cells. So, the cell type is not considered and the probaility is 
	calculated over all cells. Notice that counts are negletted in the sense that counts are considered as 1 if the gene is expressed, 
	0 if it is not expressed """

def generalExpressionProbabilies(f1):
    FILE = open(f1,"r")                 #open the table with the counts
    table = FILE.readlines()
    numberOfCells = len(table[0].strip("\n").split("\t"))
    output = {}
    for line in table[1:]:
        L = line.strip("\n").split("\t")
        gene = L[0]
        S = sum([1 for i in L[1:] if float(i) != 0.0])/numberOfCells    #probaility of that gene to be expressed in general.
        output[gene]=S
    FILE.close()
    ratioCells = pd.DataFrame.from_dict(output,orient="index",dtype=float,columns=["probs"])
    ratioCells.to_csv("global_probabilities.tsv",sep="\t",index=True,header=True)
    return ratioCells

if __name__ == "__main__":
    t = sys.argv[1]         #table with the selected cells per cell type
    a = sys.argv[2]         #annotation of the table with the selected cells per cell type
    G = sys.argv[3]         #table with the order succession genes
    T = pd.read_csv(t,sep="\t",header=0,index_col=0)
    ANNO = pd.read_csv(a,sep="\t",header=0,index_col=False)
    general = generalExpressionProbabilies(t)                   #expression of a gene in a cell in general
    ###read_general = pd.read_csv("config/background_probs.tsv",sep="\t",header=0,index_col=0)
    ###general = pd.DataFrame.from_dict({"probs":list(read_general["probs"])},dtype=float)
    ###general.index = list(read_general[N])
    ###general.to_csv("global_probabilities.tsv",sep="\t",index=True,header=True)
    singleCellType = cellTypeSpecificExpressionProbabilities(T,ANNO,G)     #expression of a cell in a cell type
    RATIO_PROBABILITY_MATRIX_CREATION(singleCellType,general)           #ratio of probiblities (expressed|cell/expressed in general)

end_time = datetime.now()
print('Duration probabilities tables generation: {}'.format(end_time - start_time))


