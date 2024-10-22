#!/usr/bin/python

""" These are the libraries required  """

import sys
import os
import pandas as pd
import numpy as np 
from numpy.linalg import norm
from datetime import datetime
import concurrent.futures
import json
import warnings
#import statistics
import itertools
import seaborn as sns

warnings.filterwarnings("ignore")

start_time = datetime.now()

""" The function rankedMostExpressed(t,g,pr,HKG) has the role of tranforming the ranking of genes per each cell type in actual list of genes
	which will be ordered from the more specific for that cell type to the less specific. Only genes having a probability ratio greater than 1.0 will be considerd."""

def rankedMostExpressed(t,g,pr,HKG):
    GENES_DF = pd.read_csv(g,sep="\t",header=0,index_col=0)      #gene list
    table = pd.read_csv(t,sep=" ",header=0,index_col=None)         #ranking table 
    ratios = pd.read_csv(pr,sep="\t",header=0,index_col=0)          #probabilities ratios 
    cellTypes = list(table.columns)
    result = {}
    header = [o.replace(" ",".").replace("–","-").replace("/","_").replace("(","_").replace(")","_").upper() for o in list(ratios.columns)]   #correct columns
    ratios.columns = header
    cutoff = 100
    for c in cellTypes:
        currCellMEG_indexes = list(table.loc[:,c])
        ORDERED_GENES_ENS_ID = list(GENES_DF.iloc[currCellMEG_indexes,0])   #get the ordered list of genes from the ranking and the corresponding probability ratio
        prob_ratios_cellType = ratios.loc[ORDERED_GENES_ENS_ID,c]
        cell_name = c+"_mostExpressedGenes.txt" 
        DF = pd.DataFrame.from_dict({"genes":ORDERED_GENES_ENS_ID,"probRatios":prob_ratios_cellType},orient="columns")   #generate the final table for this cell type
        DF.drop(HKG,inplace=True)              #eliminate the genes validated by the program "entropy_calculator.py"
        selected = DF[DF["probRatios"]>=2.0]    ####VERIFICA CON LISTE CON RATIO >= 2
        if selected.shape[0] <= cutoff:
            cutoff = selected.shape[0]
        result[c+"_mostExpressedGenes.txt"]=selected
    resultsAfterCutoff = {}
    for r in result:
        resultsAfterCutoff[r]=result[r].head(cutoff)
    return resultsAfterCutoff

''' The function cell_type_list_creator(df,N,p) properly generates the cell type specific list of genes derived either from the naive or annotation table.'''

def cell_type_list_creator(df,N,p):
    df4heatmap = {}
    if p == "naive":
        name_arrival_dir = "./naive"     #collect lists of genes when the "naive" procedure is used
    else:
        name_arrival_dir = "./custom"	#collect lists of genes when the "anno" procedure is used
    for cell_type in df:
        selected = df[cell_type].head(N)
        df4heatmap[cell_type.split("_most")[0].replace("_",".")]=list(selected.iloc[:,0])   #dictionary required for the creation of the heatmap reporting the overlaps between each pair or list of genes
        selected.to_csv(cell_type,sep="\t",header=False, index=False)
        try:
            os.system("mv "+cell_type+" "+name_arrival_dir)   #move the table in the proper directory
        except:
            pass
    return heatmap_table(df4heatmap,p)    #plot the heatmap

''' The function heatmap(df,p) is resposible of creating an heatmap reporting the percentage of overlap between each couple of list of genes'''

def heatmap(df,p):
    if p == "naive":                #naive heatmap
        sns.set(font_scale=0.1)
        HMAP = sns.clustermap(df,cmap="Greens",annot=True,annot_kws={"size": 0.5},fmt=".1f")   #cmap="mako"
        HMAP.savefig("cellTypes_fromNaiveHeatmap.png",dpi=500)
    else:
        sns.set(font_scale=0.1)		#annotation heatmap
        HMAP = sns.clustermap(df,cmap="Greens",annot=True,annot_kws={"size": 0.5},fmt=".1f")   #cmap="mako"
        HMAP.savefig("cellTypes_fromAnnotationHeatmap.png",dpi=500)

''' The function heatmap_table(d,p) has the role of creating the dataframe containing the interesection between all possible combination of cell
	type specific lists of genes. The output is a dataframe required for the creation of the heatmap plot.'''

def heatmap_table(d,p):
    combi = list(itertools.combinations_with_replacement(list(d.keys()),2))   #all possible combination of cell types
    D = {}
    for c in combi:
        MATCH = len(list(set(d[c[0]]) & set(d[c[1]])))    #intersection between the current pair of lists of genes
        D[c]=MATCH
    DF = pd.DataFrame(np.zeros((len(d.keys()),len(d.keys()))),index=d.keys(),columns=d.keys())    #fill the dataframe with the intersection values
    for k in D:
        DF.loc[k[0],k[1]]=D[k]
        DF.loc[k[1],k[0]]=D[k]
    return heatmap(DF,p)

if __name__=="__main__":
    T = sys.argv[1]   #index ranked table
    G = sys.argv[2]   #genes
    PR = sys.argv[3]    #probabilities ratios
    hkg = sys.argv[4]    #genes considered houskeeping and having a max probability among all cell types smaller then 10%
    n = int(sys.argv[5]) #number of genes that will be present in each final list
    path = sys.argv[6]   #either "naive" or "anno" depending on the origin of the annotation
    df_hkg = pd.read_csv(hkg,sep="\t",index_col=0,header=0)
    hgk_genes = list(df_hkg["genes"])
    ranking = rankedMostExpressed(T,G,PR,hgk_genes)
    cell_type_list_creator(ranking,n,path)

end_time = datetime.now()
print('Duration retrieval of the most specific genes ranked by probabilities ratio for a cell type: {}'.format(end_time - start_time))
