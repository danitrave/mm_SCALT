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

''' The function filter_table(selected_genes,path_ratios_table) is responsible of generating a temporary file called "tmp_ratios.tsv" which contains
    the probability ratios of those genes for which there was a probability value in the list of pre-compiled probabilities (background_probs.tsv)'''

def filter_table(selected_genes,path_ratios_table):
    t = pd.read_csv(path_ratios_table,sep="\t",header=0,index_col=0)  #load the table with the ratios
    genes2keep = list(selected_genes["genes"])
    filtered_t = t.loc[genes2keep,:]                                    #extract genes to be used
    filtered_t.to_csv("tmp_ratios.tsv",sep="\t",index=True,header=True)   #save the new temp file
    return "tmp_ratios.tsv"

''' The function remove_low_prob_genes_and_variable(path_probs,selected_genes) has to collect the genes to remove during the generation of the lists.
    Specifically, genes lowly expressed and those not respecting the entropy thresholds'''

def remove_low_prob_genes_and_variable(path_probs,selected_genes):
    black_list_genes = []
    t = pd.read_csv(path_probs,sep="\t",header=0,index_col=0)  #table with probabilities per gene per cell type
    genes2keep = list(selected_genes["genes"])
    filtered_t = t.loc[genes2keep,:]               #collect only genes for which there is a pre-computed background probability avaialable
    indexes = list(filtered_t.index)
    for i in indexes:
        L = filtered_t.loc[i,:]
        M_PROB = max(L)
        if M_PROB < 0.10:                 #remove genes with low probability
            if i not in black_list_genes:
                black_list_genes += [i,]
    entropy = pd.read_csv("AnnolistsBuilder_results/genes_entropy.tsv",sep="\t",header=0,index_col=0)   #entropy values calculated per gene
    f_entropy = entropy.loc[genes2keep,:]
    f_entropy.drop(black_list_genes,inplace=True)
    lower_bound = np.quantile(f_entropy["entropies"],0.05)#lower and upper bound entropy thresholds (calculated only with genes expressed in at least 10%
    upper_bound = np.quantile(f_entropy["entropies"],0.95)
    genes_names = list(f_entropy.index)
    for x in genes_names:
        curr_entropy = f_entropy.loc[x,"entropies"]
        if curr_entropy < lower_bound:               #remove genes with low entropy
            if x not in black_list_genes:
                black_list_genes+=[x,]
        elif curr_entropy > upper_bound:             #remove gene with high entropy
            if x not in black_list_genes:
                black_list_genes+=[x,]
    return black_list_genes

''' The function boostrap_lists_generator(boos,path,genes_n,A) has to create the bootstrap lists for each bootstrap sample using only probabilities and
    ratios relative to that bootstrap sample'''

def boostrap_lists_generator(boos,path,genes_n,A):
    try:
        os.system("mkdir ./mi")   #the directoty ./mi wil contain the bootstrap lists
    except:
        print("Error: unable to create the directory mi/.")
    for b in boos:
        try:
            os.system("mkdir mi/"+b)   #create current bootstrap lists directory
        except:
            print("Error: unable to create the directory mi/"+b)
        GENES_DF = pd.read_csv("AnnolistsBuilder_results/TABLE_OF_GENES.tsv",sep="\t",header=0,index_col=0) #list with the genes to consider  
        genes2remove = remove_low_prob_genes_and_variable(path+"/"+b+"/probabilities_tables/cell_type_probabilities.tsv",GENES_DF) #genes to be removed
        adj_table = filter_table(GENES_DF,path+"/"+b+"/probabilities_tables/probabilities_ratio_matrix.tsv") 
        try:
            os.system("Rscript --vanilla scripts/mi_ranker.R tmp_ratios.tsv")  #rank current lists based on the ratios relative to the current bootstrap
        except:
            print("Error: unable to run the script mi_ranker.R.")
        table = pd.read_csv("ranking_genes.tsv",sep=" ",header=0,index_col=None)         #ranking table
        ratios = pd.read_csv("tmp_ratios.tsv",sep="\t",header=0,index_col=0)
        try:
            os.system("mkdir ./boostrap_lists")             #directory that will contain the list of genes per each cell type for the current bootstrap
        except:
            print("Error: unable to create the directory boostrap_lists")
        header = [o.replace(" ",".").replace("–","-").replace("/","_").replace("(","_").replace(")","_").upper() for o in list(ratios.columns)]   #correct columns
        ratios.columns = header
        table.columns = header
        cellTypes = list(table.columns)
        for c in cellTypes:
            currCellMEG_indexes = list(table.loc[:,c])
            ORDERED_GENES_ENS_ID = list(GENES_DF.iloc[currCellMEG_indexes,0])   #get the ordered list of genes from the ranking and the corresponding probability ratio
            prob_ratios_cellType = ratios.loc[ORDERED_GENES_ENS_ID,c]
            cell_name = c.replace(".","_")+"_boostrap_list.txt"
            DF = pd.DataFrame.from_dict({"genes":ORDERED_GENES_ENS_ID,"probRatios":prob_ratios_cellType},orient="columns")   #generate the final table for this cell type
            DF.drop(genes2remove,inplace=True)  #remove genes previously mentioned
            selected = DF.head(genes_n)
            selected.to_csv(cell_name,sep="\t",header=False, index=False)
            os.system("mv "+cell_name+" ./boostrap_lists")   #move the table in the proper directory
        try:
            os.system("rm tmp_ratios.tsv")    #remove temp files and organize outpts in the correct directories
        except:
            print("Error: unable to remove file tmp_ratios.tsv")
        try:
            os.system("mv ranking_genes.tsv boostrap_lists/ mi/"+b+"/")
        except:
            print("Error: unable to move ranking_genes.tsv file in the boostrap_lists/ mi/"+b+"/ directory.")
    return "mi"

''' The function probabilitiesEstimator(L,db,I) has to estimate the probabilities of finding a gene across N bootstrap lists of the same cell type'''

def probabilitiesEstimator(L,db,I):   #L is the list of bootstrap samples
    ls = os.listdir(db+"/0_boostrap_sample/boostrap_lists")
    cts = [n.split("_boos")[0] for n in ls]    #get the name of cell types avaialble from the bootstrap
    D = {c:{} for c in cts}
    for b in L:
        current = os.listdir(db+"/"+b+"/boostrap_lists")   #see the lists of the current bootstrap and count the occurence of each gene in each cell type
        for l in current:
            txt = open(db+"/"+b+"/boostrap_lists/"+l)
            name = l.split("_boos")[0]
            for line in txt:
                g = line.strip("\n").split("\t")
                if g[0] not in D[name]:
                    D[name][g[0]]=1
                else:
                    D[name][g[0]]+=1
            txt.close()
    for cell_type in D:
        for gene in D[cell_type]:
            D[cell_type][gene]=D[cell_type][gene]/I   #estimate probability per gene per cell type
    return D

''' The function probabilitiesTables(ls,pathdb) simply loads the probabilities previosly estimated (probabilitiesEstimator(L,db,I)) in a python dictionary
    reporting the cell type in the keys and the dataframe of probabilites as corresponding value'''

def probabilitiesTables(ls,pathdb):
    K = {}
    for j in ls:
        name = j.split(".ts")[0]
        df = pd.read_csv(pathdb+"/"+j,sep="\t",header=0,index_col=0)   #load df of probabilities
        K[name]=df
    return K

''' The funciton mutual_information(L1,L2,tables,TYPE) calculates the mutual information between two lists'''

def mutual_information(L1,L2,tables,TYPE):
    t = tables[TYPE]     #get probabilities relative to the cell type represented by both L1 and L2
    common = []
    not_common = []
    for y in L1:           #split common and not common genes
        if y not in L2:
            not_common += [y,]
        else:
            common += [y,]

    #### MI shared genes between two lists ####
    MI_common = sum(np.array(t.loc[common,"p"])*(np.log10(np.array(t.loc[common,"p"])/0.0001)))
    #### MI not shared genes between two lists ####
    n = np.array([0.0001]*len(not_common))
    MI_notcommon = sum(n*(np.log10(n/0.0001)))
    MI = MI_common+MI_notcommon
    return MI

''' The function listsComparator(boos,path,midb) sets the imput in order to calculculate the mutual information between all combination of bootstrap lists
    referring to the same cell type but coming from different bootstraps excluding self and redundant comparisons'''

def listsComparator(boos,path,midb):
    ls = os.listdir("./"+midb)
    cts = []
    for a in ls:
        cts += [a.split(".ts")[0],] 
    D = {c:{} for c in cts}
    for i in boos:
        boos_CTs = os.listdir(path+"/"+i+"/boostrap_lists")  #collect the bootstrap from the current sample organizing lists based on the cell type 
        for x in boos_CTs:
            cellType_name = x.split("_boos")[0]
            txt = open(path+"/"+i+"/boostrap_lists/"+x,"r")
            genes_list = []
            for line in txt:
                l = line.strip("\n").split("\t")
                genes_list += [l[0],]
            txt.close()
            D[cellType_name][i]=genes_list
    
    #### Calculate the mutual information per each couple of lists excluding repetitions and self comparison ####
    combi = list(itertools.combinations(boos, 2))
    prob_tables = probabilitiesTables(ls,midb)
    MIs_dict = {}
    for cellType in D:
        current = []
        for com in combi:
            mi = mutual_information(D[cellType][com[0]],D[cellType][com[1]],prob_tables,cellType)
            current += [mi,]
        MIs_dict[cellType]=current
    MIs = pd.DataFrame.from_dict(MIs_dict)   #db with all MIs per combaination
    return MIs

if __name__ == "__main__":
    path = sys.argv[1]   #path to the boostrap directory
    genes = int(sys.argv[2])   #number of genes to extract from each cell type specific lists of genes from the boostrap
    iterations = int(sys.argv[3])   #number of bootstraps samples generated
    annot = sys.argv[4]
    boos = os.listdir("./"+path)
    buildLists = boostrap_lists_generator(boos,path,genes,annot)
    dirCheck = os.path.exists("./miProbTables")
    if dirCheck == False:
        try:
            os.system("mkdir ./miProbTables")
        except:
            print("Error: unable to generate the directory miProbTables.")
        p = probabilitiesEstimator(boos,buildLists,iterations)
        for c in p:
            df = pd.DataFrame.from_dict(p[c],orient='index')
            df.columns = ["p"]
            df.to_csv(c+".tsv",sep="\t",header=True,index=True)
            os.system("mv "+c+".tsv miProbTables/")
        comparison = listsComparator(boos,buildLists,"miProbTables")
        comparison.to_csv("boos_mutual_information.tsv",sep="\t",header=True,index=False)
    else:
        print("Error: the directory miProbTables already exists.")
    
end_time = datetime.now()
print('Generation of bootstrap lists and MIs values for all combinations (binomial coeffient) of lists of the same cell type: {}'.format(end_time - start_time))
