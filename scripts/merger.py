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

''' The function probabilitiesDictionary(DIR) is required to collect the probability tables needed for the calculation of the mutual information.
    Basically, using the program MI.py in the corresponding directory, we have extraxcted the probability of a gene to be found in 100 boostra lists
    from the same cathegory. Still using the program MI.py, we estimated the thresholds that quantify the value of MI that two lists must have to be 
    considered of the same kind. In this function, we simply create a dictionary were we collect the probabilities of each gene per each cell type.
    At the end, we totally have 512 {cell tyep: df} occurencies. Each dataframe will contain the estimates for all the genes that have emerged from at
    least one boostrap list of the same cell type.'''

def probabilitiesDictionary(DIR):
    ls = os.listdir(DIR)
    K = {}
    for j in ls:
        if j == "PLATELETS.tsv" or j == "MIXED_CELL_TYPES.tsv" or j == "MIXED_IMMUNE_CELLS.tsv":   #exclude these cell types from the HPA
            continue
        name = j.split(".ts")[0]
        df = pd.read_csv(DIR+"/"+j,sep="\t",header=0,index_col=0)      #create the {cell type:df} occurence for the correspoding cell type
        K[name]=df
    return K

''' The function collectData(t,path) has two roles: 1) collecting the gene-boostrap-probabilities making use of the function probabilitiesDictionary(DIR);
    2) collect the mutual information thresholds (quartile, min, max and mean) per each cell type. Then, it merges these data in a unique dictionary
    presenting the following structure: {cell type:[probabilities dataframe, MI thresholds circunscribed to the cell type]}'''

def collectData(t,path):    #t is the df with the mi thresholds while path is the path towards the gene-boostrap-probabilities tables 
    mergedDATA = {}
    ps = probabilitiesDictionary(path)                   #collect gene-boostrap-probabilities tables
    df = pd.read_csv(t,sep="\t",header=0,index_col=0)       #threshold mi
    newRowNames = [i for i in list(df.index)]            #homogenize name of the cell types to those modified in probabilitiesDictionary(DIR) 
    df.index = newRowNames
    for e in ps:
        mergedDATA[e]=[ps[e],df.loc[e,:]]     #merge data
    return mergedDATA

''' The function collectCellTypeDefiningLists(path,dbn) is used to collect the final cell type defining list emerging from the extensive re-analysis
    of the databases i.e HPA and DISCO. At the end, we get a dictionary of this kind: {cell type name: [list of the 100 defying genes]}'''

def collectCellTypeDefiningLists(path):
    ls = os.listdir(path)      #path of the lists
    d = {}
    for x in ls:
        n = x.replace(" ",".").replace(".","_").replace("–","-").replace("/","_").replace("(","_").replace(")","_").upper().split("_TSV")[0]
        F = open(path+"/"+x)         #set the name properly, add the name of the database to avoid key substitution and open the file
        l = []
        for line in F:
            a = line.strip("\n").split("\t")    #collect the genes defining that cell type
            if a[0] =="gene_symbol" or a[0] =="ensembl_id":
                continue
            l += [a[0],]
        F.close()
        d[n]=l             #create the dictionary item
    return d

''' The function removeGenesNotAvailable(R,Gs) is resposible of filtering the list of genes removing those genes for which there is no probability estimate
    in the current cell tyep. At the level of MI, these genes have zero contribution so they should be removed becuase they do not add any information of 
    of one lists with respect to the other and viceversa'''

def removeGenesNotAvailable(R,Gs):    #R is the lits of genes while Gs are the genes for which there is a probability estimate (gene-boostrap-probabilities)
    RR = []
    for o in R:
        if o not in Gs:     #if the gene is not present among those available, skip it
            continue
        else:
            RR+=[o,]         #add the gene to the filtered one if there is a probability estimate in the current cell type
    return RR

''' The function mutual_information(unfilteredL1,unfilteredL2,tables,TYPE) is required to calculate the mutual information between two lists. Let's 
    clarify the inputs: 1) unfilteredL1 is the list of cell type one ; 2) unfilteredL2 is the list of cell type two; 3) tables contain all the 
    gene-boostrap-probabilities tables; 4) TYPE is the name of the cell type one. In this context, we calculate the MI between the two lists using
    the gene-boostrap-probabilities from cell type one. In fact, the MI(A,B) is different with respect to MI(B,A) since the probabilities used are
    different eventhought the information coincindes. Therefore, we have to calculate the mutual information between all permutations of cell type 
    specific lists of genes'''

def mutual_information(unfilteredL1,unfilteredL2,tables,TYPE):
    t = tables[TYPE][0]                    #select the gene-boostrap-probabilities table of cell type TYPE i.e. the first list
    availableGenesProbs = list(t.index)
    L1 = removeGenesNotAvailable(unfilteredL1,availableGenesProbs)   #remove genes for which there is no probability from the gene-boostrap-probabilities
    L2 = removeGenesNotAvailable(unfilteredL2,availableGenesProbs)
    common = []
    not_common = []
    for y in L1:               #subdivde genes between those shared by the two lists and those not shared by the lists
        if y not in L2:
            not_common += [y,]
        else:
            common += [y,]

    #### MI shared genes between two lists ####
    MI_common = sum(np.array(t.loc[common,"p"])*(np.log10(np.array(t.loc[common,"p"])/0.0001)))
    #### MI not shared genes between two lists ####
    n = np.array([0.0001]*len(not_common))
    MI_notcommon = sum(n*(np.log10(n/0.0001)))  #if a gene is not shared, its contribution to the MI is ZERO
    MI = MI_common+MI_notcommon    #final mutual information
    return MI

''' The function MI_permutation(cts,probs) has to generate all the permutation of cell type names and calculating the relative mutual information 
invoking the function mutual_information(unfilteredL1,unfilteredL2,tables,TYPE). Additionally, the different corresponding thresholds are added'''

def MI_permutation(cts,probs):    #cts is the dictionary previously described containg all the information necessary for MI calculation
    permu = list(itertools.permutations(list(cts.keys()), 2))    #permutuations
    finalDict = {"CT_miProbs":[],"CT_toCompare":[],"mi":[],"quantile":[],"min":[],"max":[],"mean":[]}  #final table
    for i in permu:
        stats = list(probs[i[0]][1])                #thresholds
        mi = mutual_information(cts[i[0]],cts[i[1]],probs,i[0])    #calculate the MI between these two lists
        finalDict["CT_miProbs"]+=[i[0],]
        finalDict["CT_toCompare"]+=[i[1],]
        finalDict["mi"]+=[mi,]
        finalDict["quantile"]+=[stats[0],]
        finalDict["min"]+=[stats[1],]
        finalDict["max"]+=[stats[2],]
        finalDict["mean"]+=[stats[3],]

    df = pd.DataFrame.from_dict(finalDict)                       #final dataframe with all the MIs of all permutuations
    df.to_csv("experiment_MIs.tsv",sep="\t",header=True,index=True)
    return df

''' The function recursiveMerging(s) has the role of groupping the cell types that share high mutual information among each others. In light of that,
    the strategy adopted is the following: if list A has high MI with B and C, B has high mutual information with A, C and D and C has high mutual 
    information with E and F, the "transitive information" is exploited meaning that all the lists A,B,C,D,E and F will be merged. The strategy is 
    recursive in the sense that we start from the list presenting high MI with more lists and we continue untill all possible "transitive information"
    are exhauseted'''

def recursiveMerging(s):     #s is a list of lists containg the lists to merge sorted on the basis of the numebr of lists involved (transitive excluded)
    merged = []
    already_considered = []   #list that will containg the indexes of the lists already merged 
    guard = 0                    #this variable keeps track if other transitve information are still to be evaluated
    for i in range(len(s)):
        if i in already_considered:       #if these collection of list was already considered, do not further add it to a list
            continue
        already_considered += [i,] 
        current = s[i]                 #current lists to add
        s[i] = []
        for j in range(len(s)):              #compare the current group of lists to all the remaining ones 
            if j not in already_considered:
                c = 0
                for x in s[j]:              #evaluate each single cell type
                    if x in current:
                        c += 1
                        current += s[j]     #if at least one list is shared with the currently analyzed, there is "transitive infromation" so merge all
                if c != 0:                #if "transitive infromation" is found, update the guard and remove the lists merged
                    guard += 1
                    already_considered += [j,]
                    s[j] = []
        #NB: there will be no more "transitive information" if guard is 0 
        singleMerged = current 
        merged += [singleMerged,]

    mergedWithoutReps = []
    for x in merged:                #remove eventual cell type repetitions inside each merged list
        remove_current_reps = []
        for y in x:
            if y not in remove_current_reps:
                remove_current_reps += [y,]
        mergedWithoutReps += [remove_current_reps,]
    return mergedWithoutReps,guard

''' The function merge(s,g) is required to call actual merging function recursivelly. As long as "transitive information" are present, keep calling
    the merging function recursivelly'''

def merge(s,g):
    mergedRecursivelly_base = recursiveMerging(s)   #base case: initial merging
    g = mergedRecursivelly_base[1]
    new_s = mergedRecursivelly_base[0]
    while g != 0:
        mergedRecursivelly = recursiveMerging(new_s)   #keep calling the merging function until convergence i.e. no more "transitive information" 
        g = mergedRecursivelly[1]
        new_s = mergedRecursivelly[0]
    return new_s

''' The function filteringFunction(miTable) simply filters all the MIs calculated among all permutation and collect only those that have a MI above 
    the minimum threshold'''

def filteringFunction(miTable):
    groupped = {}
    allSingle = []
    for i in list(miTable.index):    #filter by the minimum
        if miTable.iloc[i,2] >= miTable.iloc[i,4]:  #index 3 is quartile while index 4 is min. 6 is the mean and 7 in the corrected min 
            allSingle += [list(miTable.iloc[i,:]),]
    for e in allSingle:
        if e[0] not in groupped:      #group cell types 
            groupped[e[0]]=[e[1],]
        else:
            groupped[e[0]]+=[e[1],]

    ITEMS = [(i[0],i[1],len(i[1])) for i in list(groupped.items())]
    sorted_items = sorted(ITEMS,key=itemgetter(2),reverse=True) #sort on the basis of the number of cell types present in each group (decreasing order)
    new_sorted_items = [[x[0]]+x[1] for x in sorted_items]
    return merge(new_sorted_items,groupped)    #merge

''' The function output_report(data) generates the report reporting the results from the mutual information experiment'''

def output_report(data,output_name,table):
    groupped = table.groupby(['CT_miProbs',"CT_toCompare"])
    result = {}
    n = 1
    for i in data:
        permu = list(itertools.permutations(i, 2))    #permutuations
        report_ske = {}
        new_cell = "Cell "+str(n)
        n+=1
        for e in permu:
            for y in groupped:
                if y[0][0]==e[0] and y[0][1]==e[1]:
                    report_ske[(e[0],e[1])]=[float(y[1]["mi"]),float(y[1]["min"]),float(y[1]["quantile"]),float(y[1]["mean"]),float(y[1]["max"])]
        result[new_cell]=report_ske
    
    #### Build the report ####
    report = open(output_name,"w")
    report.write("Merge custom cell type specific lists of genes through mutual information."+"\n")
    report.write("Results"+"\n")
    for cell in result:
        report.write("Evaluating MIs values, hypothetically, cell type "+cell+" should include:"+"\n")
        m = []
        metrics = {}
        for t in result[cell]:
            metrics[t] = "\t".join([str(a) for a in result[cell][t]])
            if t[0] not in m:
                m += [t[0],]
            if t[1] not in m:
                m += [t[1],]
        for x in m:
            ct =  x.replace("_",".").lower().capitalize()
            report.write(ct+"\n")

        report.write("Mutual information values and metrics corresponding thresholds are reported below:"+"\n")
        report.write("\t".join(["Type 1","Type 2","MI","Min","Quantile","Mean","Max"])+"\n")
        for q in metrics:
            report.write(q[0].replace("_",".").lower().capitalize()+" "+q[1].replace("_",".").lower().capitalize()+" "+metrics[q]+"\n")
    
    report.close()
if __name__== "__main__":
    t = sys.argv[1]    #file with MIs thresholds
    pathProbs = sys.argv[2]  #path to the directory containing the probabilities from DISCO
    pathLists = sys.argv[3]   #path to the lists to be merged by mi
    out = sys.argv[4]        #name for the output file
    #### Data collection ####
    data = collectData(t,pathProbs)
    
    #### Collect cell type defining lists ####
    cell_types = collectCellTypeDefiningLists(pathLists)
    #### Finally, calculate the mutual information ####
    miTable = MI_permutation(cell_types,data)
    miTable = pd.read_csv("experiment_MIs.tsv",sep="\t",header=0,index_col=0)
    lsts = filteringFunction(miTable)
    report = output_report(lsts,out,miTable)

end_time = datetime.now()
print('Duration lists merging using mutual information: {}'.format(end_time - start_time))

