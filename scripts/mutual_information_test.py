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

'''The function json2dictionaryConverter(JSON) has the role of loading the json configuration file and transform it into a dictionry. The json contains
   all the necessary pieces of information to perform the mutual information test between the lists just created and those available in SCALT.'''

def json2dictionaryConverter(JSON):
    j = open(JSON,"r")       #open the JSON file
    l = json.load(j)              #load the JSON as a dictionary
    d = {}
    for k in l:
        ensembl_id = pd.DataFrame(l[k][2],index=l[k][0],columns=["ensembl_id"],dtype=float)  #table with the probabilities expressed as ensembl ids
        gene_symbol = pd.DataFrame(l[k][2],index=l[k][1],columns=["gene_symbol"],dtype=float)  ##table with the probabilities expressed as gene symbols
        ts = [i[1] for i in l[k][3]]
        thresholds =  pd.DataFrame(ts,index=["quartile","min","max","mean"],columns=["threshold"],dtype=float)  #mi thresholds
        inner_d = {"ensembl_id":ensembl_id,"gene_symbol":gene_symbol,"thresholds":thresholds}
        d[k]=inner_d
    j.close()
    return d      #return the json transformed into a dictionary

''' mi(L,config,nota) is responsible of calculating the mutual information between one list and each available in SCALT. Nota can be either "ensembl_id"
    or "gene_symbol" depending on the kind of gene notation employed in the lists. L is the list to controle while config is the parameters file.'''

def mi(L,config,nota):
    similar = []
    for ct in config:     #split genes between common and not common with respect to the current cell type list from SCALT
        common = []
        not_common = []
        to_compare = list(config[ct][nota].index)
        for g in L:
            if g in to_compare:
                common += [g,]
            else:
                not_common += [g,]
        #### MI shared genes between two lists ####
        MI_common = sum(np.array(config[ct][nota].loc[common,nota])*(np.log10(np.array(config[ct][nota].loc[common,nota])/0.0001)))
        #### MI not shared genes between two lists ####
        n = np.array([0.0001]*len(not_common))
        MI_notcommon = sum(n*(np.log10(n/0.0001)))
        MI = MI_common+MI_notcommon
        if MI >= config[ct]["thresholds"].loc["min","threshold"]:  #if the MI between the two lists is above the min threshold, the list is not new
            similar += [(ct,MI,list(config[ct]["thresholds"]["threshold"])),]
        else:
            similar += ["NONE",]

    return similar

''' The function listsSimilarityValidator(params,notation,mode) has to manage all the mi tests between all custom lists and those availble in SCALT.
    Additionally, it writes a report with the results.'''

def listsSimilarityValidator(params,notation,mode):
    check_path = os.path.exists("./"+mode)
    if check_path == False:       #check if the directory custom or naive exists
        print("Message: the directory ./"+mode+" does not exist")
        return ""
    ls = os.listdir("./"+mode)
    if len(ls) == 0:
        print("Message: the directory ./"+mode+" is empty!")    #check is the directory custom or naive is not empty
        return ""
    report = open("MI_REPORT.txt","w")    #report
    report.write("Comparison by mutual information of lists derived from either annotation or user-defined cell type specific list of genes."+"\n")
    report.write("Results:"+"\n")
    paths_lists = [mode+"/"+l for l in ls]
    for e in paths_lists:
        f = open(e,"r")
        compare = []
        for line in f:     #collect the genes of the current custom lisst
            a = line.strip("\n").split("\t")
            compare += [a[0],]
        f.close()
        test = mi(compare,params,notation)    #mi test
        outcome = []
        c = 0
        LIST = e.split("/")[1]
        for r in test:
            if r == "NONE":
                continue
            else:
                phrase = "The cell type retrieved from custom list "+LIST+" has MI of "+str(r[1])+" with list "+r[0]+" above the minimum MI threshold. Thresholds: min="+str(r[2][1])+" first quartile="+str(r[2][0])+" mean="+str(r[2][3])+" max="+str(r[2][2])
                outcome += [phrase,]       #phrase to write in the report if a list is similar to another in SCALT
                c += 1
        if c == 0:      #new potential list
            report.write("The cell type obtained from cutsom list "+LIST+" is potentially new in respect to those available in SCALT!"+"\n")
        else:
            for p in outcome:
                report.write(p+"\n")
    report.close()

if __name__ == "__main__":
    JSON = sys.argv[1]       #json configuartion file with all necessary elements for MI calculation
    notation = sys.argv[2]     #either "ensembl_id" or "gene_symbol"
    mode = sys.argv[3]           #either "custom" or "naive"
    params = json2dictionaryConverter(JSON)
    listsSimilarityValidator(params,notation,mode)

end_time = datetime.now()
print('Duration mutual information test: {}'.format(end_time - start_time))
