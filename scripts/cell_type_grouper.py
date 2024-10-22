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
from operator import itemgetter
import random
warnings.filterwarnings("ignore")

start_time = datetime.now()

''' The function cellType_grouper(INPUT) is used to group cells in ".tsv" files based on the annotation. The input is a list of tuples. 
	Each tuple has three elements: the cell type; the sample of origin and the bash-indexes of cells having that cell type annotation in that
	sample. At the end, the cells referred to the cell type present in th list will be groupped by cell type.'''

def cellType_grouper(INPUT):      #INPUT is a list with tuples each reporting the cell type, the sample and the corrected indexes of cells from that cell type
    global column_genes             #name of the tsv file reporting the column with the gene names
    temp = []
    for SET in INPUT:
        cell_types = {SET[0]:SET[2]}    #dictionary {cell type: [indexes of cells from that cell type]}
        sample = SET[1]
        CC = SET[0]   					#cell type
        for i in cell_types:
            if len(cell_types[i]) < 10000:                   #if I have less than 10000 cells, extract them in one step
                indexes = ",".join([w for w in cell_types[i]])        #join the indexes (strings) by colon
                temp += [CC+"_temp_"+sample,]
                try:
                    os.system("cut -f"+indexes+" "+sample+" >> "+CC+"_temp_"+sample)   #use the "cut" command and the indexes to subselect the cells of that cell type
                except:
                    print("Error: The number of cells for cell type "+CC+" in "+sample+" its too long.")
            else:       	#if there are more than 10000 cells, "cut" cannot handle the number of characters multiple, so the extraction requires more steps
                temp += [CC+"_temp_"+sample,]      #name of the temporary files containing the split counts from the current cell type from the current sample. These files will be pasted.
                j = 0
                slicing = []
                while j < len(cell_types[i]): 	#subdivide the indexes of cells in multiple lists each used to generate a tsv of cells of that cell type that will be merged
                    if j+10000 > len(cell_types[i]):            #check when we arrive at the end of the table
                        slicing += [cell_types[i][j:len(cell_types[i])],]
                        j += 10000
                    else:
                        slicing += [cell_types[i][j:j+10000],]
                        j += 10000
                        
                tables2merge = []              #list that will contain the names of the tsv files reporting the slices of the current cell type from the current sample
                for s in range(len(slicing)):
                    indexes = ",".join([d for d in slicing[s]])  #join the indexes (strings) by colons
                    name = CC+"_"+sample+".table"+str(s)+".tsv"
                    tables2merge += [name,]
                    try:
                        os.system("cut -f"+indexes+" "+sample+" > "+name)   #extract the group of cells of this cell type
                    except:
                        print("Error: The number of cells for cell type "+CC+" in "+sample+" its too long.")

                merged = " ".join(tables2merge)
                try:
                    os.system("paste "+merged+" > "+CC+"_temp_"+sample)    #paste the tsv of cell of that cell type in a unique tsv
                except:
                    print("Error: Cannot merge temporary tables for "+CC)

                try:
                    os.system("rm -r "+merged)     #remove all the tsv file having cells of the same kind that were merged before.
                except:
                    print("Error: Cannot remove temporary files.")
        
    for t in temp:
        try:
            os.system("paste column_one.tsv "+t+" > "+t.split("_temp_")[0]+".tsv")    #paste the column with the genes at the beginning of the counts table
        except:
            print("Error: cannot past the column with genes on the table "+t.split("_temp_")[0]+".tsv")
        try:
            os.system("mv "+t.split("_temp_")[0]+".tsv groupped_cell_types/")    #move the groupped cells in the proper directory
        except:
            print("Error: cannot move the file "+t.split("_temp_")[0]+".tsv in the directory groupped_cell_types/")
        try:
            os.system("rm -r "+t)
        except:
            print("Error: cannot remove file "+t)

''' The function index_correction(D) is responsible of adjusting the index of cells in order to make it suitable for bash coding.
	Python starts from index 0 while bash from 1 but since we have to exclude the column of the genes (first one) 2 will be summed to each index. '''

def index_correction(D):
    DD = {}
    for x in D:
        L = [str(j+2) for j in D[x]]     #indexes must be converted to strings because the call on bash using the command line requires strings.
        DD[x]=L
    return DD

''' The function parallelizer(I) has the role of parallelize the cell type grouping procedure on multiple processors. '''

def parallelizer(I):
    global cpus          #number of CPUs
    L = []
    for key in I:
        l = []
        for ele in I[key]:
            l += [(key,ele[0],ele[1]),]    #organize inputs such to be a list of tuples reportin [(cell type,sample,corrected indexes),]
        L += l
    if cpus == 1:
        return cellType_grouper(L)      #one CPU
    else:
        cells = len(I)
        guard = len(L)
        slicer = int(len(L)/cpus)
        slicing = []
        s = 0
        if cells <= cpus:                #check: the number of cell types is smaller than the CPUs inserted
            slicer = 1
        while s < guard:                  #if multiple CPUs are used, slice the input pre-organized such to group a maximum number of cell types defined by the number of precessors at each itaration
            if s+slicer > guard:
                slicing+=[[L[s:guard]],]
                s += slicer
            else:
                slicing+=[[L[s:s+slicer]],]
                s += slicer
        for process in slicing:
            with concurrent.futures.ProcessPoolExecutor() as executor:    #parallelizing
                output = executor.map(cellType_grouper,process)

if __name__ == "__main__":
    counts = sys.argv[1]     #counts table
    annotation = sys.argv[2] #annotation table. It can come from either annotation or naive
    cpus = int(sys.argv[3])  #number of CPUs
    groupped = {}
    #### Generate the genes column ####
    try:
        os.system("cut -f1 "+counts+" >> column_one.tsv")
    except:
        pass
    #### Group cells based on the annotation ####
    anno = pd.read_csv(annotation,sep="\t",header=0,index_col=None)
    cell_types = anno["CELL_ANNOTATION"].unique()
    index_dict = {}
    for c in cell_types:         #per each cell type, collect the index of the cells that are annotated as such
        Is = anno.index[anno['CELL_ANNOTATION'] == c].tolist()
        index_dict[c]=Is
    #### Filter cells based on aboundance of cell types and correct indexes for bash coding ####
    if "unknown" in index_dict:
        index_dict.pop("unknown")     #remove unknon cells from the hypergeometric test
    filtered_index_dict = {}
    for K in index_dict:             #eliminate those cell types found only once
        if len(index_dict[K])>1:
            filtered_index_dict[K]=index_dict[K]   #meno brutto --> –
    second_filtered_index_dict = {}
    for X in filtered_index_dict:               #adjust cell type name to avoid special characters that might disturb the computation
        correct_k = X.replace(" ",".").replace("–","-").replace("/","_").replace("(","_").replace(")","_").upper()
        second_filtered_index_dict[correct_k]=filtered_index_dict[X]
        
    reindexed = index_correction(second_filtered_index_dict)     #correct indexes such to be suitable for bash coding
    for ct in reindexed:
        if ct not in groupped:
            groupped[ct]=[(counts,reindexed[ct]),]    #if the same cell type is present in more than one samples, add it
        else:
            groupped[ct]+=[(counts,reindexed[ct]),]
    #### Run the program. Before, check if parallelization is required ####
    listBuilding_run = parallelizer(groupped)
    os.system("rm -r column_one.tsv")

end_time = datetime.now()
print('Duration cell type cells grouping: {}'.format(end_time - start_time))

