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

''' The function random_cell_type_table_creation(INPUT) is used to generate a counts table presenting an equal number of cells from different
	cell types picked randomply. The input is a list of tuples. Each tuple has three elements: the cell type; the sample of origin and 
	the bash-indexes of cells having that cell type annotation in that sample.'''

def random_cell_type_table_creation(INPUT):      #INPUT is a list with tuples each reporting the cell type, the sample and the corrected indexes of cells from that cell type
    global numberOfCells								#number of cells per each cell types to pick randomly
    sample = "groupped_cell_types/"+INPUT[0]+".tsv"       #name of the sample from the groupped by cell type tsv counts
    temp = []
    equi_annotation = []							#annotation of the current randomly picked cells
    cell_types = {INPUT[0]:INPUT[1]}           #{cell type: [corrected indexes]}
    for i in cell_types:							#this for loop is required is one processor will handle multiple cell types
        equi_annotation = equi_annotation + [i]*len(cell_types[i])
        if len(cell_types[i]) < 10000:                   #if I have less than 10000 cells, extract them in one step
            CC = i.replace(" ","_").replace("/","_")
            temp += [CC+".table.tsv",]
            indexes = ",".join([w for w in cell_types[i]])        #join the indexes (strings) by colon 
            try:
                os.system("cut -f"+indexes+" "+sample+" >> "+CC+".table.tsv")   #use the "cut" command and the indexes to subselect the cells of that cell type
            except:
                print("Error: The number of columns is too large for cut to manage.")

        else:          #if there are more than 10000 cells, "cut" cannot handle the number of characters, so the extracted cells must be divided in multiple tsv and then merged
            CC = i.replace(" ","_").replace("/","_")
            temp += [CC+".table.tsv",]
            j = 0
            slicing = []
            while j < len(cell_types[i]):         #subdivide the indexes of cells in multiple lists each used to generate a tsv of cells of that cell type that will be merged
                if j+10000 > len(cell_types[i]):            #check when we arrive at the end of the table
                    slicing += [cell_types[i][j:len(cell_types[i])],]
                    j += 10000
                else:
                    slicing += [cell_type[i][j:j+10000],]    #slicing 
                    j += 10000
            tables2merge = []
            for s in range(len(slicing)):
                indexes = ",".join(d for d in slicing[s])  #join the indexes (strings) by colons
                name = CC+".table"+str(s)+".tsv"
                tables2merge += [name,]
                try:
                    os.system("cut -f"+indexes+" "+sample+" > "+name)   #extract the group of cells of this cell type
                except:
                    print("Error: The number of columns is too large for cut to manage.")

            merged = " ".join(tables2merge)
            try:
                os.system("paste "+merged+" > "+CC+".table.tsv")    #paste the tsv of cell of that cell type in a uniq tsv
            except:
                print("Error: Cannot merge the tables.")

            try:
                os.system("rm -r "+merged)     #remove all the tsv file having cells of the same kind that were merged before.
            except:
                print("Error: Cannot remove the temporary table files.")

    return equi_annotation

""" The function index_correction(D) is responsible of adjusting the index of cells in order to make it suitable for bash coding.
Python starts from index 0 while bash from 1 but since we have to exclude the column of the genes (first one) two will be summed
to each index. """

def index_correction(D):
    DD = {}
    for x in D:
        L = [str(j+2) for j in D[x]] 
        DD[x]=L
    return DD
    
''' The function parallelizer(I) has the role of parallelize the cell type grouping procedure on multiple processors. '''

def parallelizerPrior_run(reindexed,cpus):
    if cpus == 1:                              #one CPU
        final_annotation = []
        for k in reindexed:
            sample = (k,reindexed[k])           #(cell type,corrected indexes)
            run = random_cell_type_table_creation(sample)   #sample the cells and return the corresponding annotation 
            for o in run:
                final_annotation += [o,]       #put the annotation of all cells in a unique table respecing the order
        return final_annotation
    else:
        cells = list(reindexed.keys())         #more than CPUs
        slicer = int(len(cells)/cpus)
        slicing = []
        s = 0
        if len(cells) <= cpus:                #check: the number of cell types is smaller than the CPUs inserted
            slicer = 1
        while s < len(cells):                  #slice the cells types such to have the same number of cell types (equal to the number of CPUs) in each iteration 
            if s+slicer > len(cells):
                slicing+=[cells[s:len(cells)],]   #check if we are at the end of the cell types 
                s += slicer
            else:
                slicing+=[cells[s:s+slicer],]
                s += slicer
        to_parallel = []
        for e in slicing:
            d2_insert = []
            for c in e:
                d2_insert+=[(c,reindexed[c]),]    #prepare the inputs to be suitable for the random picking function
            to_parallel += [d2_insert,]

        final_annotation = []
        for process in to_parallel:
            with concurrent.futures.ProcessPoolExecutor() as executor:    #parallelizing
                run = executor.map(random_cell_type_table_creation,process)
                for o in run:
                    final_annotation += o
        return final_annotation

if __name__ == "__main__":
    ls = os.listdir("./groupped_cell_types")     #tsv from which the random picking will be performed
    numberOfCells = int(sys.argv[1])             #number of cells to pick from each cell type
    CPUs = int(sys.argv[2])						 #number of CPUs
    #### Generate the genes column ####
    try:
        os.system("cut -f1 groupped_cell_types/"+ls[0]+" >> column_one.tsv")    #genes column
    except:
        print("Error: column with genes cannot be created.")
    index_dict = {}
    #### Get the number of cells groupped per each cell type ####
    for ty in ls:
        cell_name = ty.split(".tsv")[0]
        test = pd.read_csv("groupped_cell_types/"+ty,chunksize = 1,sep="\t",index_col=0)
        df = next(test)
        cells = [i for i in range(df.shape[1])]    #python indexes of cells from the current cell type present in the tsv presenting the cells groupped by type
        index_dict[cell_name]=cells
    #### Filter those cell types found just once and correct indexes such to be suitable for bash coding ####
    filtered_index_dict = {}
    for K in index_dict:             #eliminate those cell types that the naive classifier found only once
        if len(index_dict[K]) > 1:
            filtered_index_dict[K]=index_dict[K]
    #### Random sampling #####
    random_sample = {}
    for s in filtered_index_dict:
        if len(filtered_index_dict[s]) < numberOfCells:
            random_sample[s]=filtered_index_dict[s]
        else:
            S=random.sample(filtered_index_dict[s], k=numberOfCells)
            random_sample[s]=S
    
    reindexed = index_correction(random_sample)     #correct indexes such to be suitable for bash coding

    #### Parallelizer and final tables generation####
    results = parallelizerPrior_run(reindexed,CPUs)
    cells = list(reindexed.keys())
    temp = [i.replace(" ","_").replace("/","_")+".table.tsv" for i in cells]      #name of the files containg the randomly samples cells per cell type
    
    try:
        os.system("paste  column_one.tsv "+" ".join(temp)+" > random_equilibrate_cell_type_table.tsv") #table with counts
    except:
        print("Error: the table random_equilibrate_cell_type_table.tsv cannot be created.")
    try:
        os.system("rm -r "+" ".join(temp))
    except:
        print("Error: cannot remove the temporary files.")
    
    finalEquiAnnotation = pd.DataFrame(results)
    finalEquiAnnotation.columns = ["CELL_ANNOTATION"]
    finalEquiAnnotation.to_csv("random_equilibrate_annotation.tsv",sep="\t",header=True,index=False)  #table with annotation
    try:
        os.system("cut -f1 random_equilibrate_cell_type_table.tsv > genes_column.tsv")   #save table of genes to be used in the following steps
    except:
        pass
    
    dfGENES = pd.read_csv("genes_column.tsv",sep="\t",header=0,index_col=False)
    GENES = list(dfGENES["genes/cells"])  #genes
    DF = pd.DataFrame(GENES)
    DF.columns = ["genes"]
    DF.to_csv("TABLE_OF_GENES.tsv",sep="\t")  #save the table of genes
    os.system("rm -r genes_column.tsv")
    os.system("rm column_one.tsv")

end_time = datetime.now()
print('Duration equilibrate read counts and annotation table generation: {}'.format(end_time - start_time))
    
    
