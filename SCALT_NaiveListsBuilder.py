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

''' The function custumList_tsvGenerator(custom,notation) format the lists of genes inserted by the user in the input to be suitable for the naive classification. 
    The inputs are: the ".txt" file inserted at the beginning and the gene notation to use'''
    
def customList_tsvGenerator(custom,notation):
    
    #### Generate the directory that will contain the final lists generated from the naive classification ####
    try:
        os.system("mkdir ./naive")
    except:
        pass

    #### Generate the directory that will contain the processed user-custom lists of genes ####
    try:
        os.system("mkdir ./user_lists")
    except:
        print("Error: cannot create the directory user_lists/")
    #### Collect all the custom cell types in a dictionary ####
    d = {}
    txt = open(custom,"r")      #open the ".txt" file
    for line in txt:
        if line.startswith("-"):
            cell_type = line.strip("\n")[1:]
            d[cell_type]=[]
        else:
            gene = line.strip("\n")      #collect the genes in the corresponding cell type name of the dictionary
            d[cell_type]+=[gene,]
    txt.close()

    #### Generate a tsv per each custom cell type ####
    for ct in d:
        df = pd.DataFrame.from_dict({notation:d[ct]})
        name_file = ct.replace(" ",".").replace("–","-").replace("/","_").replace("(","_").replace(")","_").upper()+".tsv"   #name of the file for the corresponding list
        df.to_csv(name_file,sep="\t",header=True,index=True)        #generate the ".tsv" for each list and save in the directory userLists/
        try:
            os.system("mv "+name_file+" user_lists/")
        except:
            print("Cannot move "+name_file+" in the user_lists/ directory!") 
            
''' The function naive_lists_generator(counts,custom,boo,cells,genes,notation,cpus) is responsible of running the generation of the cell-type specific lists of genes
    starting from a collection of user-defined lists and employing a naive-based strategy. The inputs are: the counts table; the ".txt" file reporting the lists of 
    the user-defined cell types with corresponding cell type name; the number of boostraps sample; the number of cells to sample per each cell
    type during the boostrap process; the number of genes to include in the final cell-type specific lists; the kind of gene notation; the number of processors used. '''

def naive_lists_generator(counts,custom,boo,cells,genes,notation,cpus):
    #### Before everything, put the custom list of genes in the proper format and run the naive classifier ####
    custom_lists = customList_tsvGenerator(custom,notation)
    #### Pre-process the inputs such to be suitable for the following steps ####
    try:
        os.system("python3 scripts/inputPreparation.py "+counts+" None naive "+notation+" None")  #adjust the input counts
    except:
        print("Error: the input preparation program could not be executed!")

    name_counts = counts.split(".tsv")[0]
    counts = name_counts+"_adj.tsv"
    try:
       os.system("Rscript --vanilla scripts/hypergeometric.R "+counts+" "+genes+" "+notation)   #run the hypergeometric test
    except:
        print("Error: cannot run the naive classifier!")
   
    #### Group cells by naive-defined cell type annotation ####
    groupped_dir = os.path.exists("groupped_cell_types")
    if groupped_dir == False:
        try:
            os.system("mkdir ./groupped_cell_types")   #directory that will contain the counts of cells groupped by cell type
        except:
            print('The directory "groupped_cell_types/" already exists')
    try:
        os.system("python3 scripts/cell_type_grouper.py "+counts+" "+counts.split(".tsv")[0]+"_naive_annotation.tsv"+" "+cpus) #group cells following the annotation defined by the hypergeometric test
    except:
        print("Error in calling python3 scripts/cell_type_grouper.py!")

    #### Boostrap samples calculation ####
    boostrap_dir = os.path.exists("./boostraps_samples")  #directory that will contain the boostrap samples
    if boostrap_dir == False:
        try:
            os.system("mkdir ./boostraps_samples")
        except:
            print("Error: cannot create boostraps_samples/ directory!")
        for b in range(0,int(boo)):
            #### Generate current boostrap sample directory ####
            try:
                os.system("mkdir ./boostraps_samples/"+str(b)+"_boostrap_sample")
            except:
                print("Error: cannot generate sample number "+str(b)+" in the boostraps_samples/ directory!")

            #### Random cell picking for the current boostrap####
            try:
                os.system("python3 scripts/randomCells_selector.py "+cells+" "+cpus)
            except:
                print("Error: cannot run python3 scripts/randomCells_selector.py "+cells+"!")

            #### Probabilities tables creation calculated from the current boostrap random cells selection ####
            try:
                os.system("python3 scripts/probabilities_tables_generator.py random_equilibrate_cell_type_table.tsv random_equilibrate_annotation.tsv TABLE_OF_GENES.tsv")
            except:
                print("Error: cannot generate probabilities tables for boostrap number "+str(b)+"!")

            #### Collect all probabilites tables and insert them in a directory called "probabilities_tables" ####
            try:
                os.system("mkdir ./probabilities_tables")
            except:
                print("Error: cannot generate probabilities_tables/ directory!")

            try:
                os.system("mv cell_type_probabilities.tsv global_probabilities.tsv probabilities_ratio_matrix.tsv probabilities_tables/")
            except:
                print("Error: cannot move probabilities table in the proer directory probabilities_tables/!")

            #### Move all data relative to the boostrap sample in the corresponding directory ####
            try:
                os.system("mv probabilities_tables/ TABLE_OF_GENES.tsv random_equilibrate_cell_type_table.tsv random_equilibrate_annotation.tsv boostraps_samples/"+str(b)+"_boostrap_sample")
            except:
                print("Error: cannot move data of boostrap number "+str(b)+" in the corresponding directory!")

    #### Generate the mean propabilities calculated from the probabilities table of the boostrap paradigm ####
    try:
        os.system("python3 scripts/statistics_calculator.py")
    except:
        print("Error: cannot run calculate the probabilities tables!")

    #### Rank genes based on the probabilities ratios i.e. P(G|CT)/P(G) ####
    try:
        os.system("Rscript --vanilla scripts/genesRanker_byRatio.R")
    except:
        print("Error: cannot perform ranking of genes based on the probabilities ratios!")

    #### Gene entropy calculation and extraction of genes to remove before lists creation ####
    try:
        os.system("python3 scripts/entropy_calculator.py genesCellTypes_probabilities.tsv "+notation)
    except:
        print("Error: cannot calculate entropy of genes!")

    #### Generation of the new lists of genes and heatmap reporting the intersection among each pair of cell-type specific list of genes ####
    try:
        os.system("python3 scripts/CTs_lists_generator.py genesRanking.tsv TABLE_OF_GENES.tsv genesProbabilities_ratios.tsv genes2remove.tsv "+genes+" naive")
    except:
        print("Error: cannot generate new list of genes specific for the cell types theoretically present in the sample!")
        
    #### Refine lists of genes reporting probabilities values for each gene ####
    try:
        os.system("python3 scripts/customLists_refinement.py "+notation+" naive")
    except:
        print("Error: cannot run the scripts/customLists_refinement.py")

    #### Check if the custom lists are similar to those available in SCALT using a mutual information strategy ####
    try:
        os.system("python3 scripts/mutual_information_test.py config/parametersForMIcalculation.json "+notation+" naive")
    except:
        print("Error: cannot run scripts/mutual_information_test.py")

    #### Group all the files in a directory ####
    try:
        os.system("mkdir ./NaivelistsBuilder_results")
    except:
        pass

    try:
        os.system("mv MI_REPORT.txt originalTables_zipped.zip user_lists/ "+counts.split(".tsv")[0]+"_naive_annotation.tsv "+counts.split(".tsv")[0]+"_FDR_table.tsv boostraps_samples/ cellTypes_fromNaiveHeatmap.png genes2remove.tsv genesCellTypes_probabilities.tsv genes_entropy.tsv genesGeneral_probabilities.tsv genesProbabilities_ratios.tsv genesRanking.tsv groupped_cell_types/ TABLE_OF_GENES.tsv NaivelistsBuilder_results/")
    except:
        print("Error: cannot move intermediate files in the ./NaivelistsBuilder_results")

''' Positional arguments '''
parser = argparse.ArgumentParser(description='SCALT: Build the list of genes specific for every cell type starting from a collection of user-defined lists of genes.')
parser.add_argument("Counts",metavar="Sample",help="Sample counts")
parser.add_argument("Custom",metavar="CustomLists",help='''Name of the file containing the user-defined cell types. The file must present a ".txt" extension and 
														   must present only one column. Each cell type must be defined as follows: one raw defines the name of the 
														   cell type; the following raws will be the genes defining that cell type. If multiple cell types lists
														   are provided, they must be insert subsequentlially without leaving empty raws between two cell type lists.''')

''' Optional arguments '''
parser.add_argument('-Boo',metavar='--Boostraps',default="100",help="Number of boostrap samples required to estimate probabilities. The default value is 100.")
parser.add_argument('-Cells',metavar='--Cells',default="100",help="Number of cells to collect randomly per each cell type in each boostrap sample. The default value is 100.")
parser.add_argument('-Genes',metavar='--Genes',default="100",help="Number of genes that the lists will have at the end. The default value is 100.")
parser.add_argument("-Notation",metavar="--Notation",default="ensembl_id",help='Type of gene notation to use. The defaul is "ensembl_id". Write "gene_symbol" to switch to gene symbol notation.')
parser.add_argument("-CPUs",metavar="--CPUs",default="1",help='Number of processors employed. The default value is 1.')

args = vars(parser.parse_args())

naive_lists_generator(args["Counts"],args["Custom"],args["Boo"],args["Cells"],args["Genes"],args["Notation"],args["CPUs"])

end_time = datetime.now()
print('Duration generation of cell type specific list of genes starting from custom list of genes defined by the user: {}'.format(end_time - start_time))

