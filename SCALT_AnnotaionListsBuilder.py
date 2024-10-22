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

''' The function lists_generator(counts,anno,boo,cells,genes,notation,cpus) is resposible of calling the utility of SCALT aimed to generate the custom collection 
    of cell-type specific lists of genes given the row counts table and the annotation of each cell. The function requires the following inputs: the row counts table;
    the table with one column reporting the annotation of each cell present in the counts; the number of boostraps sample; the number of cells to sample per each cell
    type during the boostrap process; the number of genes to include in the final cell-type specific lists; the kind of gene notation; the number of processors used.'''

def lists_generator(counts,anno,boo,cells,genes,notation,cpus):
    
    #### Create the directory that will contain the final lists derived from the annotation ####
    try:
        os.system("mkdir ./custom")
    except:
        pass

    #### Pre-process the inputs such to be suitable for the application run ####
    try:
        os.system("python3 scripts/inputPreparation.py "+counts+" "+anno+" anno "+notation+" None")   #input adjustment 
    except:
        print("Error: the input preparation program could not be executed!")
    
    name_counts = counts.split(".tsv")[0]
    counts = name_counts+"_adj.tsv"
    name_anno = anno.split(".tsv")[0]
    anno = name_anno+"_adj.tsv"
    
    #### Group cells by cell type annotation  ####
    groupped_dir = os.path.exists("groupped_cell_types")    #directory that will contain the counts of cells groupped by cell type
    if groupped_dir == False:
        try:
            os.system("mkdir ./groupped_cell_types")
        except:
            print('The directory "groupped_cell_types/" already exists')
    try:
        os.system("python3 scripts/cell_type_grouper.py "+counts+" "+anno+" "+cpus)   #group cells following the annotation table given as input
    except:
        print("Error in calling python3 scripts/cell_type_grouper.py!")
    
    #### Boostrap samples calculation ####
    boostrap_dir = os.path.exists("./boostraps_samples")    #directory that will contain the boostrap samples
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
                print("Error: cannot move probabilities table in the proper directory probabilities_tables/!")

            #### Move all data relative to the boostrap sample in the corresponding directory ####
            try:
                os.system("mv probabilities_tables/ TABLE_OF_GENES.tsv random_equilibrate_cell_type_table.tsv random_equilibrate_annotation.tsv boostraps_samples/"+str(b)+"_boostrap_sample")
            except:
                print("Error: cannot move data of boostrap number "+str(b)+" in the corresponding directory!")

    #### Generate the mean propabilities calculated from the probabilities table of the boostrap paradigm ####
    try:
        os.system("python3 scripts/statistics_calculator.py "+notation)
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
        os.system("python3 scripts/CTs_lists_generator.py genesRanking.tsv TABLE_OF_GENES.tsv genesProbabilities_ratios.tsv genes2remove.tsv "+genes+" anno")
    except:
        print("Error: cannot generate new list of genes specific for the cell types theoretically present in the sample!")
	
	#### Refine lists of genes reporting probabilities values for each gene ####
    try:
        os.system("python3 scripts/customLists_refinement.py "+notation+" anno")
    except:
        print("Error: cannot run the scripts/customLists_refinement.py")

    #### Check if the custom lists are similar to those available in SCALT using a mutual information strategy ####

    try:
        os.system("python3 scripts/mutual_information_test.py config/parametersForMIcalculation.json "+notation+" custom")
    except:
        print("Error: cannot run scripts/mutual_information_test.py")

    #### Group all the files in a directory ####
    try:
        os.system("mkdir ./AnnolistsBuilder_results")
    except:
        pass

    try:
        os.system("mv "+counts+" originalTables_zipped.zip MI_REPORT.txt boostraps_samples/ cellTypes_fromAnnotationHeatmap.png genes2remove.tsv genesCellTypes_probabilities.tsv genes_entropy.tsv genesGeneral_probabilities.tsv genesProbabilities_ratios.tsv genesRanking.tsv groupped_cell_types/ TABLE_OF_GENES.tsv "+anno+" AnnolistsBuilder_results/")
    except:
        print("Error: cannot move intermediate files in the ./AnnolistsBuilder_results directory")

''' Positional arguments '''
parser = argparse.ArgumentParser(description='SCALT: build the cell-type specific lists of genes starting from a counts matrix and correspoding annotation of each cell.')
parser.add_argument("Counts",metavar="Sample",help="Sample counts")
parser.add_argument("Anno",metavar="Annotation",help="Annotation of each cell present in the sample. It must be a tab separated file (.tsv) having N rows equal to the number of columns in the counts matrix (cells) and 1 column reporting the correspoding annotation.")
''' Optional arguments '''
parser.add_argument('-Boo',metavar='--Boostraps',default="100",help="Number of boostrap samples required to estimate probabilities. The default value is 100.")
parser.add_argument('-Cells',metavar='--Cells',default="100",help="Number cells to collect randomly per each cell type in each boostrap sample. The default value is 100.")
parser.add_argument('-Genes',metavar='--Genes',default="100",help="Number of genes that the lists will have at the end. The default value is 100.")
parser.add_argument("-Notation",metavar="--Notation",default="ensembl_id",help='Type of gene notation to use. The defaul is "ensembl_id". Write "gene_symbol" to switch to gene symbol notation.')
parser.add_argument("-CPUs",metavar="--CPUs",default="1",help='Number of processors employed. The default number is 1.')

args = vars(parser.parse_args())

lists_generator(args["Counts"],args["Anno"],args["Boo"],args["Cells"],args["Genes"],args["Notation"],args["CPUs"])

end_time = datetime.now()
print('Duration generation of cell type specific list of genes from annotation: {}'.format(end_time - start_time))

