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

""" The function application_run(counts,gene_threshold,notation,dirCellTypes,cpus) is resposible of running the main program of SCALT which is the single cell
    classification based on a likelihood ratio method. The function requires the following inputs: the raw data counts table; the minimum number of genes that a
    cell must presnt to be classifed; the kind of gene notation used; the cell types used in the annotation process; the number of processors used. """

def application_run(counts,gene_threshold,notation,dirCellTypes,cpus,granularity):
    #### Set the level of granilarity of the cell type specific lists of genes ####
    if granularity == "low":
        dirCellTypes = "low_granularity_cell_types"

    #### Pre-process the input such to be suitable for the application run ####
    if counts.endswith("_adj.tsv"):
        name_counts = counts.split("_adj.tsv")[0]     #if the user is generating custum cell-type lists from either annotation or external lists, there is no need of addional adjustment 
        print("Note: counts do not need conversion!")
    else:
        try:
            os.system("python3 scripts/inputPreparation.py "+counts+" None likelihood "+notation+" "+dirCellTypes)   #adjust the input file
        except:
            print("Error: the input preparation program could not be executed!")

        name_counts = counts.split(".tsv")[0]
        counts = name_counts+"_adj.tsv"
    
    #### Run the program for the single cell annotation adopting the likelihood-ratio test ####
    try:
        os.system("python3 scripts/likelihood_ratio_test.py "+counts+" "+gene_threshold+" "+notation+" "+dirCellTypes+" "+cpus)
    except:
        print("Error: Run failed!")

    #### Generate the final report ####
    
    try:
        os.system("python3 reportGenerator.py p_values.tsv "+name_counts+"_adj_genesExpressed_filter.tsv "+counts+" "+notation+" "+dirCellTypes)
    except:
        print("Error: the report could not be generated!")

    #### Group all the results in a unique directry except the report ####

    try:
        os.system("mkdir ./results_directory")
    except:
        print("Error: could not generate the directory ./results_directory")
    
    try:
        os.system("mv barplot_cellTypesAboundance.html barplot_survivedCells.html deltas.tsv originalTables_zipped.zip p_values.tsv UMAP_2D.html UMAP_3D.html "+name_counts+"_adj_genesExpressed_filter.tsv "+counts+" results_directory/")
    except:
        print("Error: cannot move the files in the final results_directory/")


''' Positional arguments '''
parser = argparse.ArgumentParser(description='''SCALT: Single Cell Annotation Likelihood Tool. SCALT introduces a paradigm-shift for the analysis of scRNAseq data. 
								Cells are annotated to a specific type at individual level, by using a simple but elegant method based on maximum likelihood, without the need
								for clustering, dimensionality reduction or manual annotation. SCALT leverages a collection of 471 lists of cell-type specific genes, 
								constructed by extensive re-analysis of comprehensive and expert curated catalogues (HPA and DISCO)''')
parser.add_argument("Counts",metavar="Sample",help="Sample counts")     

'''' Optional arguments '''
parser.add_argument('-Min',metavar='--Threshold',default="250",help="Minimum number of genes that a cell must express to be classified. The default value is 250.")   
parser.add_argument("-Notation",metavar="--Notation",default="ensembl_id",help='Type of gene notation to use. The defaul is "ensembl_id". Write "gene_symbol" to switch to gene symbol notation. Note: gene notation must coincide with the one present in the counts table')
parser.add_argument("-Types",metavar="--Types",default="cell_types",help='Directory name containg the lists of the cell types to be used in the likelihood test. By default, only the pre-compiled lists (DISCO, HPA) are used. If the user wants to use only the custom ones generated from annotation, insert "custom". Finally, if the users wants to use only the custom ones from user-defined lists, insert "naive".')
parser.add_argument("-CPUs",metavar="--CPUs",default="1",help='Number of processors employed.')
parser.add_argument("-Granularity",metavar="--Granularity",default="high",help='Level of granularity of the annotation.')

args = vars(parser.parse_args())
application_run(args["Counts"],args["Min"],args["Notation"],args["Types"],args["CPUs"],args["Granularity"])

end_time = datetime.now()
print('Duration complete tool run: {}'.format(end_time - start_time))

