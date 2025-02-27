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

def merge_run(D,notation,genes,boos,out):
    ####  MI calculation and threshold estimate ####
    try:
        os.system("python3 scripts/bootstrap_lists_mi.py AnnolistsBuilder_results/boostraps_samples "+genes+" "+boos+" "+notation)
    except:
        print("Error: unable to run the script bootstrap_lists_mi.py.")
    
    #### Estimate threshold for mutual information ####
    try:
        os.system("Rscript --vanilla scripts/thresholds_estimate_mi.R boos_mutual_information.tsv")
    except:
        print("Error: unable to run the script thresholds_estimate_mi.R.")

    #### Merge lists through mutual information ####
    try:
        os.system("python3 scripts/merger.py thresholdsMI.tsv miProbTables "+D+" "+out)
    except:
        print("Error: unable to run the script merger.py")
    
    #### Save all the results in a directory ####
    try:
        os.system("mkdir ./mi_merge_results")
    except:
        print("Error: unable to create the directory ./mi_merge_results.")
    try:
        os.system("mv boos_mutual_information.tsv thresholdsMI.tsv miProbTables/ mi/ experiment_MIs.tsv mi_merge_results/")
    except:
        print("Error: unable to move all temporary files in the directory mi_merge_results/")

''' Positional arguments '''
parser = argparse.ArgumentParser(description='''SCALT: Single Cell Annotation Likelihood Tool. Merge lists by mutual information.''')
parser.add_argument("Dir",metavar="Lists",help="Path of the cell type specific lists of genes to be merged by mutual information.")

''' Optional arguments '''
parser.add_argument("-Notation",metavar="--Notation",default="ensembl_id",help='Type of gene notation to use. The defaul is "ensembl_id". Write "gene_symbol" to switch to gene symbol notation. Note: gene notation must coincide with the one present in the counts table')
parser.add_argument("-Genes",metavar="--Genes",default="100",help='Number of genes that the lists will have at the end. The default value is 100.')
parser.add_argument("-Boo",metavar="--Boo",default="100",help='Number of bootstraps samples. The default value is 100.')
parser.add_argument("-out",metavar="--out",default="merging_diagnose.txt",help='Name for the output report.')
args = vars(parser.parse_args())
merge_run(args["Dir"],args["Notation"],args["Genes"],args["Boo"],args["out"])

end_time = datetime.now()
print('Duration merging lists by mutual information: {}'.format(end_time - start_time))

