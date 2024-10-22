#!/usr/bin/python

''' Libraries required '''

import sys
import os
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
from datetime import datetime
warnings.filterwarnings("ignore")

start_time = datetime.now()

''' The function h5ad2TSV(FILE) is responsible of extracting the counts in ".tsv" format from the ".h5ad" file. '''

def h5ad2TSV(FILE):
    name = FILE.split(".h5ad")[0]
    data = sc.read_h5ad(FILE)          #extract all the files from the h5ad file in a directory as ".csv" files
    data.T.write_csvs("H5DA", skip_data=False,sep="\t") 
    csv = open("H5DA/X.csv","r")							#open the counts table 
    
    df_genes = pd.read_csv("H5DA/obs.csv",sep="\t",index_col=0)      #put the gene ids and the cell ids in the counts table
    genes = list(df_genes.index)

    df_cells = pd.read_csv("H5DA/var.csv",sep="\t",index_col=0)
    cells = list(df_cells.index)

    f = open(name+".tsv","w")
    f.write("\t".join(["gene\cell_id"]+cells)+"\n")       #write the new ".tsv" reporting counts, gene ids and cell ids
    guard = 0
    for i in csv:
        new = "\t".join([genes[guard]]+i.strip("\n").split("\t"))
        f.write(new+"\n")
        guard += 1
    f.close()
    csv.close()

    try:
        os.system("rm -r H5DA/")    #remove the directory with all ".csv" files
    except:
        pass

if __name__ == "__main__":
    arg1 = sys.argv[1]
    if arg1.split(".")[-1]=="h5ad":    #check if the file has ".h5ad" format 
        h5ad2TSV(arg1)
    else:
        print('Error: the file is not in the ".h5ad" file.')

end_time = datetime.now()
print('Transformaton h5ad format to TSV: {}'.format(end_time - start_time))
