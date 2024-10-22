#!/usr/bin/python

''' Required libraries '''

import sys
import os
import warnings
import argparse
from datetime import datetime
warnings.filterwarnings("ignore")

start_time = datetime.now()

''' The function converter(f) is responsible of converting the input file from either ".rds", ".RData" or "h5ad" to ".tsv" ''' 

def converter(f):
    if f.split(".")[-1]=="h5ad":        #if the file has the ".h5ad" format, call the program that converts it into ".tsv"
        try:
            os.system("python3 scripts/h5ad_format_converter.py "+f)
        except:
            pass
    elif f.split(".")[-1]=="RData" or f.split(".")[-1]=="rds":   #elif the file has either ".rds" or ".RData" format, call the program that converts it into ".tsv"
        try:
            os.system("Rscript --vanilla scripts/RDS_RData_formats_converter.R "+f)
        except:
            pass 

''' Positional arguments '''
parser = argparse.ArgumentParser(description='SCALT: auxiliary program for file format conversion. The program is able to convert files in ".h5ad", ".rds" and ".RData" to ".tsv" format.')
parser.add_argument("Counts",metavar="Sample",help="Sample counts")

args = vars(parser.parse_args())
converter(args["Counts"])

end_time = datetime.now()
print('Conversion to TSV format: {}'.format(end_time - start_time))
