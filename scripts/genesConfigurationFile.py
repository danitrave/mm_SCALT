#!/usr/bin/python

""" Libraries required  """

import sys
import os
import pandas as pd
import numpy as np
from numpy.linalg import norm
from datetime import datetime
import concurrent.futures
import itertools
import json
from operator import itemgetter

start_time = datetime.now()

''' The function geneNames_configurationFile(f) is used to generate a configuration file reporting the ensembl ids and the gene names of the genes+
	used in the tool. The input is the gene code annotation version 109 "Homo_sapiens.GRCh38.109.gtf" directly downloaded from the ftp fo GENECODE.'''

def geneNames_configurationFile(f):
    dd = {}
    gtf = open(f,"r")      #open the gtf
    for line in gtf:
        if line.startswith("#"):   #skip the header. Keep all kind of annotations (havana, ensembl, havana_ensembl)
            continue
        else:
            l = line.strip("\n").split("\t")    #remove "\n" and split by tab
            if l[2]=="gene":
                info = [i.split('"') for i in l[8].replace(' ',"").split(";")]
                d_info = {e[0]:e[1] for e in info if len(e)>=2}
                gene_biotype = d_info["gene_biotype"]
                if gene_biotype == "protein_coding":    #collect only protein coding genes for which there is a gene symbol (those genes having only the ensembl id are not collected)
                    if "gene_name" in d_info:
                        if d_info["gene_id"] not in dd:
                            dd[d_info["gene_id"]]=d_info["gene_name"]
                        else:
                            print(d_info["gene_id"])
    gtf.close()

    df = pd.DataFrame.from_dict({"ensembl_id":list(dd.keys()),"gene_symbol":list(dd.values())})   #generate the configuration file
    df.to_csv("GRCh38.109_ensembID_geneSymbol.tsv",sep="\t",header=True, index=True)

    configDir_check = os.path.exists("./config")     #if the config/ directory is not available yet, generate it
    if configDir_check == False:
        try:
            os.system("mkdir ./config")
        except:
            print("Error: cannor generate the directory for configuration files!")

        try:
            os.system("mv GRCh38.109_ensembID_geneSymbol.tsv config/")     #move the newly created configuration file in the proper directory
        except:
            print("Error: cannot move the configuration file GRCh38.109_ensembID_geneSymbol.tsv in the config/ directory")
    else:
        try:
            os.system("mv GRCh38.109_ensembID_geneSymbol.tsv config/")    #move the newly created configuration file in the proper directory
        except:
            print("Error: cannot move the configuration file GRCh38.109_ensembID_geneSymbol.tsv in the config/ directory")
    try:
        os.system("rm Homo_sapiens.GRCh38.109.gtf")    #remove the original gtf file
    except:
        print("Error: cannot remove the file Homo_sapiens.GRCh38.109.gtf")

''' The function download() simply downloades the gtf containing the annotation version 109 from GENECODE ftp site'''

def download():
    try:
        os.system("wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz")  #download
    except:
        print("Error: cannot download the Homo_sapiens.GRCh38.109.gtf.gz form the Ensembl Archive FTP directory")
    try:
        os.system("gunzip Homo_sapiens.GRCh38.109.gtf.gz")				#unzip
    except:
        print("Error: cannot unzip the file Homo_sapiens.GRCh38.109.gtf.gz")

if __name__ == "__main__":
    configExists = os.path.exists("config/GRCh38.109_ensembID_geneSymbol.tsv")
    if configExists == False:
        check_gtf = os.path.exists("Homo_sapiens.GRCh38.109.gtf")
        if check_gtf == False:					#if the gtf is not present in the directory, download it and generate the configuration file
            file_download = download()
        geneNames_configurationFile("Homo_sapiens.GRCh38.109.gtf")
    else:
        print('The configuration file "config/GRCh38.109_ensembID_geneSymbol.tsv" already exists')

end_time = datetime.now()
print('Time required for the creation of the configuration file: {}'.format(end_time - start_time))

