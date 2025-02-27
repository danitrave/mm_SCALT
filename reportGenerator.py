#!/usr/bin/python

""" Libraries required """

import plotly
import plotly.express as px
import sys
import os
import pandas as pd
import numpy as np
from numpy.linalg import norm
from datetime import datetime
from operator import itemgetter
import random
import umap
from umap import UMAP
import matplotlib.pyplot as plt
import seaborn as sns

''' Define the drawing stype and the size of the plots '''

sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})

start_time = datetime.now()

''' The function uncertainty_validator(ds,t) verifies wherever the most likely annotation is significantly categorizing the cell with respect to that
    ranking second in terms of significance. This happens if the first and the second differ for a threshold defined (t) likelihood'''

def uncertainty_validator(ds,t):
    deltaOfdeltas = abs((ds[0])-(ds[1]))
    if deltaOfdeltas > t: #9.23:
        return "unbiased"
    else:
        return "uncertain"

''' 'The function significance_validator(R) is used to collect the annotations that resulted significant after likelihood-ratio test. The function takes a list as input.
     The latter presents a series of tuples reporting the annotation as first element and the corresponding p-value from the likelihood-ratio test as second. Only those
     cell types having p-value lower than 1 are collected.'''

def significance_validator(R):
    L = []
    for j in R:
        if j[1] < 0.05:     #p-value significance validation
            L += [j,]
    return L

''' The function filter_barPlot(f) is used to generate a barplot reporting the number of cells that survived or not to the number of expressed genes filter.''' 

def filter_barPlot(f,name_counts):
    table = open(f,"r")
    d = {"Pass":0,"Not Pass":0}     #the number of cells that survived or not is collected in dictionary. If the cell is marked as "PASS", it passed the filter. "EXCLUDE" otherwise.
    for line in table:
        l = line.strip("\n").split("\t")
        if l[1]=="PASS":
            d["Pass"]+=1
        elif l[1]=="EXCLUDE":
            d["Not Pass"]+=1
        else:
            continue
    table.close()
    df = pd.DataFrame.from_dict(d,orient="index")    #dataframe for plotting
    df.columns=["Counts"]
    df.insert(0, "Class", list(df.index), True)
    fig = px.bar(df,x="Class",y="Counts",text_auto='.2s',color="Class",color_discrete_map={'Pass': '#2ca25f','Not Pass': '#fdbb84'},  #plot
    title="Number of cells that survived and not survived the genes expression filter")  
    fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
    fig.write_html(name_counts+"_barplot_survivedCells.html")        #save as an html file 
    return fig

""" The function cellTypes_barplot(T1,F1,DELTAS,THRESH) generates a barplot reporting the number of cell types identified per each cell type class."""

def cellTypes_barplot(T1,F1,DELTAS,THRESH,counts_data,mode,name_counts,output):
    tsv = open(counts_data,"r")
    stop = 0
    for line in tsv:
        if stop == 1:
            break
        cells_ids = line.strip("\n").split("\t")[1:]
        stop = 1
    tsv.close()
    #### columns of the table containing results ####
    primary_label = []
    multi_label = []
    alternative_labels = []
    #################################################
    annotation = []
    cell_ontology_annotation = []
    cell_ontology_color = ["Unclassified","Multiassigned"]
    df = pd.read_csv(T1,sep="\t",header=0,index_col=0)    #table reporting the p-values from the likelihood-ratio test
    INDEXES = list(df.index)
    COLS = list(df.columns)
    deltas = pd.read_csv(DELTAS,sep="\t",header=0,index_col=0)  #likelihood difference bewteen avrg type and cell type model
    surv = pd.read_csv(F1,sep="\t",header=0,index_col=0)   #table reporting the PASS-EXCLUDE notation based on the survival or not to the genes-expressed filter
    cellTypes_counts = {"Unclassified":0,"Multiassigned":0}
    ##### Cell ontology config loading #####
    COD = {}
    CO = open("config/ontologies_config_summary.tsv","r")
    for line in CO:
        if line.startswith("cell_type_name"):
            continue
        else:
            k = line.strip("\n").split(",")
            if k[1] not in COD:
                COD[k[1]]=k[2]
            else:
                continue
    CO.close()
    for i in range(len(INDEXES)):
        if surv.iloc[i,0]=="EXCLUDE":              #cells not passing the genes-expressed filer are annotated as "unclassified"
            cellTypes_counts["Unclassified"]+=1
            annotation+=["Unclassified",]
            cell_ontology_annotation += ["Unclassified",]
            primary_label+=["Unclassified",]
            multi_label+=["NA",]
            alternative_labels+=["NA",]
            continue
        else:
            alt = list(df.iloc[i,:])
            minimum = min(alt)
            deltas_cell = sorted(list(deltas.iloc[i,:]))
            if minimum > 0.05:
                cellTypes_counts["Unclassified"]+=1     #cells having the lowest p-value greater than 0.05 are annotated as "unclassified"
                annotation+=["Unclassified",]
                cell_ontology_annotation += ["Unclassified",]
                primary_label+=["Unclassified",]
                multi_label+=["NA",]
                alternative_labels+=["NA",]
                continue

            uncertainty = uncertainty_validator(deltas_cell,THRESH)   #cells which potentially have multiple plausible (significant)--> Multiassigned 
            #if uncertainty == "uncertain":
                #cellTypes_counts["Multiassigned"]+=1
                #annotation+=["Multiassigned",]
                #continue
            zipped = list(zip(COLS,alt))
            sortedZipped = sorted(zipped,key=itemgetter(1),reverse=False)  #sort annotations by p-values in increasing order
            upperSortedZipped = [(e[0].replace("_",".").replace(" ","."),e[1]) for e in sortedZipped]
            retain_significant = significance_validator(upperSortedZipped)       #collect only significant annotations
            retain_significant_annotation = [k[0] for k in retain_significant]
            if uncertainty == "uncertain":
                #### SPECIAL CASE: according to the deltas is multi assigned BUT only one annotation is significant ####
                if len(retain_significant_annotation)==1:
                    annotation+=[retain_significant_annotation[0].replace("."," "),]
                    if mode == "cell_types":
                        ontology_na = COD["CL"+retain_significant_annotation[0].split("-CL")[1]]+"-CL"+retain_significant_annotation[0].split("-CL")[1]
                        cell_ontology_annotation+=[ontology_na,]
                        if ontology_na not in cell_ontology_color:
                            cell_ontology_color += [ontology_na,]
                    multi_label+=["NO",]
                    primary_label+=[retain_significant_annotation[0].replace("."," ")]
                    alternative_labels+=["NA",]
                    if retain_significant_annotation[0].replace("."," ") not in cellTypes_counts:   #count the number of cells annotated to a specific cell type and collect the numbers into a dictionary having the cell type cathegory as keys
                        cellTypes_counts[retain_significant_annotation[0].replace("."," ")]=1
                    else:
                        cellTypes_counts[retain_significant_annotation[0].replace("."," ")]+=1
                    continue
                else:
                    cellTypes_counts["Multiassigned"]+=1
                    annotation+=["Multiassigned",]
                    cell_ontology_annotation += ["Multiassigned",]
                    multi_label+=["YES",]
                    primary_label+=[retain_significant_annotation[0].replace("."," "),]
                    alphabetic = sorted(retain_significant_annotation)
                    alternative_labels+=[",".join(alphabetic),]
                    continue
            annotation+=[retain_significant_annotation[0].replace("."," "),]
            if mode == "cell_types":
                ontology_na = COD["CL"+retain_significant_annotation[0].split("-CL")[1]]+"-CL"+retain_significant_annotation[0].split("-CL")[1]
                cell_ontology_annotation+=[ontology_na,]
                if ontology_na not in cell_ontology_color:
                    cell_ontology_color += [ontology_na,]
            multi_label+=["NO",]
            primary_label+=[retain_significant_annotation[0].replace("."," ")]
            alternative_labels+=["NA",]
            if retain_significant_annotation[0].replace("."," ") not in cellTypes_counts:   #count the number of cells annotated to a specific cell type and collect the numbers into a dictionary having the cell type cathegory as keys
                cellTypes_counts[retain_significant_annotation[0].replace("."," ")]=1
            else:
                cellTypes_counts[retain_significant_annotation[0].replace("."," ")]+=1

    df = pd.DataFrame.from_dict(cellTypes_counts,orient="index")    #set the dataframe for plotting
    df.columns = ["Counts"]
    df.insert(0,"Cell Type",list(df.index), True)
    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(df.shape[0])]   #set colors for each cathegory
    color_ontologies = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(cell_ontology_color))]
    df.insert(2,"Color",color,True)
    descrete_coloring = {}
    descrete_coloring_ontologies = {}
    for y in range(len(cell_ontology_color)):
        descrete_coloring_ontologies[cell_ontology_color[y]]=color_ontologies[y]
    descrete_coloring_ontologies["Unclassified"]="black"                        #unclassified cells are "black"
    descrete_coloring_ontologies["Multiassigned"]="#d9d9d9"

    for x in list(df.index):
        descrete_coloring[df.loc[x,"Cell Type"]]=df.loc[x,"Color"]
    descrete_coloring["Unclassified"]="black"                        #unclassified cells are "black"
    descrete_coloring["Multiassigned"]="#d9d9d9"
    fig = px.bar(df,x="Cell Type",y="Counts",text_auto='.2s',color="Cell Type",color_discrete_map=descrete_coloring,  #plot
    title="Number of cells annotated per each cell type class")  
    fig.update_traces(textfont_size=15, textangle=0, textposition="outside", cliponaxis=False)
    fig.update_layout(xaxis_tickangle=-45,legend_font_size=15,legend_title_font_size=15,legend_itemsizing='constant')
    fig.update_xaxes(tickfont_size=7)
    fig.write_html(name_counts+"_barplot_cellTypesAboundance.html")     #save the plot in an html file
    TABLE = pd.DataFrame.from_dict({"Cell_ID":cells_ids,"Prim.Lab":primary_label,"Multi_label":multi_label,"Alt.Lab":alternative_labels},orient='columns')
    TABLE.to_csv(output+"_annotation_table.tsv",sep="\t",header=True,index=True)
    return fig,annotation,descrete_coloring,cells_ids,cell_ontology_annotation,descrete_coloring_ontologies

''' The function umapPlot(t,a,note,color_dict) is resposinble of generating two umaps to show the annotation outcome in a grafical fashion. 
	Two kind of umaps will be generated: 2D and 3D maps. The colors of each cells are the same of those reported in the "barplot_cellTypesAboundance.html" file.
	Cells will be plotted following the expression defiving from the union (without replacement) of the cell types lists used for the likelihood-based annotation.'''

def umapPlot(t,a,note,color_dict,M,cells_ids,onto_anno,onto_color,name_counts):
    df = pd.read_csv(t,sep="\t",header=0,index_col=0)   #counts table
    anno = pd.DataFrame(a)                              #annotation of each cell from the likelihood-ratio annotation
    anno.columns=["CELL_ANNOTATION"]
    uniqueAnno = list(anno["CELL_ANNOTATION"].unique())
    color = [color_dict[f] for f in list(anno["CELL_ANNOTATION"])]   #assing each cell to the corresponding color based on the annotation outcome
    
    #### Save as before but for the ontology ####
    onto_cell_type = pd.DataFrame(onto_anno)
    onto_cell_type.columns=["CELL_ANNOTATION"]
    onto_unique_anno = list(onto_cell_type["CELL_ANNOTATION"].unique())
    color_sc_onto_anno = [onto_color[f] for f in list(onto_cell_type["CELL_ANNOTATION"])]
    
    #### consider only the genes extracted from the cell types that have at leat 50 cells annototated ####
    occurencesCounting = {}
    for o in list(anno["CELL_ANNOTATION"]):
        if o not in occurencesCounting:
            occurencesCounting[o]=1
        else:
            occurencesCounting[o]+=1
    cell_typesForUMAP = []
    for k in occurencesCounting:
        if occurencesCounting[k] >= 50:
            cell_typesForUMAP+=[k.replace(" ",".")+".tsv",]
    
    genes = []
    cells = os.listdir(M+"/")     #genes to consider to the umap creation
    for c in cells:
        if c in cell_typesForUMAP:
            table = pd.read_csv(M+"/"+c,sep="\t",header=0,index_col=None)
            g = list(table[note])
            for i in g:
                if i not in genes:    #union (without replacement) of the genes from the lists used for the likelihood-based classification
                    genes += [i,]
                else:
                    continue
        else:
            continue
    
    filtDF = df.loc[genes,:].T 
    umap_2d = UMAP(n_components=2)     #calculate 2D UMAP coordinates 
    results_data = filtDF.values
    proj_2d = pd.DataFrame(umap_2d.fit_transform(results_data))    #dataframe for the plot
    proj_2d["Annotation"]=list(anno["CELL_ANNOTATION"])
    proj_2d.columns  = ["UMAP_1","UMAP_2","Annotation"]
    fig_2d = px.scatter(proj_2d, x="UMAP_1", y="UMAP_2",color=proj_2d.Annotation,color_discrete_map=color_dict,title="UMAP-2D")#labels={color_dict[k]:k for k in color_dict},title="UMAP-2D")
    fig_2d.update_layout(legend_font_size=15,legend_title_font_size=15,legend_itemsizing='constant')   #plot
    fig_2d.update_traces(marker_size=6)
    fig_2d.write_html(name_counts+"_UMAP_2D.html")    #save the plot in an html file

    umap_3d = UMAP(n_components=3)      #calculate 3D UMAP coordinates
    proj_3d = pd.DataFrame(umap_3d.fit_transform(results_data))  #dataframe for the plot
    proj_3d["Annotation"]=list(anno["CELL_ANNOTATION"])
    proj_3d.columns  = ["UMAP_1","UMAP_2","UMAP_3","Annotation"]
    fig_3d = px.scatter_3d(proj_3d, x="UMAP_1", y="UMAP_2", z="UMAP_3",color=proj_3d.Annotation,color_discrete_map=color_dict,title="UMAP-3D")#labels={color_dict[k]:k for k in color_dict},title="UMAP-3D")
    fig_3d.update_layout(legend_font_size=15,legend_title_font_size=15,legend_itemsizing='constant')   #plot
    fig_3d.update_traces(marker_size=3.5)
    fig_3d.write_html(name_counts+"_UMAP_3D.html")   #save the plot in an html file
    proj_2d.index=cells_ids
    proj_3d.index=cells_ids
    proj_2d.to_csv(name_counts+"_umap_2d_coords.tsv",sep="\t",header=True, index=True)
    proj_3d.to_csv(name_counts+"_umap_3d_coords.tsv",sep="\t",header=True,index=True)

    #### Save the umap but using the cell ontology annotation ####
    if M == "cell_types":
        proj_2d["Annotation"]=list(onto_cell_type["CELL_ANNOTATION"])
        proj_2d.columns  = ["UMAP_1","UMAP_2","Annotation"]
        fig_2d_onto = px.scatter(proj_2d, x="UMAP_1", y="UMAP_2",color=proj_2d.Annotation,color_discrete_map=onto_color,title="UMAP-2D cell ontologies")
        fig_2d_onto.update_layout(legend_font_size=15,legend_title_font_size=15,legend_itemsizing='constant')   #plot
        fig_2d_onto.update_traces(marker_size=6)
        fig_2d_onto.write_html(name_counts+"_UMAP_2D_ONTO.html")
    else:
        fig_2d_onto="NULL"
    return fig_2d,fig_3d,fig_2d_onto
    
''' The function html_FinalReport(fig1,fig2,fig3,fig4) collects all the figures previously described and generate a final html report. '''

def html_FinalReport(fig1,fig2,fig3,fig4,fig5,output):
    output_name = output+"_report.html"	
    if fig5 == "NULL":
        with open(output_name, 'w') as html:
            html.write(fig1.to_html(full_html=False, include_plotlyjs='cdn'))
            html.write(fig2.to_html(full_html=False, include_plotlyjs='cdn'))
            html.write(fig3.to_html(full_html=False, include_plotlyjs='cdn'))
            html.write(fig4.to_html(full_html=False, include_plotlyjs='cdn'))
        html.close()
    else:
        with open(output_name, 'w') as html:
            html.write(fig1.to_html(full_html=False, include_plotlyjs='cdn'))
            html.write(fig2.to_html(full_html=False, include_plotlyjs='cdn'))
            html.write(fig3.to_html(full_html=False, include_plotlyjs='cdn'))
            html.write(fig4.to_html(full_html=False, include_plotlyjs='cdn'))
            html.write(fig5.to_html(full_html=False, include_plotlyjs='cdn'))
        html.close()

if __name__ == "__main__":
    p_values = sys.argv[1]       #table with all p-values
    expression = sys.argv[2]     #table that says if a cell has passed the genes-expression filter or not
    counts = sys.argv[3]         #counts table adjusted
    notation = sys.argv[4]       #either "ensembl_id" or "gene_symbol" 
    mode = sys.argv[5]           #either "cell_type", "custom" or "naive"
    deltas = sys.argv[6]         #table with delta values
    lik_threshold = float(sys.argv[7])  #likelihood difference threshold 
    output = sys.argv[8]    #name of the output table, report and directory	
    name_counts = counts.split(".tsv")[0]	
    survival_Barplot = filter_barPlot(expression,name_counts)
    cellTypeQuantity_Barplot = cellTypes_barplot(p_values,expression,deltas,lik_threshold,counts,mode,name_counts,output)
    umaps = umapPlot(counts,cellTypeQuantity_Barplot[1],notation,cellTypeQuantity_Barplot[2],mode,cellTypeQuantity_Barplot[3],cellTypeQuantity_Barplot[-2],cellTypeQuantity_Barplot[-1],name_counts)
    html_FinalReport(survival_Barplot,cellTypeQuantity_Barplot[0],umaps[0],umaps[1],umaps[2],output)
