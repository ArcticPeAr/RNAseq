#########################################################################################################
##This script takes a list of up and down regulated genes and uses the fold change to determine the top 150 genes.
##Also the fold change will be kept in a matrix as "mFC" in the R-package "metaLINCS". This for analysis of perturgagens.  
#> head(mFC)
#        Resistant.vs.Sensitive WithaferinA.vs.Untreated UTR.vs.UTS WAR.vs.WAS
#CDKN2A                   6.200                   -0.036      6.150      6.250
#CTAG2                    5.655                   -0.232      5.716      5.594
# HMOX1                    0.602                    6.479      0.669      0.535
# NLRP11                   4.493                   -0.548      4.318      4.668
# MBL2                    -4.896                    0.092     -4.850     -4.943
# SLC43A3                  4.113                   -0.422      4.222      4.003
#         WAR.vs.UTR WAS.vs.UTS
# CDKN2A       0.013     -0.086
# CTAG2       -0.293     -0.170
# HMOX1        6.412      6.547
# NLRP11      -0.373     -0.723
# MBL2         0.046      0.139
# SLC43A3     -0.531     -0.312
#########################################################################################################

import pandas as pd
from biomart import BiomartServer

def reverseString(string):
    '''Accepts a string and returns the string in reverse order.'''
    parts = string.split("-")
    return "-".join(reversed(parts))


#From https://autobencoder.com/2021-10-03-gene-conversion/
def getEnsemblMappings():
    '''Connects to Ensembl Biomart and retrieves the mapping between Entrez Gene IDs and HGNC symbols. Returns a dict with the mapping.'''
    # Set up connection to server                                               
    server = BiomartServer("http://www.ensembl.org/biomart/martservice")
    mart = server.datasets["hsapiens_gene_ensembl"]                            
    # List the types of data we want                                            
    attributes = ["entrezgene_id", "hgnc_symbol", "ensembl_gene_id"]
    # Get the mapping between the attributes                                    
    response = mart.search({'attributes': attributes})                          
    data = response.raw.data.decode('ascii')
    entrezID2Name = {}
    entrezID2Ensembl = {}
    # Store the data in a dict                                                  
    for line in data.splitlines():
        line = line.split('\t')
        # The entries are in the same order as in the `attributes` variable
        entrezgene_id = line[0]
        hgnc_symbol = line[1]
        # Some of these keys may be an empty string. If you want, you can 
        # avoid having a '' key in your dict by ensuring the 
        # transcript/gene/peptide ids have a nonzero length before
        # adding them to the dict
        entrezID2Name[entrezgene_id] = hgnc_symbol
    for line in data.splitlines():
        line = line.split('\t')
        # The entries are in the same order as in the `attributes` variable
        entrezgene_id = line[0]
        ensembl_gene_id = line[2]
        # Some of these keys may be an empty string. If you want, you can 
        # avoid having a '' key in your dict by ensuring the 
        # transcript/gene/peptide ids have a nonzero length before
        # adding them to the dict
        entrezID2Ensembl[entrezgene_id] = ensembl_gene_id
    return entrezID2Name, entrezID2Ensembl 


#function to find Novogenes new csv-files to find the fold change.
def topGeneFoldGeneName(versus, mapfile):
    '''Accepts the versus as string and a file with mapping of and returns the logFCfile, upregulated genes as upRegDF and downregulated genes as downRegDF.'''
    #open the files
    #logFCfile
    logFCfile = pd.read_feather("/home/petear/Documents/TippyTopGeneDF_ALL.feather")
    #versus is a string with the comparison. Remove those not part of the versus.
    if versus not in logFCfile["Versus"].unique():
        versus = reverseString(versus)
    logFCfile = logFCfile[logFCfile["Versus"].str.contains(versus) == True]
    #Change the ENTREZID to string from int
    logFCfile["ENTREZID"] = logFCfile["ENTREZID"].astype(str)
    #Add gene name column to logFCfile
    logFCfile["GeneID"] = logFCfile["ENTREZID"].map(mapfile[0])
    #sort the logFCfile by the logFC
    logFCfile = logFCfile.sort_values(by="logFC", ascending=False)
    #split the logFCfile in up and down regulated genes
    upRegDF = logFCfile[logFCfile["logFC"] >= 0]
    downRegDF = logFCfile[logFCfile["logFC"] < 0]
    return logFCfile , upRegDF, downRegDF

def topGeneFoldEnsembl(versus, mapfile):
    '''Accepts the versus as string and a file with mapping of and returns the logFCfile, upregulated genes as upRegDF and downregulated genes as downRegDF.'''
    #open the files
    #logFCfile
    logFCfile = pd.read_feather("/home/petear/Documents/TippyTopGeneDF_ALL.feather")
    #versus is a string with the comparison. Remove those not part of the versus.
    if versus not in logFCfile["Versus"].unique():
        versus = reverseString(versus)
    logFCfile = logFCfile[logFCfile["Versus"].str.contains(versus) == True]
    #Change the ENTREZID to string from int
    logFCfile["ENTREZID"] = logFCfile["ENTREZID"].astype(str)
    #Add gene name column to logFCfile
    logFCfile["Ensemble"] = logFCfile["ENTREZID"].map(mapfile[1])
    #sort the logFCfile by the logFC
    logFCfile = logFCfile.sort_values(by="logFC", ascending=False)
    #split the logFCfile in up and down regulated genes
    upRegDF = logFCfile[logFCfile["logFC"] >= 0]
    downRegDF = logFCfile[logFCfile["logFC"] < 0]
    return logFCfile, upRegDF, downRegDF



def extractJustGenes(DF):
    '''Accepts a dataframe and returns a list of the genes.'''
    #Extract the genes
    genes = DF["GeneID"].tolist()
    return genes

def ExistInList(dataframe, list):
    '''Accepts a dataframe and a list and returns a dataframe with the genes in the list.'''
    #Extract the genes
    genes = extractJustGenes(dataframe)
    #Check if the genes are in the list
    genesInList = []
    for gene in genes:
        if gene in list:
            genesInList.append(gene)
    #Make a dataframe with the genes in the list
    genesInListDF = dataframe[dataframe["GeneID"].isin(genesInList)]
    return genesInListDF

def Export2Excel(dataframe, filename):
    '''Accepts a dataframe and a filename and exports the dataframe to an excel file.'''
    #Export the dataframe to an excel file
    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    dataframe.to_excel(writer, sheet_name='Sheet1')
    writer.save()
    return