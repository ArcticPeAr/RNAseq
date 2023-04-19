 #This module is used to search for specific genes in compana and return the found genes with the corresponding information such as versus, fold change and p-value.

import pandas as pd

def readGeneFile(geneFile):
    '''This function reads in users gene list'''
    geneList = []
    for line in geneFile:
        line = line.strip()
        geneList.append(line)
    return geneList

def searchCompana(geneList):
    '''This function searches for the genes in the gene list and returns a pandas dataframe with the found genes as columns and the corresponding information as rows'''
    geneDict = {}
    for gene in geneList:
        geneDict[gene] = []
        
