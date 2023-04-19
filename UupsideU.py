 
#This module is used to find the set and intersection of two sets of specific GO terms.

import pandas as pd

#read in GO groups from txt file
def readExtGoList(GOtxt):
    '''Reads a text file with one GO term per line and adds those to a list.'''
    GOlistFile = []
    with open(GOtxt, 'r') as f:
        for line in f:
            GOlistFile.append(line.strip())
    return GOlistFile

def unclusteredGo(GOgroups):
    '''This function takes a dataframe with grouped GO terms where columns are the names of the groups or clusters of groups and the values are the GO terms. It then returns a list of GO groups that are not a part of a cluster'''
    df = pd.read_excel(GOgroups, dtype = "str")
    unclusteredDict = {}
    clusteredDict = {}
    for col in df.columns:
        if df[col].isin(df.columns).any():
            groups = df[col].tolist()
            groups = [x for x in groups if x == x]
            clusteredDict[col] = groups
        else:
            groups = df[col].tolist()
            groups = [x for x in groups if x == x]
            unclusteredDict[col] = groups
    return unclusteredDict, clusteredDict

def removeNotInGoTXT(GOlistFile, DictOfGOterms):
    '''This function takes a dictionary of GO terms and a list of GO terms and returns a dictionary of GO terms that are in the list.'''
    DictOfGOterms2 = DictOfGOterms.copy()
    for key in DictOfGOterms2.keys():
        # If the key is not in the list of allowed keys, remove the key-value pair
        if key not in GOlistFile:
            del DictOfGOterms[key]
    return DictOfGOterms


#Dict 
# Godict {GOterm1: {sample1: (up+down), sample2: (up+down), sample1+sample2 (up+up+down+down)}, GOterm2: {sample1: (up+down), sample2: (up+own),sample1+sample2 (up+up+down+down)}}



def unionizedDict(DictOfGoTerms):
    '''This function accepts a dict where keys are GO terms and values are dicts of versus and tuple with two lists of up and down genes. It then returns a dict where the keys are the GO terms and the values are a dict where the keys are the versuses and the values are a list of mergend up and down genes. In addition the versuses will be a collected as a last key in the dictionary called "Union" and the values will be a list of all up and down regulated genes for that GO term.'''
    unionDict = {}
    for key in DictOfGoTerms.keys():
        unionDict[key] = {}
        for versus in DictOfGoTerms[key].keys():
            unionDict[key][versus] = DictOfGoTerms[key][versus][0] + DictOfGoTerms[key][versus][1]
        unionDict[key]["Union"] = unionDict[key][versus]
    return unionDict

