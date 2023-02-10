from topGeneFold import * 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def posUpNegDownRound(x):
    if x >= 0:
        return int(round(x + 0.5))
    else:
        return int(round(x - 0.5))    

testList = ["T3","T5","T6","T9","T11","T12","T8", "C1"]
#testList = ["C1","T2","C7","T8"]

def sampCombos(samples):
    '''Accepts a list of samples and returns a list of all possible combinations of samples.'''
    combinations = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            string1 = samples[i]
            string2 = samples[j]
            combination1 = f"{string1[0]}_{string1[1:]}-{string2[0]}_{string2[1:]}"
            combination2 = f"{string2[0]}_{string2[1:]}-{string1[0]}_{string1[1:]}"
            combinations.append(combination1)
            combinations.append(combination2)
    return combinations

testList = sampCombos(testList)

def wantedVersuses(testList):
    '''Accepts a list of versuses and returns a heatmap with the logFC of the genes in the list'''
    #Read in the known risk genes
    Risk = pd.read_excel("KnownRiskGenes.xlsx")
    riskGeneList = Risk["riskGene"].tolist()
    TippyTop = pd.read_feather("/home/petear/Documents/TippyTopGeneDF_ALL.feather")
    #Add gene name column to tippytop
    mapfile = getEnsemblMappings()
    FullDF = pd.DataFrame()
    #remove the versuses not in the list from TippyTop
    TippyTop = TippyTop[TippyTop["Versus"].isin(testList)]#
    #remove the versuses not in the risk gene list from TippyTop
    for unq in TippyTop["Versus"].unique():
        df = topGeneFoldGeneName(unq, mapfile)
        df = df[0]
        df = df[df['GeneID'].isin(riskGeneList)]
        #merge df with FullDF
        FullDF = pd.concat([FullDF, df])   
        #make a dictionary with the values of the heatmap
    pivDF = FullDF.pivot("GeneID", "Versus", "logFC")
    #Finding max and min value of the pifDF and convert to int rounding up for max and down for min
    max = posUpNegDownRound(pivDF.max().max())
    min = posUpNegDownRound(pivDF.min().min())    
    #make a heatmap with nan colored white and the rest of the values colored by the color map
    cmap = sns.color_palette("magma", n_colors=128)
    sns.heatmap(pivDF, cmap=cmap, cbar_kws={'label': 'logFC'},annot_kws={"style": "italic"}, vmin=min, vmax=max, xticklabels=1, yticklabels=1, linewidths=0.5, linecolor="grey", square=True) 
    plt.show()

