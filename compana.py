 ####Library imports
import pandas as pd
import matplotlib.pyplot as plt
import BarPlotModule as bpp
import sys
import os

######################################################################################################
####Welcome message:
######################################################################################################
print("Welcome to the Compana!\n This program is still under development, but here is a quick rundown:\n")
print("Type 'help' for a list of commands.\n")
######################################################################################################
####Create folder hierarchy: 
######################################################################################################

# get the current working directory
cwd = os.getcwd()
# get the parent directory
parentDir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

# point to input folder in parent directory
inputFolder = os.path.join(parentDir, "CompFiles")

#pd cus both for pandas and parental directory
pdFolder = os.path.join(parentDir)

# point to output folder in parent directory if it exists
outputFolder = os.path.join(parentDir, "Output")

######################################################################################################
####Functions in Main:
######################################################################################################

def readExternalList(txtFile):
    '''Reads a text file with one GO term per line and adds those to a list.'''
    listFromFile = []
    with open(txtFile, 'r') as f:
        for line in f:
            listFromFile.append(line.strip().lower())
    return listFromFile

def readExternalTSV(tsvFile):
    '''Accept a csv file with multiple columns and returns a list of all values from the columns.'''
    df = pd.read_csv(tsvFile, dtype = "str", sep='\t')
    data = df.values.tolist()
    flatList = [item for sublist in data for item in sublist]
    return flatList

def defineGroupsFromTSV(tsvFile):
    '''Accepts a TSV file and returns a dictionary where the column names are keys and the corresponding rows are values as a list..'''
    df = pd.read_csv(tsvFile, sep='\t')
    sampGr = df.to_dict('list')
    return sampGr


import re
def add_underscore(s):
    s = re.sub(r'([a-zA-Z])([0-9])', r'\1_\2', s)
    return s

def pairItems(clist):
    pairs = []
    for i in range(len(clist)):
        for j in range(i+1, len(clist)):
            # Add underscores if needed
            if "_" not in clist[i]:
                clist[i] = clist[i][:1] + "_" + clist[i][1:]
            if "_" not in clist[j]:
                clist[j] = clist[j][:1] + "_" + clist[j][1:]
            # Create pair and add to list
            pairs.append(clist[i] + "-" + clist[j])
            pairs.append(clist[j] + "-" + clist[i])
    return pairs

def reverseString(s):
    parts = s.split('-')
    return parts[1] + '-' + parts[0]
######################################################################################################
####Find samples names DOWN:
######################################################################################################

bddFile = pdFolder + "/BeDF_DOWN_perGO.csv"
cddFile = pdFolder + "/CeDF_DOWN_perGO.csv"
mddFile = pdFolder + "/MeDF_DOWN_perGO.csv"

bdd = pd.read_csv(bddFile, dtype = "str")
bdd.columns = bdd.columns.str.lower()
cdd = pd.read_csv(cddFile, dtype = "str")
cdd.columns = cdd.columns.str.lower()
mdd = pd.read_csv(mddFile, dtype = "str")
mdd.columns = mdd.columns.str.lower()

framesDown = [bdd, cdd, mdd]
down = pd.concat(framesDown)
downSamp = list(down)

downStrList = []

for st in downSamp:
    tt = st.split(' ', 1)[0]
    downStrList.append(tt)

downStrList2 = []
for st2 in downStrList:
    tt = st2.split('-', 1)
    downStrList2.append(tt)

flatDownList = [item for sublist in downStrList2 for item in sublist]
optionsVerDown = list(dict.fromkeys(flatDownList))


######################################################################################################
####Find samples names UP:
######################################################################################################

bduFile = pdFolder + "/BeDF_UP_perGO.csv"
cduFile = pdFolder + "/CeDF_UP_perGO.csv"
mduFile = pdFolder + "/MeDF_UP_perGO.csv"

bdu = pd.read_csv(bduFile, dtype = "str")
bdu.columns = bdu.columns.str.lower()
cdu = pd.read_csv(cduFile, dtype = "str")
cdu.columns = cdu.columns.str.lower()
mdu = pd.read_csv(mduFile, dtype = "str")
mdu.columns = mdu.columns.str.lower()

framesUp = [bdu, cdu, mdu]
up = pd.concat(framesUp)
upSamp = list(up)

upStrList = []

for st in upSamp:
    tt = st.split(' ', 1)[0]
    upStrList.append(tt)

upStrList2 = []
for st2 in upStrList:
    tt = st2.split('-', 1)
    upStrList2.append(tt)

flatupList = [item for sublist in upStrList2 for item in sublist]
optionsVerUp = list(dict.fromkeys(flatupList))

######################################################################################################
####Find common for up and down:
######################################################################################################
downStrListSansDUP = list(set(downStrList))
upStrListSansDUP = list(set(upStrList))
StrListStart = upStrListSansDUP + downStrListSansDUP
StrList = [valstr for valstr in StrListStart if StrListStart.count(valstr) > 1]
StrList = list(set(StrList))
#StrList = [x.upper() for x in StrList]

optionsVerDup = optionsVerUp + optionsVerDown
optionsVer = [val for val in optionsVerDup if optionsVerDup.count(val) > 1]
optionsVer = list(set(optionsVer))

#add "ALL" to optionsVer
optionsVer.append("all")

#
######################################################################################################

# bda = pd.read_csv("BeDF_FULL_perGO.csv")
# cda = pd.read_csv("CeDF_FULL_perGO.csv")
# mda = pd.read_csv("MeDF_FULL_perGO.csv")

# framesAll = [bdd, cdd, mdd]
# down = pd.concat(framesAll)

######################################################################################################
#### Starting Versus selection:
######################################################################################################
versList = []

# New version will consist of a TSV file with versus names in different columns
sampListFile = inputFolder + "/Versuses.tsv"
sampList = readExternalTSV(sampListFile)

#If user wants to make specific versus selection, but mistook the order of the samples
for item in sampList:
    if len(item) > 4:
        versList.append(item)
        item = reverseString(item)
        versList.append(item)

readVersList = pairItems(sampList)

    
for item in readVersList:
    if item in StrList:
        versList.append(item)

#lowercaps for every item in versList:
versList = [x.lower() for x in versList]

versList = [x for x in versList if x in StrList]


######################################################################################################
#### Flaggies!
######################################################################################################

if __name__ == '__main__':
    lists = False
    GOterm = False
    GOlists = False
    FC = False
    Perts = False
    Venn = False
    gene = False
    FCperVersus = False
    FCperGOterm = False
    UnionIntersect = False
    All = False
    help = False

for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-lists':
        lists = True
    elif sys.argv[i] == '-GO':
        GOterm = True
    elif sys.argv[i] == '-GOlists':
        GOlists = True
    elif sys.argv[i] == '-FC':
        FC = True
    elif sys.argv[i] == '-FClists':
        FClists = True
    elif sys.argv[i] == '-Perts':
        Perts = True
    elif sys.argv[i] == '-Venn':
        Venn = True
    elif sys.argv[i] == '-gene':
        gene = True
    elif sys.argv[i] == '-FCperVersus':
        FCperVersus = True
    elif sys.argv[i] == '-FCperGOterm':
        FCperGOterm = True
    elif sys.argv[i] == '-UnIsect':
        UnionIntersect = True
    elif sys.argv[i] == '-all':
        All = True
    elif sys.argv[i] == '-help' or sys.argv[i] == '-h':
        help = True

if len(sys.argv) == 1:
    help = True

if UnionIntersect == True:
    lists = True
    GOlists = True
if GOlists == True or Venn == True:
    GOterm = True
if Perts == True:
    FC = True
if gene == True:
    FC = True

   
############################################################################################
#Help  
############################################################################################
if help == True:
    print("Usage: python3 compana.py -[flag]")
    print("Flags:")
    print("-lists: Create lists of genes up and down regulated for each versus")
    print("-GO: Create GO term analysis for each versus")
    print("-GOlists: Create lists of genes up and down regulated for each GO term")
    print("-FC: Add fold change for each versus and GO term")
    print("-Perts: Create FC analysis for each perturbation")
    print("-Venn: Create Venn diagrams for each versus")
    print("-gene: Search for specific genes in the results")
    print("-all: Run all the above")
    print("-help: Print this help")
    print('ATTENTION: Because "Perts" is a flag that requires the use of external resources (https://maayanlab.cloud/L1000CDS2/query), it is not included in the "all" flag. If you want to use it along the "all" flag you can simply add it as well (E.g. python3 compana.py -all -Perts')
    sys.exit()
if help == False:
    print(sampList)
    print("Versuses selected: ", versList)
    print(StrList)

############################################################################################
#Versus list
############################################################################################

versListDict = {}
versPDDict = {}
for versus in versList:
    down2 = down.filter(regex=versus)
    up2 = up.filter(regex=versus)
    geneListUp = list()
    geneListDown = list()
    versPDDict[versus] = [up2, down2]       #Dict createad to be used in GO term analysis
#DOWN
    for (columnName, columnData) in down2.items():
        da = columnData.tolist()
        list_no_nan = [x for x in da if pd.notnull(x)]
        geneListDown.append(list_no_nan)
    flatDown = [item for sublist in geneListDown for item in sublist]
#UP
    for (columnName, columnData) in up2.items():
        da = columnData.tolist()
        list_no_nan = [x for x in da if pd.notnull(x)]
        geneListUp.append(list_no_nan)
    flatUp = [item for sublist in geneListUp for item in sublist]
    flatUP1sansNA = list(dict.fromkeys(flatUp))
    flatDown1sansNA = list(dict.fromkeys(flatDown))
    versListDict[versus] = [flatUP1sansNA, flatDown1sansNA]
print("Versuses selectedAAA: ", versListDict.keys())
print("Versuses selectedBBB: ", versListDict)
######################################################################################################
#### Printing lists
######################################################################################################
for versus in versListDict:
    print("Versus: ", versus)
    print("UP: ", versListDict[versus][0])
    print("DOWN: ", versListDict[versus][1])

if lists == True or All == True:
    for key in versListDict:
        filename = f"{key.upper()}-UPandnDOWN.xlsx"
        dfU = pd.DataFrame(versListDict[key][0]) # to load less packages use pandas
        dfD = pd.DataFrame(versListDict[key][1]) # to load less packages use pandas
        with pd.ExcelWriter(f"{outputFolder}/{filename}", engine='xlsxwriter') as writer:  
            dfU.to_excel(writer, sheet_name='UP', index=False, header=False)
            dfD.to_excel(writer, sheet_name='DOWN', index=False, header=False)
        print("Versus-files have been written")

######################################################################################################
#### Searching for wanted GO terms
######################################################################################################
import GoSearch as gs

if GOterm == True or All == True:
    inputWantedGOFile = inputFolder + "/wantedGo.xlsx"
    #Any clusters looked for?
    GoClusterDictFile = inputFolder + "/GOGroupings.xlsx"
    GoClusterDict = gs.IdentifyClusteredGO(GoClusterDictFile)
    #User provided list of GO terms they want to search for
    ExtGoListFile = inputFolder + "/GoList.txt" 
    #this creates a list with the GO terms in the file. Some might be single GO terms, some might be clusters
    ExtGoList = readExternalList(ExtGoListFile)
    versGODict = {}
    for key in versPDDict.items():
        versus = key[0]
        up2 = key[1][0]
        down2 = key[1][1]
        listOfGoTerms = gs.goSearch(inputWantedGOFile, down2, up2)
        #Create a dict from the GO terms in ExtGoList with up and down genes lists as values
        GODict = {}
        for item in ExtGoList:
            ExtGoTup = gs.GoTermUpDownListMaker(item, GoClusterDict, down2, up2) #this creates a tuple with the up and down genes for the GO term: Up is returned first
            GODict[item]=ExtGoTup
        versGODict[versus] = GODict 
print("VERSUS GO DICT: ", versGODict)
######################################################################################################
#### Sorting for risk genes
######################################################################################################

# #sorting for risk genes using the risk gene list
# riskGenes = pd.read_excel("KnownRiskGenes.xlsx")
# riskList = riskGenes["riskGene"].tolist()
# riskUp = list(set(riskList) & set(flatUPsansNA2))
# print(riskUp)
# riskDown = list(set(riskList) & set(flatDOWNsansNA2))
# print(riskDown)
# dfRisk = pd.DataFrame()

# dfRisk["riskUp"] = riskUp
# dfRisk["riskDown"] = riskDown

# print(dfRisk)
# filenameRisk = f"{versus}_SortedRisk.xlsx"
# dfRisk.to_excel(filenameRisk, index=False)
# print(f"Up-regulated genes for {userinputGO2} in {versus.upper()} that are known risk genes are:\n")

######################################################################################################
#### Printing out  for GO
######################################################################################################


if GOlists == True or All == True:
    dicttt = gs.create_dataframe_dict(versGODict)
    for key, value in dicttt.items():
        pathGOlist = f"{outputFolder}/{key}-GO.xlsx"
        df = pd.DataFrame(value)
        df.to_excel(pathGOlist, index=False)
        print("File has been written")

######################################################################################################
#### Fold change/Foldchange
######################################################################################################
print(f"VERSSDOWN: {versListDict[versus][1]}")

if FCperGOterm == True or FCperVersus == True or All == True:
    from topGeneFold import *
    print("Please wait. Catching fold change..")
    mapDict = getEnsemblMappings()
    print("Mappings have been loaded.")


if FCperVersus == True or All == True: 
    FCversDict = {}
    for vers in versListDict.keys():
        versUpCase=vers.upper()
        foldTuple = topGeneFoldGeneName(versUpCase, mapDict)
        dfUp = foldTuple[1]
        dfDown = foldTuple[2]
        #print("UP: ", dfUp)
        #print("DOWN: ", dfDown)
        #remove rows in dfUp["geneID"] that is not in the list, versListDict[vers][0]
        dfUpSorted = dfUp[dfUp['GeneID'].isin(versListDict[vers][0])]
        dfDownSorted = dfDown[dfDown['GeneID'].isin(versListDict[vers][1])]
        with pd.ExcelWriter(f"{outputFolder}/{vers}-FC.xlsx", engine='xlsxwriter') as writer:  
            dfUpSorted.to_excel(writer, sheet_name='UP', index=False)
            dfDownSorted.to_excel(writer, sheet_name='DOWN', index=False)

if (FCperGOterm == True and GOterm == True) or All == True:
    FCGOtermDict = {}
    for vers in versGODict.keys():
        versUpCase=vers.upper()
        geneFoldTuple = topGeneFoldGeneName(versUpCase, mapDict)
        fcFileNameWOGO = f"{versUpCase}-GO_FC"
        fcFileName1WOGO = fcFileNameWOGO.upper()
        fcFileName2WOG = f"{outputFolder}/{fcFileName1WOGO}.xlsx"
        print(f"Fold change has been calculated\n The dataframes are too big to display here, but the excel file has been saved as {fcFileName2WOG}")
        dfUpWOGO = geneFoldTuple[1]
        dfDownWOGO = geneFoldTuple[2]
        with pd.ExcelWriter(fcFileName2WOG) as writer:
                dfUpWOGO.to_excel(writer, sheet_name='UP', index=False)
                dfDownWOGO.to_excel(writer, sheet_name='DOWN', index=False) 
        #To find regulated genes for the GO terms
        for key in versGODict[vers].keys():
            UpList = versGODict[vers][key][0] 
            DownList = versGODict[vers][key][1]
            dfUpWithGO = dfUpWOGO.query('GeneID in @UpList')
            dfDownWithGO = dfDownWOGO.query('GeneID in @DownList')
            fcFileName = f"{vers}_{key}_FC"
            fcFileNamePath = f"{outputFolder}/{fcFileName}.xlsx"
            with pd.ExcelWriter(fcFileNamePath) as writer:
                    dfUpWithGO.to_excel(writer, sheet_name='UP', index=False)
                    dfDownWithGO.to_excel(writer, sheet_name='DOWN', index=False) 
            UpGenes = extractJustGenes(dfUpWithGO)
            #remove nan values
            UpGenes = [x for x in UpGenes if str(x) != 'nan']
            UpGenes = [i for i in UpGenes if i != '']
            DownGenes = extractJustGenes(dfDownWithGO)
            #remove nan values
            DownGenes = [x for x in DownGenes if str(x) != 'nan']
            DownGenes = [i for i in DownGenes if i != '']
        #print(f"Downgenes are {DownGenes}, upgenes are {UpGenes}")
            # if Perts == True:
            #     import requests
            #     import json
            #     genDictPerts = {}
            #     genDictPerts["upGenes"] = UpGenes
            #     genDictPerts["dnGenes"] = DownGenes
            #     url = 'https://maayanlab.cloud/L1000CDS2/query'
            #     config = {"aggravate":True,"searchMethod":"geneSet","share":False,"combination":True,"db-version":"latest"}
            #     metadata = [{"key":"Tag","value":fcFileName1}]
            #     payload = {"data":genDictPerts,"config":config,"meta":metadata}
            #     headers = {'content-type':'application/json'}
            #     r = requests.post(url,data=json.dumps(payload),headers=headers)
            #     resGeneSet = r.json()
            #     dfFromDict = pd.DataFrame.from_dict(resGeneSet['topMeta'])
            #     fcFileName3 = f"{fcFileName1}_perts.xlsx"
            #     dfFromDict.to_excel(fcFileName3, sheet_name='Sheet1', index=True, header=True, startrow=0, startcol=0, engine='xlsxwriter') 
            #     print(f"The results have been saved as {fcFileName3}")

#####################################################################################################
### Venn
#####################################################################################################
#Create the dictionary for Venn:

# if (Venn == True or All == True) and len(list(goDict)) > 1:
#    import VennModule as vennMod
    






#####################################################################################################
### Union/intersection for up/down regulated genes   
#####################################################################################################

#Opptak 
import UupsideU as uuu

GoListTXT = uuu.readExtGoList("inputFolder/SilicoScreening.txt")
GoDictSilicoScreen = uuu.unclusteredGo("inputFolder/GOGroupings.xlsx")

#GoDict[0] = unclustered GO groups and consists of Opptak as key and a list of GO terms as value
GoDictSilicoScreen0 = GoDictSilicoScreen[0]
CorrectGOdict = removeNotInGoTXT(GoListTXT, GoDictSilicoScreen0) 
#VersGODict = {versus: {GOterm: [upregulated genes, downregulated genes]}}


# Loop over VersGODict and for each key (versus) loop over the values (GOterm) an collect the GO terms that correspond to the key in GOdict and then collect the upregulated and downregulated genes for each GO term.


# SampleDict = {}
# geneDict = {}
# UpList = []
# DownList = []
# for key, val in versGODict.items():
#     #key = versus
#     #val = {GOterm: [upregulated genes, downregulated genes]}
#     #SubDict now
#     for secKey, secVal  in versGODict[key].items():
#         #secKey = GOterm
#         #secVal = [upregulated genes, downregulated genes]
#         for Gokey, Goval in CorrectGOdict.items():
#             #Gokey = GOgroup
#             #Goval = [GOterms]
#             for term in Goval:
#                 if term == secKey:
#                     print(f"term is {term} and secKey is {secKey}")
#                     UpList.append(versGODict[key][secKey][0])
#                     DownList.append(versGODict[key][secKey][1])
#     SampleDict[key] = (UpList, DownList)



SampleDict = {}
geneDict = {}
UpList = []
DownList = []
for Gokey, Goval in CorrectGOdict.items():
    #Gokey = GOgroup
    #Goval = [GOterms]
    for term in Goval:
        for key, val in versGODict.items():
        #key = versus
        #val = {GOterm: [upregulated genes, downregulated genes]}
        #SubDict now
            for secKey, secVal  in versGODict[key].items():
                #secKey = GOterm
                #secVal = [upregulated genes, downregulated genes]
                if term == secKey:
                    UpList.append(versGODict[key][secKey][0])
                    DownList.append(versGODict[key][secKey][1])
                    for i in UpList:
                        if i == []:
                            UpList.remove(i)
                    for i in DownList:
                        if i == []:
                            DownList.remove(i)
                    UpList2 = [item for sublist in UpList for item in sublist]
                    DownList2 = [item for sublist in DownList for item in sublist]
                    UpList3 = [item for i, item in enumerate(UpList2) if item not in UpList2[:i]]
                    DownList3 = [item for i, item in enumerate(DownList2) if item not in DownList2[:i]]
                    print(f"UpList3 = {UpList3}")
                    print(f"DownList3 = {DownList3}")
                #make UpList into a single list
            SampleDict[key] = (UpList3, DownList3)
            print(f"SampleDict = {SampleDict}")
        UpList = []
        DownList = []
    geneDict[Gokey] = SampleDict            

