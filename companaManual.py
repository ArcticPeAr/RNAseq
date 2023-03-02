 ####Library imports
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns



######################################################################################################
####Premade GO groups from excel file and make dictionary of each column with column name as key   
######################################################################################################

#Function accepts a dataframe of GO terms and returns a dictionary of GO terms with the column name as the key. If the first value of a column is the name of another column, the function will add the values of that column with the values of the other columns in that column to the dictionary.
#Returns a dictionary of GO terms with the column name as the key

def IdentifyClusteredGO(GOgroups):
    #check first value of each column contains a string that corresponds to a previous column name. If so, add the first column name to dictionary as key
    multiClustDict = {}
    FinalDict = {}
    for col in GOgroups.columns:
        colDF = GOgroups.loc[0, col]
        if colDF in GOgroups.columns:
            colList = GOgroups[col].tolist()
            colList = [x for x in colList if x == x]
            multiClustDict[col] = colList
        else:
            goList = GOgroups[col].tolist()
            ColumnValueList1 = []
            for item in goList:
                ColumnValueList1.append(item)
            FinalDict[col] = ColumnValueList1
    #Get the values in the columns corresponding to the lists in the dictionary
    ColumnValueList = []
    for key in multiClustDict:
        for value in multiClustDict[key]:
            tempList = GOgroups[value].tolist()
            tempList = [x for x in tempList if x == x]
            for item in tempList:
                ColumnValueList.append(item)
            FinalDict[key] = ColumnValueList
    #remove nan values from lists
    for key in FinalDict:
        FinalDict[key] = [x for x in FinalDict[key] if x == x]
    return FinalDict





######################################################################################################
####Find samples names DOWN:
######################################################################################################

bdd = pd.read_csv("BeDF_DOWN_perGO.csv", dtype = "str")
bdd.columns = bdd.columns.str.lower()
cdd = pd.read_csv("CeDF_DOWN_perGO.csv", dtype = "str")
cdd.columns = cdd.columns.str.lower()
mdd = pd.read_csv("MeDF_DOWN_perGO.csv", dtype = "str")
mdd.columns = mdd.columns.str.lower()

framesDown = [bdd, cdd, mdd]
down = pd.concat(framesDown   )
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

bdu = pd.read_csv("BeDF_UP_perGO.csv", dtype = "str")
bdu.columns = bdu.columns.str.lower()
cdu = pd.read_csv("CeDF_UP_perGO.csv", dtype = "str")
cdu.columns = cdu.columns.str.lower()
mdu = pd.read_csv("MeDF_UP_perGO.csv", dtype = "str")
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

optionsVerDup = optionsVerUp + optionsVerDown
optionsVer = [val for val in optionsVerDup if optionsVerDup.count(val) > 1]
optionsVer = list(set(optionsVer))

#
######################################################################################################

# bda = pd.read_csv("BeDF_FULL_perGO.csv")
# cda = pd.read_csv("CeDF_FULL_perGO.csv")
# mda = pd.read_csv("MeDF_FULL_perGO.csv")

# framesAll = [bdd, cdd, mdd]
# down = pd.concat(framesAll)

######################################################################################################
#### Starting input
######################################################################################################
import re
def add_underscore(s):
    s = re.sub(r'([a-zA-Z])([0-9])', r'\1_\2', s)
    return s


print("Your available samples are:")
for item in optionsVer:
    item = item.upper()
    print(item)

usVerInput1 = input('Please select the first group to compare. E.g "C_1 or T_12"\n')
usVerInput1 = usVerInput1.lower()
usVerInput1 = add_underscore(usVerInput1)
while usVerInput1 not in optionsVer:
    if len(usVerInput1) > 4:
        splitted = usVerInput1.split("-")
        if splitted[-1] in optionsVer:
            usVerInput2 = splitted[-1]
            usVerInput2 = usVerInput2.lower()
            usVerInput1 = splitted[0]
            usVerInput1 = usVerInput1.lower()
        else:
            usVerInput1 = input('Error: Correct syntax is: "C_1" or "C_1-T_8"\n')
            usVerInput1 = usVerInput1.lower()
            usVerInput1 = add_underscore(usVerInput1)
    else:
        usVerInput1 = input('Error: Correct syntax is: "C_1" or "C_1-T_8"\n')
        usVerInput1 = usVerInput1.lower()
        usVerInput1 = add_underscore(usVerInput1)

if 'usVerInput2' not in locals():
    usVerInput2 = input('Please select the second group to compare. E.g "C_2 or T_11"\n')
    usVerInput2 = usVerInput2.lower()
    usVerInput2 = add_underscore(usVerInput2)
    while usVerInput2 not in optionsVer or usVerInput1 == usVerInput2:
        usVerInput2 = input('Error: remember: "T_2", not "T2 and first sample cannot be the same as second sample"\n')
        usVerInput2 = usVerInput2.lower()
        usVerInput2 = add_underscore(usVerInput2)

    
# while usVerInput2 == usVerInput1:
#     usVerInput2 = input('Error: Second sample cannot be the same as first sample!\nPlease type new sample: \n')

print(usVerInput1.upper())
print(usVerInput2.upper())
usVerStringF = f"{usVerInput1}-{usVerInput2}"
usVerStringR = f"{usVerInput2}-{usVerInput1}"

if usVerStringF in StrList:
    print(f"The versus is this: {usVerStringF.upper()}")
    print("Computing genes for this versus")
    versus = usVerStringF
elif usVerStringR in StrList:
    print(f"The versus is this: {usVerStringR.upper()}")
    print("Computing genes for this versus")
    versus = usVerStringR
else:
    print("VERSUS IS NOT AVAILABLE!")
    print("Run program again with another comparison")
    quit()

down2 = down.filter(regex=versus)
up2 = up.filter(regex=versus)


geneListUp = list()
geneListDown = list()

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

######################################################################################################
#### Printing lists
######################################################################################################

print(f"Up-regulated genes for ALL GO terms are:: {flatUp}\n")
print(f"Down-regulated genes for ALL GO terms are:: {flatDown}\n")
versUpper = versus.upper()
print(f"The comparison setup was this: {versUpper}\n")

print("Do you want to save these lists")
userinputXL = input('Please enter yes or no":\n')
userinputXL = userinputXL.lower()

optionsXL =["yes","y", "no","n"]

while userinputXL not in optionsXL:
    print('Error, please type "yes", "y", "no" or "n"')
    userinputXL = input('Please enter yes or no":\n')
    userinputXL = userinputXL.lower()

filename = f"{versUpper}_UPandDOWN_GO.xlsx"

if userinputXL == "yes" or userinputXL == "y":
    dfU = pd.DataFrame(flatUp) # to load less packages use pandas
    dfD = pd.DataFrame(flatDown) # to load less packages use pandas
    with pd.ExcelWriter(filename) as writer:  
        dfU.to_excel(writer, sheet_name='UP', index=False, header=False)
        dfD.to_excel(writer, sheet_name='DOWN', index=False, header=False)
    print("File has been written")


######################################################################################################
#### Searching for GO
######################################################################################################
#is GO in either of the dfs?
optionsGO = []

GO = pd.read_excel("wantedGo.xlsx")
colListSem = GO[GO.columns[0]].values.tolist()
optionsGO.append(colListSem)



flatGO = [item for sublist in optionsGO for item in sublist]

stringsInDFColDOWN = []

for iteo in flatGO:
    iteo = iteo.lower()
    aa = down2.columns.str.contains(iteo)
    ITEO = iteo.upper()
    AA = down2.columns.str.contains(ITEO)
    if any(aa) == True or any(AA) == True:
        stringsInDFColDOWN.append(iteo)

stringsInDFColUP = []

for ites in flatGO:
    ites = ites.lower()
    aa = up2.columns.str.contains(ites)
    ITES = ites.upper()
    AA = up2.columns.str.contains(ITES)
    if any(aa) == True or any(AA) == True:
        stringsInDFColUP.append(ites)


stringsInDFCol = stringsInDFColUP + stringsInDFColDOWN
stringsDF = [vals for vals in stringsInDFCol if stringsInDFCol.count(vals) > 1]
stringsDF = list(set(stringsDF))
stringsDF.sort()

######################################################################################################
#### Clustered GO terms
######################################################################################################
GoClusters = pd.read_excel("GOGroupings.xlsx", dtype = "str")
ClusterDict = IdentifyClusteredGO(GoClusters)

        



######################################################################################################
#### Go term loop:
######################################################################################################

goDict ={}

userInputReGO = "g"
while userInputReGO == "g":
    print(f"The terms you have available for the versus {versus.upper()} are:\n")

    for itep in stringsDF:
        print(itep)

    
    print("\n")
    print("The terms you have available as premade clusters are:")
    keyLowerList = []
    for key, value in ClusterDict.items():
        print(f"{key.upper()}")
        key = key.lower()
        keyLowerList.append(key)

    stringUserLooksfor = stringsDF + keyLowerList

    userinputGO2 = input('\nPlease enter GO term or GO-cluster you are looking for:\n')
    userinputGO2 = userinputGO2.lower()
    
    while userinputGO2 not in stringUserLooksfor:
        userinputGO2 = input('Are you sure you typed the term correctly?\n hint: You can just copy paste from above:\n')
        userinputGO2 = userinputGO2.lower()


    #Make list from userinput and columns matching that term
    down3List = []
    up3List = []

    if userinputGO2 in keyLowerList:
        userinputGO2 = userinputGO2.capitalize()
        for item in ClusterDict[userinputGO2]:
            down3 = [col for col in down2.columns if item in col]
            for tem in down3:
                dd = down2[down3].values.tolist()
                down3List.append(dd)

            up3 = [col for col in up2.columns if item in col]
            for stem in up3:
                dd = up2[up3].values.tolist()
                up3List.append(dd)

    else:
        down3 = [col for col in down2.columns if userinputGO2 in col]
        for tem in down3:
            dd = down2[down3].values.tolist()
            down3List.append(dd)
        
        up3 = [col for col in up2.columns if userinputGO2 in col]
        for stem in up3:
            dd = up2[up3].values.tolist()
            up3List.append(dd)


    flatUp3 = [item4 for sublist in up3List for item4 in sublist]
    flatUp4 = [item4 for sublist in flatUp3 for item4 in sublist]
    
    flatUPsansNA = list(dict.fromkeys(flatUp4))
    flatUPsansNA2 = [item for item in flatUPsansNA if not(pd.isnull(item)) == True]


    flatDown3 = [item4 for sublist in down3List for item4 in sublist]
    flatDown4 = [item4 for sublist in flatDown3 for item4 in sublist]
    flatDOWNsansNA = list(dict.fromkeys(flatDown4))
    flatDOWNsansNA2 = [item for item in flatDOWNsansNA if not(pd.isnull(item)) == True]

    print(f"Up-regulated genes for {userinputGO2} in {versus.upper()} are:\n")
    print(f"{flatUPsansNA2} \n")
    print(f"Down-regulated genes for {userinputGO2} in {versus.upper()} are:\n")
    print(f"{flatDOWNsansNA2} \n")


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


    print("Do you want to save these lists")
    userinputXL2 = input('Please enter yes or no":\n')
    userinputXL2 = userinputXL2.lower()


    while userinputXL2 not in optionsXL:
        print('Error, please type "yes", "y", "no" or "n"')
        userinputXL2 = input('Please enter yes or no":\n')
        userinputXL2 = userinputXL2.lower()

    filename2 = f"{versus}_UPandDOWN_{userinputGO2}.xlsx"

    if userinputXL2 == "yes" or userinputXL2 == "y":
        dfU2 = pd.DataFrame(flatUPsansNA2) # to load less packages use pandas
        dfD2 = pd.DataFrame(flatDOWNsansNA2) # to load less packages use pandas
        with pd.ExcelWriter(filename2) as writer:  
            dfU2.to_excel(writer, sheet_name='UP', index=False, header=False)
            dfD2.to_excel(writer, sheet_name='DOWN', index=False, header=False)
        print("File has been written")
        print("Thanks for using this little program!")
    elif userinputXL2 == "no" or userinputXL2 == "n":
        print("Thanks for using this little program!")
        
    ######################################################################################################
    #### Fold change/Foldchange
    ######################################################################################################
    
    print("Do you want to calculate the fold change for these lists\nWARNING this will take a while")
    
    userinputFC = input('Please enter y or yes. Anything else will abort this operation:\n')
    userinputFC = userinputFC.lower()
    versus = versus.upper()

    if userinputFC == "yes" or userinputFC == "y":
        from topGeneFold import *
        print("Please wait. Catching fold change..")
        mapDict = getEnsemblMappings()
        print("Mappings have been loaded.")
        geneFoldTuple = topGeneFoldGeneName(versus, mapDict)
        fcFileName = f"{versus}_FC_{userinputGO2}"
        fcFileName1 = fcFileName.upper()
        fcFileName2 = f"{fcFileName1}.xlsx"
        print(f"Fold change has been calculated\n The dataframes are too big to display here, but the excel file has been saved as {fcFileName}")
        dfUp = geneFoldTuple[1]
        dfDown = geneFoldTuple[2]
        dfUp = dfUp.query('GeneID in @flatUPsansNA2')
        dfDown = dfDown.query('GeneID in @flatDOWNsansNA2')
        with pd.ExcelWriter(fcFileName2) as writer:
                dfUp.to_excel(writer, sheet_name='UP', index=False)
                dfDown.to_excel(writer, sheet_name='DOWN', index=False) 
        UpGenes = extractJustGenes(dfUp)
        #remove nan values
        UpGenes = [x for x in UpGenes if str(x) != 'nan']
        UpGenes = [i for i in UpGenes if i != '']
        DownGenes = extractJustGenes(dfDown)
        #remove nan values
        DownGenes = [x for x in DownGenes if str(x) != 'nan']
        DownGenes = [i for i in DownGenes if i != '']
        userinputGetPerts = input('Do you want to send these genes to MaayanLabs L1000CDS²?\n WARNING this is an external resource and should not be used indiscriminately\n y or yes to send, anything else to abort:\n')
        userinputGetPerts = userinputGetPerts.lower()
        
        if userinputGetPerts == "yes" or userinputGetPerts == "y":
            import requests
            import json
            genDictPerts = {}
            genDictPerts["upGenes"] = UpGenes
            genDictPerts["dnGenes"] = DownGenes
            url = 'https://maayanlab.cloud/L1000CDS2/query'
            config = {"aggravate":True,"searchMethod":"geneSet","share":False,"combination":True,"db-version":"latest"}
            metadata = [{"key":"Tag","value":fcFileName1}]
            payload = {"data":genDictPerts,"config":config,"meta":metadata}
            headers = {'content-type':'application/json'}
            r = requests.post(url,data=json.dumps(payload),headers=headers)
            resGeneSet = r.json()
            dfFromDict = pd.DataFrame.from_dict(resGeneSet['topMeta'])
            fcFileName3 = f"{fcFileName1}_perts.xlsx"
            dfFromDict.to_excel(fcFileName3, sheet_name='Sheet1', index=True, header=True, startrow=0, startcol=0, engine='xlsxwriter') 
            print(f"The results have been saved as {fcFileName3}")
    
    
    ######################################################################################################
    #### Venn
    ######################################################################################################
    #Create the dictionary for Venn:
    listLists =[]
    listLists.append(flatUPsansNA2)
    listLists.append(flatDOWNsansNA2)
    goDict[userinputGO2] = listLists

    vennList1 = list(goDict)[-1]
    if len(list(goDict)) > 1:
        vennList2 = list(goDict)[-2]

    print(f"{userinputGO2} and its up and down lists have been saved for later use")

#   #Plotting venn and getting 
    if len(goDict) > 1:
        userinputVenn = input('Press "v" you want to make a venn diagram of the following terms:\nPress anything else to forfeit plotting:\n')
        userinputVenn = userinputVenn.lower()
        if userinputVenn == "v":
            print(f"the venn diagram of {vennList1} and {vennList2} is being made")
            #Find the intersection of the two lists for UP:
            intersectionUP = list(set(goDict[vennList1][0]) & set(goDict[vennList2][0]))
            #find genes in first list but not in second for UP:
            firstNotSecondUP = list(set(goDict[vennList1][0]) - set(goDict[vennList2][0]))
            #find genes in second list but not in first for UP:
            secondNotFirstUP = list(set(goDict[vennList2][0]) - set(goDict[vennList1][0]))
            plt.figure(figsize=(15,5))
            vdUP = venn2(subsets = (len(firstNotSecondUP), len(secondNotFirstUP), len(intersectionUP)), set_labels = (vennList1, vennList2))
            titleUp = f'Up for "{vennList1}" and "{vennList2} in {versus}"'
            vdUP = plt.title(titleUp)
            vdUP = plt.savefig(titleUp)
            plt.close(vdUP)
            #Find the intersection of the two lists for DOWN:
            intersectionDOWN = list(set(goDict[vennList1][1]) & set(goDict[vennList2][1]))
            #find genes in first list but not in second for DOWN:
            firstNotSecondDOWN = list(set(goDict[vennList1][1]) - set(goDict[vennList2][1]))
            #find genes in second list but not in first for DOWN:
            secondNotFirstDOWN = list(set(goDict[vennList2][1]) - set(goDict[vennList1][1]))
            plt.figure(figsize=(15,5))
            vdDOWN = venn2(subsets = (len(firstNotSecondDOWN), len(secondNotFirstDOWN), len(intersectionDOWN)), set_labels = (vennList1, vennList2))
            titleDown = f'Down for "{vennList1}" and "{vennList2} in {versus}"'
            vdDOWN = plt.title(titleDown)
            vdDOWN = plt.savefig(titleDown)
            plt.close(vdDOWN)
            #Now printing excel files.
            filenameVenn = f"{versus}_{vennList1}-{vennList2}.xlsx"
            intersectionUPDF = pd.DataFrame(intersectionUP, columns = ["Shared "])
            firstNotSecondUPDF = pd.DataFrame(firstNotSecondUP, columns = ["Only_in " + vennList1])
            secondNotFirstUPDF = pd.DataFrame(secondNotFirstUP, columns = ["Only_in " + vennList2])
            intersectionDOWNDF = pd.DataFrame(intersectionDOWN, columns = ["Shared"])
            firstNotSecondDOWNDF = pd.DataFrame(firstNotSecondDOWN, columns = ["Only_in " + vennList1])
            secondNotFirstDOWNDF = pd.DataFrame(secondNotFirstDOWN, columns = ["Only_in " + vennList2])
            combinedDF_UP = pd.concat([intersectionUPDF, firstNotSecondUPDF, secondNotFirstUPDF], axis=1)
            combinedDF_DOWN = pd.concat([intersectionDOWNDF, firstNotSecondDOWNDF, secondNotFirstDOWNDF], axis=1)
            with pd.ExcelWriter(filenameVenn) as writer:
                combinedDF_UP.to_excel(writer, sheet_name='UP', index=False)
                combinedDF_DOWN.to_excel(writer, sheet_name='DOWN', index=False)  
        else:
            pass

    

    
    userInputReGO = input('Press "g" you want to search for another GO term.\n')
    userInputReGO = userInputReGO.lower()
    
print( "Thanks for using this little program!")
print("Producing bar plot of number of up and down genes belonging for the GO terms you have searched for")



#plot horizontal bar plot of the number of up and down genes for each GO term with wrapped labels with gene count on top of bar
#make a list of the number of up and down genes for each GO term
upList = []
downList = []
for key in goDict:
    upList.append(len(goDict[key][0]))
    downList.append(len(goDict[key][1]))
#make a list of the GO terms
goList = list(goDict)
#make a dataframe of the lists
df = pd.DataFrame(list(zip(goList, upList, downList)), columns =['GO term', 'Up', 'Down'])
#make a bar plot
plt.figure(figsize=(18,10))
ax = sns.barplot(y="GO term", x="value", hue="variable", data=pd.melt(df, ['GO term']))
ax.set_yticklabels(ax.get_yticklabels(), rotation=70, ha="right")
plt.xlabel('Gene count')
plt.ylabel('GO term')
plt.title(f'Number of up and down genes for each GO term searched for in {versus}')
plt.tight_layout()
#add gene count on top of bar
for p in ax.patches:
    ax.annotate(str(p.get_width()), (p.get_width(), p.get_y()+0.55*p.get_height()))
plt.show()

