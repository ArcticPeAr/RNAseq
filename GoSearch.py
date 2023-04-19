import pandas as pd

def goSearch(wantedGO, down2, up2):
    '''This function takes a list of GO terms and searches for them in the columns of a dataframe. If the GO term is found in the columns of the dataframe, the function will add the column name to a list. The function will then check the values in the columns corresponding to the list and add the values to a list. The function will then remove nan values from the list and return the list. The list is GO terms that are found in the dataframe.'''
    optionsGO = []
    Go = pd.read_excel(wantedGO)
    colListSem = Go[Go.columns[0]].values.tolist()
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
    GoList = list(set(stringsDF))
    GoList.sort()
    return GoList

#Function accepts a dataframe of GO terms and returns a dictionary of GO terms with the column name as the key. If the first value of a column is the name of another column, the function will add the values of that column with the values of the other columns in that column to the dictionary.
#Returns a dictionary of GO terms with the column name as the key

#GOgroups = pd.read_excel("GOGroupings.xlsx", dtype = "str")
def IdentifyClusteredGO(GOgroups):
    '''This function takes a dataframe with GO terms in the columns and a list of GO terms in the first row. If the first row contains a GO term that is also a column name, the function will add the column name to a dictionary as a key and the values in the column as a list. If the first row does not contain a GO term that is also a column name, the function will add the column name to a dictionary as a key and the values in the column as a list. The function will then check the values in the columns corresponding to the lists in the dictionary and add the values to the dictionary as a list. The function will then remove nan values from the lists and return the dictionary.'''
    #check first value of each column contains a string that corresponds to a previous column name. If so, add the first column name to dictionary as key
    groups = pd.read_excel(GOgroups, dtype = "str")
    multiClustDict = {}
    FinalDict = {}
    for col in groups.columns:
        colDF = groups.loc[0, col]
        if colDF in groups.columns:
            colList = groups[col].tolist()
            colList = [x for x in colList if x == x]
            multiClustDict[col] = colList
        else:
            goList = groups[col].tolist()
            ColumnValueList1 = []
            for item in goList:
                ColumnValueList1.append(item)
            FinalDict[col] = ColumnValueList1
    #Get the values in the columns corresponding to the lists in the dictionary
    ColumnValueList = []
    for key in multiClustDict:
        for value in multiClustDict[key]:
            tempList = groups[value].tolist()
            tempList = [x for x in tempList if x == x]
            for item in tempList:
                ColumnValueList.append(item)
            FinalDict[key] = ColumnValueList
    #remove nan values from lists
    for key in FinalDict:
        FinalDict[key] = [x for x in FinalDict[key] if x == x]
    return FinalDict

def readExtGoList(GOtxt):
    '''Reads a text file with one GO term per line and adds those to a list.'''
    GOlistFile = []
    with open(GOtxt, 'r') as f:
        for line in f:
            GOlistFile.append(line.strip().lower())
    return GOlistFile

def ExtTermAvailable(ClusterDict,Golist):
    '''This function takes a list of GO terms and a dictionary of clustered go terms and returns a list of common GO terms and of terms in the list that are not in the dictionary.'''
    commonTerms = []
    for term in Golist:
        if term in ClusterDict:
            commonTerms.append(term)
            uncommonTerms = [x for x in Golist if x not in commonTerms]
    return commonTerms, uncommonTerms



def GoTermUpDownListMaker(GoTerm, ClusterDict, down2, up2):
    '''This function takes a GO term, a dictionary of clustered GO terms, and two dataframes. The function will check if the GO term is in the dictionary. If it is, the function will add the values in the columns corresponding to the GO term to a list. If the GO term is not in the dictionary, the function will add the values in the columns corresponding to the GO term to a list. The function will then remove nan values from the list and return the list.'''
    down3List = []
    up3List = []
    keyLowerList = []
    for key, value in ClusterDict.items():
        key = key.lower()
        keyLowerList.append(key)
    if GoTerm in keyLowerList:
        GoTerm = GoTerm.capitalize()
        for item in ClusterDict[GoTerm]:
            down3 = [col for col in down2.columns if item in col]
            up3 = [col for col in up2.columns if item in col]
            for tem in down3:
                dd = down2[down3].values.tolist()
                down3List.append(dd)
                up3 = [col for col in up2.columns if item in col]
            for stem in up3:
                dd = up2[up3].values.tolist()
                up3List.append(dd)
    else:
        down3 = [col for col in down2.columns if GoTerm in col]
        for tem in down3:
            dd = down2[down3].values.tolist()
            down3List.append(dd)
        up3 = [col for col in up2.columns if GoTerm in col]
        for stem in up3:
            dd = up2[up3].values.tolist()
            up3List.append(dd)
    #make em flat
    flatUp3 = [item4 for sublist in up3List for item4 in sublist]
    flatUp4 = [item4 for sublist in flatUp3 for item4 in sublist]
    flatUPsansNA = list(dict.fromkeys(flatUp4))
    flatUPsansNA2 = [item for item in flatUPsansNA if not(pd.isnull(item)) == True]
    flatDown3 = [item4 for sublist in down3List for item4 in sublist]
    flatDown4 = [item4 for sublist in flatDown3 for item4 in sublist]
    flatDOWNsansNA = list(dict.fromkeys(flatDown4))
    flatDOWNsansNA2 = [item for item in flatDOWNsansNA if not(pd.isnull(item)) == True]
    return flatUPsansNA2, flatDOWNsansNA2

def createDataframeDict(versus_dict):
    '''function takes a dictionary with versus names as keys and a subdict as values. The subdict contais GO terms as keys and a tuple as values. Each tuple contains two lists of gene names. This function will create a new dictionary with the versus names as keys and a pandas dataframe as values. The dataframe will consist of columns with the GO terms + UP or DOWN. If the name is GO term + UP the rowsn will be the first list from the tuple and if its GO term + DOWN the rows will consist of the second list from the tuple'''
    dataframe_dict = {}
    for versus_name, subdict in versus_dict.items():
        dataframe = pd.DataFrame(columns=[f"{go_term}_DOWN" for go_term in subdict.keys()] + [f"{go_term}_UP" for go_term in subdict.keys()])
        for go_term, (down_genes, up_genes) in subdict.items():
            for gene_name in down_genes:
                dataframe.loc[gene_name, f"{go_term}_DOWN"] = gene_name
            for gene_name in up_genes:
                dataframe.loc[gene_name, f"{go_term}_UP"] = gene_name
        dataframe_dict[versus_name] = dataframe.fillna("")
    return dataframe_dict

