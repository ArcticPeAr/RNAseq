import pandas as pd



bdd = pd.read_csv("BeDF_DOWN_perGO.csv")
cdd = pd.read_csv("CeDF_DOWN_perGO.csv")
mdd = pd.read_csv("MeDF_DOWN_perGO.csv")

######################################################################################################
####Find samples names DOWN:
######################################################################################################

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

bdu = pd.read_csv("BeDF_UP_perGO.csv")
cdu = pd.read_csv("CeDF_UP_perGO.csv")
mdu = pd.read_csv("MeDF_UP_perGO.csv")

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

print("Your available samples are:")
for item in optionsVer:
    print(item)

usVerInput1 = input('Please select first group to compare. E.g "C_1 or T_12"\n')
usVerInput1 = usVerInput1.upper()
while usVerInput1 not in optionsVer:
    usVerInput1 = input('Error: remember: "C_1", not "C1"\n')
    


usVerInput2 = input('Please select first group to compare. E.g "C_2 or T_11"\n')
usVerInput2 = usVerInput2.upper()
while usVerInput2 not in optionsVer or usVerInput1 == usVerInput2:
    usVerInput2 = input('Error: remember: "T_2", not "T2 and first sample cannot be the same as second sample"\n')
    
# while usVerInput2 == usVerInput1:
#     usVerInput2 = input('Error: Second sample cannot be the same as first sample!\nPlease type new sample: \n')


usVerStringF = f"{usVerInput1}-{usVerInput2}"
usVerStringR = f"{usVerInput2}-{usVerInput1}"

if usVerStringF in StrList:
    print(f"The versus is this: {usVerStringF}")
    print("Computing genes for this versus")
    versus = usVerStringF
elif usVerStringR in StrList:
    print(f"The versus is this: {usVerStringR}")
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
for (columnName, columnData) in down2.iteritems():
    da = columnData.tolist()
    list_no_nan = [x for x in da if pd.notnull(x)]
    geneListDown.append(list_no_nan)

flatDown = [item for sublist in geneListDown for item in sublist]

#UP
for (columnName, columnData) in up2.iteritems():
    da = columnData.tolist()
    list_no_nan = [x for x in da if pd.notnull(x)]
    geneListUp.append(list_no_nan)

flatUp = [item for sublist in geneListUp for item in sublist]

######################################################################################################
#### Printing lists
######################################################################################################

print(f"Up-regulated genes for ALL GO terms are:: {flatUp}\n")
print(f"Down-regulated genes for ALL GO terms are:: {flatDown}\n")
print(f"The comparison setup was this: {versus}\n")

print("Do you want to save these lists")
userinputXL = input('Please enter yes or no":\n')
userinputXL = userinputXL.lower()

optionsXL =["yes","y", "no","n"]

while userinputXL not in optionsXL:
    print('Error, please type "yes", "y", "no" or "n"')
    userinputXL = input('Please enter yes or no":\n')
    userinputXL = userinputXL.lower()

filename = f"{versus}_UPandDOWN_GO.xlsx"

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


print("\n")
print("Do you want to continue to filter for GO terms")
userinputGO = input('Please enter yes or no":\n')

while userinputGO not in optionsXL:
    print('Error, please type "yes", "y", "no" or "n"')
    userinputGO = input('Please enter yes or no":\n')
    userinputGO = userinputGO.lower()

if userinputGO == "no" or userinputGO =="n":
    quit()


#is GO in either of the dfs?
optionsGO = []

GO = pd.read_excel("wantedGo.xlsx")
colListSem = GO[GO.columns[0]].values.tolist()
optionsGO.append(colListSem)
colListID = GO[GO.columns[1]].values.tolist()
optionsGO.append(colListID)


flatGO = [item for sublist in optionsGO for item in sublist]

stringsInDFColDOWN = []

for iteo in flatGO:
    aa = down2.columns.str.contains(iteo)
    ITEO = iteo.upper()
    AA = down2.columns.str.contains(ITEO)
    if any(aa) == True or any(AA) == True:
        stringsInDFColDOWN.append(iteo)

stringsInDFColUP = []

for ites in flatGO:
    aa = up2.columns.str.contains(ites)
    ITES = ites.upper()
    AA = up2.columns.str.contains(ITES)
    if any(aa) == True or any(AA) == True:
        stringsInDFColUP.append(ites)


stringsInDFCol = stringsInDFColUP + stringsInDFColDOWN
stringsDF = [vals for vals in stringsInDFCol if stringsInDFCol.count(vals) > 1]
stringsDF = list(set(stringsDF))


print(f"The terms you have available for the versus {versus} are:\n")

for itep in stringsDF:
    print(itep)

print("\n")

userinputGO2 = input('Please enter GO term you are looking for:\n')
while userinputGO2 not in stringsDF:
    userinputGO2 = input('Error, remember: "endocytosis", not "Endocytosis"\n Are you sure you typed the term correctly?\n hint: You can just copy paste from above:\n')


#Make list from userinput and columns matching that term
down3List = []
up3List = []

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



flatDown3 = [item4 for sublist in down3List for item4 in sublist]
flatDown4 = [item4 for sublist in flatDown3 for item4 in sublist]
flatDOWNsansNA = list(dict.fromkeys(flatDown4))

print(f"Up-regulated genes for {userinputGO2} are:\n")
print(f"{flatUPsansNA} \n")
print(f"Down-regulated genes for {userinputGO2} are:\n")
print(f"{flatDOWNsansNA} \n")


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
    dfU2 = pd.DataFrame(flatUPsansNA) # to load less packages use pandas
    dfD2 = pd.DataFrame(flatDOWNsansNA) # to load less packages use pandas
    with pd.ExcelWriter(filename2) as writer:  
        dfU2.to_excel(writer, sheet_name='UP', index=False, header=False)
        dfD2.to_excel(writer, sheet_name='DOWN', index=False, header=False)
    print("File has been written")
    print("Thanks for using this little program!")
elif userinputXL2 == "no" or userinputXL2 == "n":
    print("Thanks for using this little program!")
    quit()