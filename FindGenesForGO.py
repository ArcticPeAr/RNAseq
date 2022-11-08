import pandas as pd
import numpy as np 

userinputDF = input('Please enter "up", "down", or "all":\n')
optionsDF = ["up", "down", "all"]

while userinputDF not in optionsDF:
    print('Error, remember all are in small letters: "up", "down", or "all"')
    userinputDF = input('Please enter "up", "down", or "all":\n')


if userinputDF == "up":
    bd = pd.read_csv("BeDF_UP_perGO.csv")
    cd = pd.read_csv("CeDF_UP_perGO.csv")
    md = pd.read_csv("MeDF_UP_perGO.csv")
if userinputDF == "down":
    bd = pd.read_csv("BeDF_DOWN_perGO.csv")
    cd = pd.read_csv("CeDF_DOWN_perGO.csv")
    md = pd.read_csv("MeDF_DOWN_perGO.csv")
if userinputDF == "all":
    bd = pd.read_csv("BeDF_FULL_perGO.csv")
    cd = pd.read_csv("CeDF_FULL_perGO.csv")
    md = pd.read_csv("MeDF_FULL_perGO.csv")

optionsGO = []

GO = pd.read_excel("wantedGo.xlsx")
colList = GO[GO.columns[0]].values.tolist()
optionsGO.append(colList)
colList = GO[GO.columns[1]].values.tolist()
optionsGO.append(colList)

flatGO = [item for sublist in optionsGO for item in sublist]

print('The GO terms you can choose from are the following:\n')
print("\n")
print(flatGO)

userinputGO = input('Please enter GO term you are looking for:\n')

while userinputGO not in flatGO:
    print('Error, remember: "endocytosis", not "Endocytosis"\n Are you sure you typed the term correctly?\n hint: You can just copy past from above by ctrl+shift c and ctrl+shift v')
    userinputGO = input('Please enter GO term you are looking for:\n')



frames = [bd, cd, md]
down = pd.concat(frames)
endocytosis = down.columns[down.columns.str.contains(userinputGO)]
pinolist = list(endocytosis)
tt = down[pinolist]

ttlist = list()

for item in pinolist:
    da = tt[item].tolist()
    list_no_nan = [x for x in da if pd.notnull(x)]
    ttlist.append(list_no_nan)


flat_list = [item for sublist in ttlist for item in sublist]

mylist = list(dict.fromkeys(flat_list))

print("\n",mylist,"\n")
