#This script reads in Go terms found on maayanlabs with CTS-1027 regulated genes and searches for the terms in a list of GO terms from the GO groupings excel file.
#CHANGED TO PYTHON SCRIPT
import pandas as pd

ctsbpU = pd.read_csv("/home/petear/CTS/CTS_BP_Up.txt", delimiter='\t')
ctsmfU = pd.read_csv("/home/petear/CTS/CTS_MF_Up.txt", delimiter='\t')
ctsccU = pd.read_csv("/home/petear/CTS/CTS_CC_Up.txt", delimiter='\t')

ctsUp = pd.concat([ctsbpU, ctsmfU, ctsccU], axis=0, ignore_index=True)

ctsbpD = pd.read_csv("/home/petear/CTS/CTS_BP_Down.txt", delimiter='\t')
ctsmfD = pd.read_csv("/home/petear/CTS/CTS_MF_Down.txt", delimiter='\t')
ctsccD = pd.read_csv("/home/petear/CTS/CTS_CC_Down.txt", delimiter='\t')

ctsDown = pd.concat([ctsbpD, ctsmfD, ctsccD], axis=0, ignore_index=True)


cluster = pd.read_excel("/home/petear/GOGroupings.xlsx", sheet_name = "Sheet1")


OpptakList=cluster["Opptak"].tolist()
OpptakList = [x for x in OpptakList if str(x) != 'nan']

patternOpptakList = r'\b(?:{})\b'.format('|'.join(OpptakList))

DegraderingList=cluster["Degradering"].tolist()
DegraderingList = [x for x in DegraderingList if str(x) != 'nan']

patternDegraderingList = r'\b(?:{})\b'.format('|'.join(DegraderingList))

ctsUpOpptak = ctsUp[ctsUp['Term'].str.contains(patternOpptakList)]
ctsDownOpptak = ctsDown[ctsDown['Term'].str.contains(patternOpptakList)]



ctsUpDegradering = ctsUp[ctsUp['Term'].str.contains(patternDegraderingList)]
ctsDownDegradering = ctsDown[ctsDown['Term'].str.contains(patternDegraderingList)]

################################################
#Make list of both Opptak and Degradering lists
################################################


listlist = [OpptakList, DegraderingList]

#make list of lists into one list
listlist = [item for sublist in listlist for item in sublist]

#remove NaN
listlist = [x for x in listlist if str(x) != 'nan']


patternListlist = r'\b(?:{})\b'.format('|'.join(listlist))


ctsUpFiltered = ctsUp[ctsUp['Term'].str.contains(patternListlist)]
ctsDownFiltered = ctsDown[ctsDown['Term'].str.contains(patternListlist)]


with pd.ExcelWriter("CTS.xlsx") as writer:
    ctsUpOpptak.to_excel(writer, sheet_name='UpOpptak', index=False)
    ctsDownOpptak.to_excel(writer, sheet_name='DownOpptak', index=False)
    ctsUpDegradering.to_excel(writer, sheet_name='UpDegradering', index=False)
    ctsDownDegradering.to_excel(writer, sheet_name='DownDegradering', index=False)
    ctsUpFiltered.to_excel(writer, sheet_name='UpDegradOpptak', index=False)
    ctsDownFiltered.to_excel(writer, sheet_name='DownDegradOpptak', index=False)