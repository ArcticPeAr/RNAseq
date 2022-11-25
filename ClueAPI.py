 """
 This program willl be used to pass gene lists to Clues API and return the results in a JSON format
 """
import requests
import json
import pandas as pd
import numpy as np

#create up and dowm gene lists for testing
upL = ["BIRC2",
"NEDD4",
"CNOT4",
"PML",
"SIAH2",
"UBE2D1",
"UBE2D3",
"TRIM25",
"PIAS1",
"UBE2L6",
"TRIM38",
"RBCK1",
"FBXO6",
"KLHL20",
"UBE2S",
"HERC5",
"ANKIB1",
"HERC6",
"RNF114",
"SMURF1",
"HECW2",
"RNF213",
"UBE2Z",
"RNF185",
"RFFL",
"RNF19B",
"DTX3L",
"UBR1",
"RNF144B",
"RNF149",
"RNF175",
"NCCRP1"
]

downL = ["ARRB1",
"CLU",
"HSPA1B",
"PRKN",
"PCBP2",
"RPL5",
"ELOB",
"UBE2G1",
"UCHL1",
"CUL3",
"USP2",
"RPL23",
"BAG5",
"STUB1",
"ERLIN1",
"FAF1",
"ERLIN2",
"UBE2D4",
"TOLLIP",
"DET1",
"YOD1",
"NKD2",
]













#Get API key from file
with open('Clue_API_key.txt', 'r') as f:
    api_key = f.read()
    api_key = api_key.replace('\n', '')
    
        
#Get gene list from file
with open('/home/petear/Downloads/example_uptag_CRCGN009.gmt', 'r') as f:
    gene_list = f.read().splitlines()
    

#set curl command
curl = 'https://api.clue.io/api/jobs'

#change gene list to string separated by tabs
gene_list = '\t'.join(gene_list)

#set new line character at end of gene list
gene_list = gene_list + '\n'

#remove new line character from gene list
gene_list = gene_list.replace('n', '')

#save string to file
with open('ulistGMT.txt', 'w') as f:
    f.write(ulistGMT)

#send curk command to API
r = requests.post(curl, data = {'api_key': api_key, 'ulistGMT': gene_list})

response = requests.post('https://api.clue.io/api/jobs', headers=headers, data=data, auth=(api_key, ''))