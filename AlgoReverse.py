import pandas as pd

#1.1.2
df1UP_C1T5 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_1-T_5_FC_OPPTAK.xlsx", sheet_name='UP')
df2UP_C7T11 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_7-T_11_FC_OPPTAK.xlsx", sheet_name='UP')

dfUPMerged1_3 = pd.concat([df1UP_C1T5, df2UP_C7T11], axis=0)
#dfUPMerged1_3 = dfUPMerged1_3.drop_duplicates(subset=['GeneID'], keep='first')
dfUPMerged1_3 = dfUPMerged1_3.sort_values(by=['logFC'], ascending=False)
#Reverse logFC values, logCPM values and sampleValues
dfUPMerged1_3['logFC'] = dfUPMerged1_3['logFC'] * -1
dfUPMerged1_3['logCPM'] = dfUPMerged1_3['logCPM'] * -1
dfUPMerged1_3['Versus'] = dfUPMerged1_3['Versus'].apply(lambda x: '-'.join(x.split('-')[::-1]))
#Since the logFC values are reversed, what was considered upregulated for Cx-Tx is now downregulated and vice versa:
DOWNMerged1_3 = dfUPMerged1_3

#1.1.1
df1DOWN_C1T5 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_1-T_5_FC_OPPTAK.xlsx", sheet_name='DOWN')
df2DOWN_C7T11 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_7-T_11_FC_OPPTAK.xlsx", sheet_name='DOWN')

dfDOWNMerged1_3 = pd.concat([df1DOWN_C1T5, df2DOWN_C7T11], axis=0)
#dfDOWNMerged1_3 = dfDOWNMerged1_3.drop_duplicates(subset=['GeneID'], keep='first')
dfDOWNMerged1_3 = dfDOWNMerged1_3.sort_values(by=['logFC'], ascending=True)
#Reverse logFC values, logCPM values and sampleValues
dfDOWNMerged1_3['logFC'] = dfDOWNMerged1_3['logFC'] * -1
dfDOWNMerged1_3['logCPM'] = dfDOWNMerged1_3['logCPM'] * -1
dfDOWNMerged1_3['Versus'] = dfDOWNMerged1_3['Versus'].apply(lambda x: '-'.join(x.split('-')[::-1]))
#Since the logFC values are reversed, what was considered upregulated for Cx-Tx is now downregulated and vice versa:
UPMerged1_3 = dfDOWNMerged1_3

#1.3
UnionDF1_3 = pd.concat([DOWNMerged1_3, UPMerged1_3], axis=0)
#UnionDF1_3 = UnionDF1_3.drop_duplicates(subset=['GeneID'], keep='first')
UnionDF1_3 = UnionDF1_3.sort_values(by=['logFC'], ascending=False)

#IntersectDF1_3 = pd.merge(dfUPMerged, dfDOWNMerged, on='GeneID', how='inner')


############################################################################################################
#1.2.2
df1UP_C1T2 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_1-T_2_FC_OPPTAK.xlsx", sheet_name='UP')
df2UP_C7T8 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_7-T_8_FC_OPPTAK.xlsx", sheet_name='UP')


dfUPMerged1_4 = pd.concat([df1UP_C1T2, df2UP_C7T8], axis=0)
#dfUPMerged1_4 = dfUPMerged1_4.drop_duplicates(subset=['GeneID'], keep='first')
dfUPMerged1_4 = dfUPMerged1_4.sort_values(by=['logFC'], ascending=False)
#Reverse logFC values, logCPM values and sampleValues
dfUPMerged1_4['logFC'] = dfUPMerged1_4['logFC'] * -1
dfUPMerged1_4['logCPM'] = dfUPMerged1_4['logCPM'] * -1
dfUPMerged1_4['Versus'] = dfUPMerged1_4['Versus'].apply(lambda x: '-'.join(x.split('-')[::-1]))
#Since the logFC values are reversed, what was considered upregulated for Cx-Tx is now downregulated and vice versa:
DOWNMerged1_4 = dfUPMerged1_4

#1.2.1
df1DOWN_C1T2 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_1-T_2_FC_OPPTAK.xlsx", sheet_name='DOWN')
df2DOWN_C7T8 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_7-T_8_FC_OPPTAK.xlsx", sheet_name='DOWN')

dfDOWNMerged1_4 = pd.concat([df1DOWN_C1T2, df2DOWN_C7T8], axis=0)
#dfDOWNMerged1_4 = dfDOWNMerged1_4.drop_duplicates(subset=['GeneID'], keep='first')
dfDOWNMerged1_4 = dfDOWNMerged1_4.sort_values(by=['logFC'], ascending=True)
#Reverse logFC values, logCPM values and sampleValues
dfDOWNMerged1_4['logFC'] = dfDOWNMerged1_4['logFC'] * -1
dfDOWNMerged1_4['logCPM'] = dfDOWNMerged1_4['logCPM'] * -1
dfDOWNMerged1_4['Versus'] = dfDOWNMerged1_4['Versus'].apply(lambda x: '-'.join(x.split('-')[::-1]))
#Since the logFC values are reversed, what was considered upregulated for Cx-Tx is now downregulated and vice versa:
UPMerged1_4 = dfDOWNMerged1_4


#1.4
UnionDF1_4 = pd.concat([DOWNMerged1_3, DOWNMerged1_4], axis=0)
#UnionDF1_4 = UnionDF1_4.drop_duplicates(subset=['GeneID'], keep='first')
UnionDF1_4 = UnionDF1_4.sort_values(by=['logFC'], ascending=True)
#IntersectDF1_4 = pd.merge(dfDOWNMerged1_3, dfDOWNMerged1_4, on='GeneID', how='inner')

############################################################################################################

#1.5
UnionDF1_5 = pd.concat([UPMerged1_3, UPMerged1_4], axis=0)
#UnionDF1_5 = UnionDF1_5.drop_duplicates(subset=['GeneID'], keep='first')
UnionDF1_5 = UnionDF1_5.sort_values(by=['logFC'], ascending=False)

IntersectDF1_5 = pd.merge(UPMerged1_3, UPMerged1_4, on='GeneID', how='inner')
#IntersectDF1_5 = IntersectDF1_5.drop_duplicates(subset=['GeneID'], keep='first')
#IntersectDF1_5 = IntersectDF1_5.sort_values(by=['logFC'], ascending=False) #sort by logFC not working because logFC changed to logFC_x and logFC_y


############################################################################################################

#2.1.2
dfDegUP_C1T5 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_1-T_5_FC_DEGRADERING.xlsx", sheet_name='UP')
dfDegUP_C7T11 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_7-T_11_FC_DEGRADERING.xlsx", sheet_name='UP')

dfDegMerged2_1 = pd.concat([dfDegUP_C1T5, dfDegUP_C7T11], axis=0)
#dfDegMerged2_1 = dfDegMerged2_1.drop_duplicates(subset=['GeneID'], keep='first')
dfDegMerged2_1 = dfDegMerged2_1.sort_values(by=['logFC'], ascending=False)
#Reverse logFC values, logCPM values and sampleValues
dfDegMerged2_1['logFC'] = dfDegMerged2_1['logFC'] * -1
dfDegMerged2_1['logCPM'] = dfDegMerged2_1['logCPM'] * -1
dfDegMerged2_1['Versus'] = dfDegMerged2_1['Versus'].apply(lambda x: '-'.join(x.split('-')[::-1]))
#Since the logFC values are reversed, what was considered upregulated for Cx-Tx is now downregulated and vice versa:
DOWNMerged2_1 = dfDegMerged2_1

#2.1.1
dfDegDOWN_C1T5 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_1-T_5_FC_DEGRADERING.xlsx", sheet_name='DOWN')
dfDegDOWN_C7T11 = pd.read_excel("/home/petear/MEGA/TormodGroup/InputData/C_7-T_11_FC_DEGRADERING.xlsx", sheet_name='DOWN')

dfDegMerged2_2 = pd.concat([dfDegDOWN_C1T5, dfDegDOWN_C7T11], axis=0)
#dfDegMerged2_2 = dfDegMerged2_2.drop_duplicates(subset=['GeneID'], keep='first')
dfDegMerged2_2 = dfDegMerged2_2.sort_values(by=['logFC'], ascending=True)
#Reverse logFC values, logCPM values and sampleValues
dfDegMerged2_2['logFC'] = dfDegMerged2_2['logFC'] * -1
dfDegMerged2_2['logCPM'] = dfDegMerged2_2['logCPM'] * -1
dfDegMerged2_2['Versus'] = dfDegMerged2_2['Versus'].apply(lambda x: '-'.join(x.split('-')[::-1]))
#Since the logFC values are reversed, what was considered upregulated for Cx-Tx is now downregulated and vice versa:
UPMerged2_2 = dfDegMerged2_2

writer = pd.ExcelWriter('AlgoSiScreener.xlsx', engine='xlsxwriter')

UPMerged1_3.to_excel(writer, sheet_name='1.1.1')
DOWNMerged1_3.to_excel(writer, sheet_name='1.1.2')
UnionDF1_3.to_excel(writer, sheet_name='1.3')
UPMerged1_4.to_excel(writer, sheet_name='1.2.1')
DOWNMerged1_4.to_excel(writer, sheet_name='1.2.2')
UnionDF1_4.to_excel(writer, sheet_name='1.4')
UnionDF1_5.to_excel(writer, sheet_name='1.5-Uni')
IntersectDF1_5.to_excel(writer, sheet_name='1.5-Int')
UPMerged2_2.to_excel(writer, sheet_name='2.1.1')
DOWNMerged2_1.to_excel(writer, sheet_name='2.1.2')

writer.save()