import pandas as pd

df1UP_C1T5 = pd.read_excel("/home/petear/New Folder/C_1-T_5_FC_OPPTAK.xlsx", sheet_name='UP')
df2UP_C7T11 = pd.read_excel("/home/petear/New Folder/C_7-T_11_FC_OPPTAK.xlsx", sheet_name='UP')

#1.1.1
dfUPMerged1_3 = pd.concat([df1UP_C1T5, df2UP_C7T11], axis=0)
dfUPMerged1_3 = dfUPMerged1_3.drop_duplicates(subset=['GeneID'], keep='first')
#sort on logFC column in descending order
dfUPMerged1_3 = dfUPMerged1_3.sort_values(by=['logFC'], ascending=False)

df1DOWN_C1T5 = pd.read_excel("/home/petear/New Folder/C_1-T_5_FC_OPPTAK.xlsx", sheet_name='DOWN')
df2DOWN_C7T11 = pd.read_excel("/home/petear/New Folder/C_7-T_11_FC_OPPTAK.xlsx", sheet_name='DOWN')

#1.1.2
dfDOWNMerged1_3 = pd.concat([df1DOWN_C1T5, df2DOWN_C7T11], axis=0)
dfDOWNMerged1_3 = dfDOWNMerged1_3.drop_duplicates(subset=['GeneID'], keep='first')
dfDOWNMerged1_3 = dfDOWNMerged1_3.sort_values(by=['logFC'], ascending=True)
#find genes that are in both dataframes in the ENTREZID column
#1.3
UnionDF1_3 = pd.concat([dfUPMerged1_3, dfDOWNMerged1_3], axis=0)
UnionDF1_3 = UnionDF1_3.drop_duplicates(subset=['GeneID'], keep='first')
UnionDF1_3 = UnionDF1_3.sort_values(by=['logFC'], ascending=False)
#IntersectDF1_3 = pd.merge(dfUPMerged, dfDOWNMerged, on='GeneID', how='inner')


############################################################################################################
df1UP_C1T2 = pd.read_excel("/home/petear/New Folder/C_1-T_2_FC_OPPTAK.xlsx", sheet_name='UP')
df2UP_C7T8 = pd.read_excel("/home/petear/New Folder/C_7-T_8_FC_OPPTAK.xlsx", sheet_name='UP')

#1.2.1
dfUPMerged1_4 = pd.concat([df1UP_C1T2, df2UP_C7T8], axis=0)
dfUPMerged1_4 = dfUPMerged1_4.drop_duplicates(subset=['GeneID'], keep='first')
dfUPMerged1_4 = dfUPMerged1_4.sort_values(by=['logFC'], ascending=False)


df1DOWN_C1T2 = pd.read_excel("/home/petear/New Folder/C_1-T_2_FC_OPPTAK.xlsx", sheet_name='DOWN')
df2DOWN_C7T8 = pd.read_excel("/home/petear/New Folder/C_7-T_8_FC_OPPTAK.xlsx", sheet_name='DOWN')

#1.2.2
dfDOWNMerged1_4 = pd.concat([df1DOWN_C1T2, df2DOWN_C7T8], axis=0)
dfDOWNMerged1_4 = dfDOWNMerged1_4.drop_duplicates(subset=['GeneID'], keep='first')
dfDOWNMerged1_4 = dfDOWNMerged1_4.sort_values(by=['logFC'], ascending=True)

#1.4
UnionDF1_4 = pd.concat([dfDOWNMerged1_3, dfDOWNMerged1_4], axis=0)
UnionDF1_4 = UnionDF1_4.drop_duplicates(subset=['GeneID'], keep='first')
UnionDF1_4 = UnionDF1_4.sort_values(by=['logFC'], ascending=True)
#IntersectDF1_4 = pd.merge(dfDOWNMerged1_3, dfDOWNMerged1_4, on='GeneID', how='inner')

############################################################################################################

#1.5
UnionDF1_5 = pd.concat([dfUPMerged1_3, dfUPMerged1_4], axis=0)
UnionDF1_5 = UnionDF1_5.drop_duplicates(subset=['GeneID'], keep='first')
UnionDF1_5 = UnionDF1_5.sort_values(by=['logFC'], ascending=False)

IntersectDF1_5 = pd.merge(dfUPMerged1_3, dfUPMerged1_4, on='GeneID', how='inner')
IntersectDF1_5 = IntersectDF1_5.drop_duplicates(subset=['GeneID'], keep='first')
IntersectDF1_5 = IntersectDF1_5.sort_values(by=['logFC'], ascending=False)

############################################################################################################

#2.1.1
dfDegUP_C1T5 = pd.read_excel("/home/petear/New Folder/C_1-T_5_FC_DEGRADERING.xlsx", sheet_name='UP')
dfDegUP_C7T11 = pd.read_excel("/home/petear/New Folder/C_7-T_11_FC_DEGRADERING.xlsx", sheet_name='UP')

dfDegMerged2_1 = pd.concat([dfDegUP_C1T5, dfDegUP_C7T11], axis=0)
dfDegMerged2_1 = dfDegMerged2_1.drop_duplicates(subset=['GeneID'], keep='first')
dfDegMerged2_1 = dfDegMerged2_1.sort_values(by=['logFC'], ascending=False)

#2.1.2
dfDegDOWN_C1T5 = pd.read_excel("/home/petear/New Folder/C_1-T_5_FC_DEGRADERING.xlsx", sheet_name='DOWN')
dfDegDOWN_C7T11 = pd.read_excel("/home/petear/New Folder/C_7-T_11_FC_DEGRADERING.xlsx", sheet_name='DOWN')

dfDegMerged2_2 = pd.concat([dfDegDOWN_C1T5, dfDegDOWN_C7T11], axis=0)
dfDegMerged2_2 = dfDegMerged2_2.drop_duplicates(subset=['GeneID'], keep='first')
dfDegMerged2_2 = dfDegMerged2_2.sort_values(by=['logFC'], ascending=True)


writer = pd.ExcelWriter('AlgoSiScreener.xlsx', engine='xlsxwriter')

dfUPMerged1_3.to_excel(writer, sheet_name='1.1.1')
dfDOWNMerged1_3.to_excel(writer, sheet_name='1.1.2')
UnionDF1_3.to_excel(writer, sheet_name='1.3')
dfUPMerged1_4.to_excel(writer, sheet_name='1.2.1')
dfDOWNMerged1_4.to_excel(writer, sheet_name='1.2.2')
UnionDF1_4.to_excel(writer, sheet_name='1.4')
UnionDF1_5.to_excel(writer, sheet_name='1.5-Uni')
IntersectDF1_5.to_excel(writer, sheet_name='1.5-Int')
dfDegMerged2_1.to_excel(writer, sheet_name='2.1.1')
dfDegMerged2_2.to_excel(writer, sheet_name='2.1.2')

writer.save()