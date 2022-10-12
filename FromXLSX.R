library(xlsx)
library(tidyverse)


################################################################################
##For gener TABLE 2BB
#nedregulert i ‘Aβ vs TPA’       <- C7 vs C1  n1
#OG
#nedregulert i ‘LPS+Aβ vs Aβ’    <- T8 vs C7 n10
#OG
#oppregulert i ‘DHA+NS+Aβ vs Aβ’ <- T11 vs C7
################################################################################

c7c1.down <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/C_7_vs_C_1.transcripts.significant.DEA.DOWN.csv", header = TRUE, sep = "\t")
c7c1.down <- c7c1.down %>% rename(posVal = C_7.value, negVal = C_1.value)

t8c7.down <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_8_vs_C_7.transcripts.significant.DEA.DOWN.csv", header=TRUE, sep = "\t")
colnames(t8c7.down) <- colnames(c7c1.down)
 
t11c7.up <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_11_vs_C_7.transcripts.significant.DEA.UP.csv", header=TRUE, sep = "\t")
colnames(t11c7.up) <- colnames(c7c1.down)

fullDF.ddu <- rbind(c7c1.down,t8c7.down,t11c7.up)
fullDF.ddu <- fullDF.ddu[order(fullDF.ddu$qvalue),]

write.xlsx(fullDF.ddu, "/media/petear/SharedPart/Table2B.xlsx")
################################################################################
##For gener TABLE2AA
#OPP i ‘Aβ vs TPA’       <- C7 vs C1  n1
#ELLER
#OPP i ‘LPS+Aβ vs Aβ’    <- T8 vs C7 n10
#ELLER
#NED i ‘DHA+NS+Aβ vs Aβ’ <- T11 vs C7
################################################################################

c7c1.up <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/C_7_vs_C_1.transcripts.significant.DEA.UP.csv", header = TRUE, sep = "\t")
c7c1.up <- c7c1.up %>% rename(posVal = C_7.value, negVal = C_1.value)

t8c7.up <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_8_vs_C_7.transcripts.significant.DEA.UP.csv", header=TRUE, sep = "\t")
colnames(t8c7.up) <- colnames(c7c1.up)

t11c7.down <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_11_vs_C_7.transcripts.significant.DEA.DOWN.csv", header=TRUE, sep = "\t")
colnames(t11c7.down) <- colnames(c7c1.up)

fullDF.uud <- rbind(c7c1.up,t8c7.up,t11c7.down)
fullDF.uud <- fullDF.uud[order(fullDF.uud$qvalue),]
write.xlsx(fullDF.uud, "/media/petear/SharedPart/Table2A.xlsx")


################################################################################
check2a <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DofiTab2a.csv", header=TRUE, sep = ",")
checkVec2a <- c(check2a$Name, check2a$Name.1)
vecFull.uud <- fullDF.uud$gene_name

commonGenesA <- intersect(checkVec2a,vecFull.uud)
DF.uud.ordByVec <- fullDF.uud[match(commonGenesA, fullDF.uud$gene_name),]
fullDF.uud2 <- rbind(DF.uud.ordByVec, fullDF.uud)
write.xlsx(fullDF.uud2, "/media/petear/SharedPart/2ACommonGenesFirst.xlsx")


check2b <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DofiTab2b.csv", header=TRUE, sep = ",")
checkVec2b <- c(check2b$name)
vecFull.ddu <- fullDF.ddu$gene_name

commonGenesB <- intersect(checkVec2b,vecFull.ddu)
DF.ddu.ordByVec <- fullDF.ddu[match(commonGenesB, fullDF.ddu$gene_name),]
fullDF.ddu2 <- rbind(DF.ddu.ordByVec, fullDF.ddu)
write.xlsx(fullDF.ddu2, "/media/petear/SharedPart/2BCommonGenesFirst.xlsx")

################################################################################
#Check directly against Novogens files

##For gener
#nedregulert i ‘Aβ vs TPA’       <- C7 vs C1  n1
#OG
#nedregulert i ‘LPS+Aβ vs Aβ’    <- T8 vs C7 n10
#OG
#oppregulert i ‘DHA+NS+Aβ vs Aβ’ <- T11 vs C7
################################################################################

n1 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/C_7_vs_C_1.transcripts.significant.DEA.DOWN.csv", header=TRUE, sep ="\t")
n1 <- n1 %>% rename(posVal = C_7.value, negVal = C_1.value)

n2 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_8_vs_C_7.transcripts.significant.DEA.DOWN.xls", header=TRUE, sep ="\t")
colnames(n2) <- colnames(n1)

n3 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_11_vs_C_7.transcripts.significant.DEA.UP.xls", header=TRUE, sep ="\t")
colnames(n3) <- colnames(n1)

nAll1_3 <- rbind(n1,n2,n3)
nAll1_3 <- nAll1_3$gene_name
intersect(checkVec2b,nAll1_3)

################################################################################
#Check directly against Novogens files

##For gener
#OPP i ‘Aβ vs TPA’       <- C7 vs C1  n1
#ELLER
#OPP i ‘LPS+Aβ vs Aβ’    <- T8 vs C7 n10
#ELLER
#NED i ‘DHA+NS+Aβ vs Aβ’ <- T11 vs C7
################################################################################
n4 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/C_7_vs_C_1.transcripts.significant.DEA.UP.csv", header=TRUE, sep ="\t")

n5 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_8_vs_C_7.transcripts.significant.DEA.UP.xls", header=TRUE, sep ="\t")
colnames(n5) <- colnames(n4)

n6 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_11_vs_C_7.transcripts.significant.DEA.DOWN.xls", header=TRUE, sep ="\t")
colnames(n6) <- colnames(n4)

nAll4_6 <- rbind(n4,n5,n6)
nAll4_6 <- nAll4_6$gene_name
intersect(checkVec2a,nAll4_6)

################################################################################
#Check directly against Novogens files for ALL

################################################################################

n7 <- read.delim("/media/petear/SharedPart/T11_C7ANDT8C7ANDC1C7.csv", header=TRUE, sep =",")
n7Vec <- n7$gene_name
checkAll <- read.delim("/media/petear/SharedPart/CheckAll.csv", header=TRUE, sep =",")
checkAllVec <- checkAll$Name

intersect(checkAllVec,nAll1_3)

################################################################################
#Check directly against Novogens files for ALL and INSIGNIFICANT 

################################################################################
n8 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/C_7_vs_C_1.transcripts.DEA.csv", header=TRUE, sep ="\t")

n9 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_8_vs_C_7.transcripts.DEA.csv", header=TRUE, sep ="\t")
colnames(n9) <- colnames(n8)

n10 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_11_vs_C_7.transcripts.DEA.csv", header=TRUE, sep ="\t")
colnames(n10) <- colnames(n8)

fullDF.insig <- rbind(n8,n9,n10)
intersect(checkAllVec,fullDF.insig$gene_name)

