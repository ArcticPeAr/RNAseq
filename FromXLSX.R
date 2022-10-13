library(xlsx)
library(tidyverse)
library("org.Hs.eg.db")

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
#This section creates check files from the manuscript
# And then checks files against a vector of gene names from above 
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

#UPDATED FILES!!
################################################################################

n1 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/C_7_vs_C_1.transcripts.significant.DEA.DOWN.xls", header=TRUE, sep ="\t")
n1 <- n1 %>% rename(posVal = C_7.value, negVal = C_1.value)

n2 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/T_8_vs_C_7.transcripts.significant.DEA.DOWN.xls", header=TRUE, sep ="\t")
colnames(n2) <- colnames(n1)

n3 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/T_11_vs_C_7.transcripts.significant.DEA.UP.xls", header=TRUE, sep ="\t")
colnames(n3) <- colnames(n1)

nAll1_3 <- rbind(n1,n2,n3)
nAll1_3GN <- nAll1_3$gene_name
intersect(checkVec2b,nAll1_3GN)

################################################################################
#Check directly against Novogens files

##For gener
#OPP i ‘Aβ vs TPA’       <- C7 vs C1  n1
#ELLER
#OPP i ‘LPS+Aβ vs Aβ’    <- T8 vs C7 n10
#ELLER
#NED i ‘DHA+NS+Aβ vs Aβ’ <- T11 vs C7

#UPDATED FILES!!
################################################################################
n4 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/C_7_vs_C_1.transcripts.significant.DEA.UP.xls", header=TRUE, sep ="\t")

n5 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/T_8_vs_C_7.transcripts.significant.DEA.UP.xls", header=TRUE, sep ="\t")
colnames(n5) <- colnames(n4)

n6 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/T_11_vs_C_7.transcripts.significant.DEA.DOWN.xls", header=TRUE, sep ="\t")
colnames(n6) <- colnames(n4)

nAll4_6 <- rbind(n4,n5,n6)
nAll4_6GN <- nAll4_6$gene_name
intersect(checkVec2a,nAll4_6GN)

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
#NOT NEEDED!
################################################################################
n8 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/C_7_vs_C_1.transcripts.DEA.csv", header=TRUE, sep ="\t")

n9 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_8_vs_C_7.transcripts.DEA.csv", header=TRUE, sep ="\t")
colnames(n9) <- colnames(n8)

n10 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/T_11_vs_C_7.transcripts.DEA.csv", header=TRUE, sep ="\t")
colnames(n10) <- colnames(n8)

fullDF.insig <- rbind(n8,n9,n10)
intersect(checkAllVec,fullDF.insig$gene_name)

################################################################################
#Check with GO that Tormod liker
#
################################################################################
wantedGo <- read.xlsx2("/media/petear/SharedPart/wantedGo.xlsx", sheetIndex = 1)
wantGoVec <- c(wantedGo$goID)



goTermDF <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(goTermDF) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENSEMBL")

#Make a DF of all genes connected with certain Go terms
for (term in wantGoVec)
{
  rtrv <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=term, columns="ENSEMBL")
  goTermDF <- rbind(goTermDF, rtrv)
}

#nAll4_6
#nAll1_3

names(nAll1_3)[names(nAll1_3) == "gene_id"] <- "ENSEMBL"
names(nAll4_6)[names(nAll4_6) == "gene_id"] <- "ENSEMBL"

n1_3Joined <- inner_join(goTermDF, nAll1_3, by = "ENSEMBL")
n4_6Joined <- inner_join(goTermDF, nAll4_6, by = "ENSEMBL")

#For just the 80 + 63 terms 

n1_3JoinSorted <- n1_3Joined[order(n1_3Joined$pvalue),]
n4_6JoinSorted <- n4_6Joined[order(n4_6Joined$pvalue),]

write.xlsx(n1_3JoinSorted, "/media/petear/SharedPart/Ab_TPA-DownLPSAb_Ab-DownDHANSAb_Ab-Up.xlsx")
write.xlsx(n4_6JoinSorted, "/media/petear/SharedPart/Ab_TPA-UpLPSAb_Ab-UpDHANSAb_Ab-Down.xlsx")


################################################################################
#make plots
################################################################################
plotList <- list("n1_3JoinSorted",n1_3JoinSorted, "n4_6JoinSorted",n4_6JoinSorted)
i = 1
pdf("test.pdf", width = 20, height =20)
for (el in plotList) 
{
      if ((i %% 2) != 0)
      {
        identifier = el
      i=i+1
      }
      else
      {
        MeGo <- enrichGO(el$ENSEMBL, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "ENSEMBL")
        print(cnetplot(MeGo, color_category='#1b9e77',
                       color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        
        BeGo <- enrichGO(el$ENSEMBL, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENSEMBL")
        print(cnetplot(BeGo, color_category='#1b9e77',
                       color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        
        CeGo <- enrichGO(el$ENSEMBL, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "ENSEMBL")
        print(cnetplot(CeGo, color_category='#1b9e77', 
                       color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        
        AeGo <- enrichGO(el$ENSEMBL, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "ENSEMBL")
        print(cnetplot(AeGo, color_category='#1b9e77', 
                       color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        i=i+1
      }
}

dev.off()



