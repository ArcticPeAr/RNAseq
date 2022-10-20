library(xlsx)
library(tidyverse)
library("org.Hs.eg.db")
library(biomaRt)
library(rWikiPathways)
library(clusterProfiler)
library(stringr) # to split strings
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
n1 <- n1 %>% rename("posVal" = "C_7.value", "negVal" = "C_1.value")

n2 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/T_8_vs_C_7.transcripts.significant.DEA.DOWN.xls", header=TRUE, sep ="\t")
colnames(n2) <- colnames(n1)

n3 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/T_11_vs_C_7.transcripts.significant.DEA.UP.xls", header=TRUE, sep ="\t")
colnames(n3) <- colnames(n1)

nAll1_3 <- rbind(n1,n2,n3)
nAll1_3GN <- nAll1_3$gene_name
intersect(checkVec2b,nAll1_3GN)

names(nAll1_3)[names(nAll1_3) == "gene_id"] <- "ENSEMBL"

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
n4 <- n4 %>% rename("posVal" = "C_7.value", "negVal" = "C_1.value")

n5 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/T_8_vs_C_7.transcripts.significant.DEA.UP.xls", header=TRUE, sep ="\t")
colnames(n5) <- colnames(n4)

n6 <- read.delim("/media/petear/SharedPart/RNAseq/Result_X204SC21061474-Z01-F002_Homo_sapiens/7.DiffExprAnalysis/transcript/DEGlist/DEGlistUpdated/T_11_vs_C_7.transcripts.significant.DEA.DOWN.xls", header=TRUE, sep ="\t")
colnames(n6) <- colnames(n4)

nAll4_6 <- rbind(n4,n5,n6)
nAll4_6GN <- nAll4_6$gene_name
intersect(checkVec2a,nAll4_6GN)

names(nAll4_6)[names(nAll4_6) == "gene_id"] <- "ENSEMBL"

################################################################################
#'*Check directly against Novogens files for ALL*
################################################################################

n7 <- read.delim("/media/petear/SharedPart/T11_C7ANDT8C7ANDC1C7.csv", header=TRUE, sep =",")
n7Vec <- n7$gene_name
checkAll <- read.delim("/media/petear/SharedPart/CheckAll.csv", header=TRUE, sep =",")
checkAllVec <- checkAll$Name

intersect(checkAllVec,nAll1_3)

################################################################################
#'*Check with GO that Tormod liker*
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

n1_3Joined <- inner_join(goTermDF, nAll1_3, by = "ENSEMBL")
n4_6Joined <- inner_join(goTermDF, nAll4_6, by = "ENSEMBL")

#For just the 80 + 63 terms 

n1_3JoinSorted <- n1_3Joined[order(n1_3Joined$pvalue),]
n4_6JoinSorted <- n4_6Joined[order(n4_6Joined$pvalue),]

write.xlsx(n1_3JoinSorted, "/media/petear/SharedPart/Ab_TPA-DownLPSAb_Ab-DownDHANSAb_Ab-Up.xlsx")
write.xlsx(n4_6JoinSorted, "/media/petear/SharedPart/Ab_TPA-UpLPSAb_Ab-UpDHANSAb_Ab-Down.xlsx")


################################################################################
#'*Go term Kulb wanted: GO:0048156*
################################################################################
kulbGoTerm <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(kulbGoTerm) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENSEMBL")

rtrv <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0048156", columns="ENSEMBL")
kulbGoTerm <- rbind(kulbGoTerm, rtrv)

kulbGon1_3Joined <- inner_join(kulbGoTerm, nAll1_3, by = "ENSEMBL")
kulbGon4_6Joined <- inner_join(kulbGoTerm, nAll4_6, by = "ENSEMBL")

write.xlsx(kulbGon1_3Joined, "/media/petear/SharedPart/Ab_TPA-DownLPSAb_Ab-DownDHANSAb_Ab-Up.xlsx")
write.xlsx(kulbGon4_6Joined, "/media/petear/SharedPart/Ab_TPA-UpLPSAb_Ab-UpDHANSAb_Ab-Down.xlsx")

################################################################################
#'*Make plots*
################################################################################
setwd("/media/petear/SharedPart/Plots/")

#Edit plotList to you hearts content!
plotList <- list("kulbGon1_3Joined",kulbGon1_3Joined, "kulbGon4_6Joined",kulbGon4_6Joined) #Add here what was taken out last time you ran plotting: "n1_3JoinSorted",n1_3JoinSorted, "n4_6JoinSorted",n4_6JoinSorted
#Did you edit?


# Create vectors for genes
MeVec <- c()
BeVec <- c()
CeVec <- c()
#okay, PLOT AWAY!
i = 1
pdf("TauBindProt_n1_3_n4_6.pdf", width = 20, height =20)
## SYMBOL -- ENTREZID -- ENSEMBL
## gene_name -- NA -- ENSEMBL

for (el in plotList) 
{
      if ((i %% 2) != 0)
      {
        identifier = el
      i=i+1
      }
      else
      {
        MeGo <- enrichGO(el$gene_name, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL")
        Mel <- as.data.frame(MeGo)
        Mel <- c(Mel$geneID[1])
        #col <- str_split(col, "/")
        MeVec <- append(MeVec, Mel)
        print(cnetplot(MeGo, color_category='#1b9e77',
                       color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm'))) +scale_color_gradient(low="#b2df8a", high="#1f78b4")
        
        BeGo <- enrichGO(el$gene_name, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL")
        Bel <- as.data.frame(BeGo)
        Bel <- c(Bel$geneID[1])
        #col <- str_split(col, "/")
        BeVec <- append(BeVec, Bel)
        print(cnetplot(BeGo, color_category='#1b9e77',
                       color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm'))) +scale_color_gradient(low="#b2df8a", high="#1f78b4")
        
        CeGo <- enrichGO(el$gene_name, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL")
        Cel <- as.data.frame(CeGo)
        Cel <- c(Cel$geneID[1])
        #col <- str_split(col, "/")
        CeVec <- append(CeVec, Cel)
        print(cnetplot(CeGo, color_category='#1b9e77', 
                       color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
        print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm'))) +scale_color_gradient(low="#b2df8a", high="#1f78b4")
        
        AeGo <- enrichGO(el$gene_name, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "SYMBOL")
        print(cnetplot(AeGo, color_category='#1b9e77', 
                       color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm'))) +scale_color_gradient(low="#b2df8a", high="#1f78b4")
        i=i+1
      }
}

dev.off()


################################################################################
#'*Convert from ENSEMBL to ENTREZ*
################################################################################

marty <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                 dataset = "hsapiens_gene_ensembl",
                 host = "https://www.ensembl.org")

genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","entrezgene_id"),
               values = n1_3JoinSorted$ENSEMBL, 
               mart = marty)

################################################################################
#'*Pathway analysis*
################################################################################






pdf("KEGG.pdf")
eKEGG <- enrichWP(n1_3JoinSorted$ENSEMBL,keyType = "ENSEMBL")
cnetplot(eKEGG, color_category='#1b9e77', 
         color_gene='#d95f02')  + ggtitle("KEGG") + theme(plot.margin=unit(c(0.0,0.2,0.0,0.2), 'cm'))
cnetplot(eKEGG, showCategory = 22, color_category='#1b9e77', 
         color_gene='#d95f02') + ggtitle("KEGG -forced extended connections") + theme(plot.margin=unit(c(0.0,0.2,0.0,0.2), 'cm'))

