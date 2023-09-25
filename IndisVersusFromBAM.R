###TODO:
#CREATE PLOTS FOR GO TERMS COLLECTED:
# Opptak <- c(
# "phagocytosis",
# "endosome",
# "endocytosis",
# "early.endosome",
# "endocytic.vesicle",
# "transport.vesicle"
# )
# Degradering <- c(
# "proteolysis",
# "proteasome-mediated.ubiquitin-dependent.protein.catabolic.process",
# "ubiquitin-dependent.protein.catabolic.process",
# "metallopeptidase.activity",
# "lysosome",
# "endolysosome",
# "lysosomal.protein.catabolic.process"
# )
# Inflammasjon <- c(
# "inflammasome.complex",
# "neuroinflammatory.response",
# "inflammatory.response",
# "autophagy"
# )


# ######Read libraries####
library(edgeR)
library(openxlsx)
library(tidyverse)
library(clusterProfiler)
library("org.Hs.eg.db")
library(biomaRt)
library(Rsubread)
library(feather)
################################################################################
#'*Point code to where files are:*
################################################################################


######Files for analyzos ######
bambamFiles <- list.files("/media/veracrypt11/RNAseq/Bam/",pattern = "\\.bam$")
#GTF is of full
gtfFile <- "/media/veracrypt11/RNAseq/Homo_sapiens.GRCh38.109.chr_patch_hapl_scaff.gtf"

directory <- "/media/veracrypt11/RNAseq/Bam/"

#Add directory to file names for featureCounts
bamFilesDir <- file.path(directory, bambamFiles)

##featCounter <- featureCounts(bamFilesDir, annot.inbuilt = "hg38", isPairedEnd = TRUE, nthreads=18)
#featCounter but with GTF file annot.ext
#featCounter <- featureCounts(bamFilesDir, annot.ext = gtfFile, isPairedEnd = TRUE, nthreads=15, isGTFAnnotationFile = TRUE)

#Already done. Dont need to do it again for now.
featCounter <- readRDS("featCounter.RDS")

fCounterCount <- featCounter$counts


#Create vector for the different groupings for EdgeR (This also keeps them in order of appearance)
grVector <- sapply(colnames(fCounterCount), function(el) substring(el, 1, nchar(el)-5))
grVector <- factor(grVector)
design <- model.matrix(~0+grVector)   #0 is instruction to not add intercept



y <- DGEList(counts=fCounterCount, group=grVector)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

colnames(design) <-  levels(y$samples$group)
y <- estimateDisp(y,design)


#Quasi-likelyhood F-test om ALL:
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit, contrast = c())

################################################################################
#'*PAIRWISE COMPARISONS*
################################################################################

contrastDF <- data.frame(t(combn(unique(grVector),2)))
colnames(contrastDF) <- c("Y1", "Y2")

contrastDF <- contrastDF %>% unite("contCol", Y1:Y2, sep = "-", remove = FALSE)

contrastVec2 <- c(contrastDF$contCol)


contrasts <- makeContrasts(contrasts=contrastVec2, levels=design)

versuses <- colnames(contrasts)


################################################################################
#'*UPREGULATED DEGs*
################################################################################

TippyTopGeneDF_UP <- data.frame(matrix(ncol = 7, nrow = 0)) # create4 vector for list in loop
colnames(TippyTopGeneDF_UP) <- c("ENTREZID", "logFC", "logCPM", "F", "PValue", "FDR", "Versus")


for (versus in versuses)
{
  qlf <- glmQLFTest(fit, contrast=contrasts[,versus])
  TippyTopGenes <- topTags(qlf, sort.by="PValue", n=Inf, p=0.05)
  TippyTopGeneR2C <- tibble::rownames_to_column(data.frame(TippyTopGenes), "ENTREZID")
  #assign(paste0("downReg", versus), listNegPosDFs[[versus]])
  #<- tippyTopTagsR2C[tippyTopTagsR2C$logFC < 0,]
  versVec <- rep(c(versus),times=nrow(TippyTopGeneR2C))
  TippyTopGeneR2C[["Versus"]] <- versVec
  if (nrow(TippyTopGenes)> 0)
  {
    df2<-subset(TippyTopGeneR2C, logFC > 0)
    TippyTopGeneDF_UP <- rbind(TippyTopGeneDF_UP, df2)
  }
}

write_feather(TippyTopGeneDF_UP, "TippyTopGeneDF_UP.feather")

################################################################################
#'*DOWNREGULATED DEGs*
################################################################################

TippyTopGeneDF_DOWN <- data.frame(matrix(ncol = 7, nrow = 0)) # create4 vector for list in loop
colnames(TippyTopGeneDF_DOWN) <- c("ENTREZID", "logFC", "logCPM", "F", "PValue", "FDR", "Versus")


for (versus in versuses)
{
  qlf <- glmQLFTest(fit, contrast=contrasts[,versus])
  TippyTopGenes <- topTags(qlf, sort.by="PValue", n=Inf, p=0.05)
  TippyTopGeneR2C <- tibble::rownames_to_column(data.frame(TippyTopGenes), "ENTREZID")
  #assign(paste0("downReg", versus), listNegPosDFs[[versus]])
  #<- tippyTopTagsR2C[tippyTopTagsR2C$logFC < 0,]
  versVec <- rep(c(versus),times=nrow(TippyTopGeneR2C))
  TippyTopGeneR2C[["Versus"]] <- versVec
  if (nrow(TippyTopGenes)> 0)
  {
    df2<-subset(TippyTopGeneR2C, logFC < 0)
    TippyTopGeneDF_DOWN <- rbind(TippyTopGeneDF_DOWN, df2)
  }
}

write_feather(TippyTopGeneDF_DOWN, "TippyTopGeneDF_DOWN.feather")

################################################################################
#'*ALL DEGs*
################################################################################


TippyTopGeneDF_ALL <- data.frame(matrix(ncol = 7, nrow = 0)) # create4 vector for list in loop
colnames(TippyTopGeneDF_ALL) <- c("ENTREZID", "logFC", "logCPM", "F", "PValue", "FDR", "Versus")


for (versus in versuses)
{
  qlf <- glmQLFTest(fit, contrast=contrasts[,versus])
  TippyTopGenes <- topTags(qlf, sort.by="PValue", n=Inf, p=0.05)
  TippyTopGeneR2C <- tibble::rownames_to_column(data.frame(TippyTopGenes), "ENTREZID")
  #assign(paste0("downReg", versus), listNegPosDFs[[versus]])
  #<- tippyTopTagsR2C[tippyTopTagsR2C$logFC < 0,]
  versVec <- rep(c(versus),times=nrow(TippyTopGeneR2C))
  TippyTopGeneR2C[["Versus"]] <- versVec
  TippyTopGeneDF_ALL <- rbind(TippyTopGeneDF_ALL, TippyTopGeneR2C)
}

write_feather(TippyTopGeneDF_ALL, "TippyTopGeneDF_ALL.feather")

################################################################################
#'*TORMODS GO TERMS <3*
################################################################################

wantedGo <- read.xlsx("/home/petear/CompanDIr/CompFiles/wantedGo.xlsx", sheet = "Sheet1")
wantGoVec <- c(wantedGo$goID)


goTermDF <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(goTermDF) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENTREZID", "ENSEMBL")


#Make a DF of all genes connected with certain Go terms
# for (term in wantGoVec)
# {
#   print(term)
#   rtrv <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=term, columns=c("ENTREZID","ENSEMBL"))
#   goTermDF <- rbind(goTermDF, rtrv)
# }
#Above is not needed. Faster running whole wantGoVec as keys:
rtrv <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=wantGoVec, columns=c("ENTREZID","ENSEMBL"))
goTermDF <- rbind(goTermDF, rtrv)

################################################################################
#'*LOOP OVER TERMS AND DEGs for FULL (ALL)*
################################################################################

#This is for all
uniqGOTermVec <- unique(goTermDF$GOALL)
uniqVers_ALL <- unique(TippyTopGeneDF_ALL$Versus) #Done again because not all versuses might get hits

joinedDF_GO_FULL <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(joinedDF_GO_FULL) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENTREZID", "ENSEMBL", "logFC", "logCPM", "F", "PValue", "FDR", "Versus")

for (goterm in uniqGOTermVec)
{
  dfGO <- goTermDF[goTermDF$GOALL==goterm,]
  for (vers in uniqVers_ALL)
  {
    dfVERS <- TippyTopGeneDF_ALL[TippyTopGeneDF_ALL$Versus==vers,]
    GoVersJoined <- inner_join(dfGO, dfVERS, by = "ENTREZID")
    joinedDF_GO_FULL <- rbind(joinedDF_GO_FULL, GoVersJoined)
  }
}

################################################################################
#'*ADD HGNC SYMBOLS TO joinedDF_GO_FULL (ALL)*
################################################################################
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

valsFull <- joinedDF_GO_FULL$ENSEMBL

convertsFull <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = valsFull, mart = mart)
convertsFull <- convertsFull %>% rename("ensembl_gene_id" = "ENSEMBL", "hgnc_symbol" = "SYMBOL")

joinedDF_GO_FULL <- inner_join(joinedDF_GO_FULL, convertsFull, by = "ENSEMBL")

##SHOULD MAYBE USE joinedDF_GO_FULL INSTEAD OF TippyTopGeneDF_ALL FOR EXPORTING TO PYTHON INSTEAD OF USING BIOMART IN PYTHON
#write_feather(joinedDF_GO_FULL, "joinedDF_GO_FULL.feather")

################################################################################
#'*MAKE PLOTS FOR FULL!*
################################################################################
setwd("/home/petear/")


# Create DF for genes
MeDF_FULL <- data.frame(matrix(ncol = length(uniqVers_ALL), nrow = 500 ))
colnames(MeDF_FULL) <- uniqVers_ALL
BeDF_FULL <- data.frame(matrix(ncol = length(uniqVers_ALL), nrow = 500 ))
colnames(BeDF_FULL) <- uniqVers_ALL
CeDF_FULL <- data.frame(matrix(ncol = length(uniqVers_ALL), nrow = 500 ))
colnames(CeDF_FULL) <- uniqVers_ALL
#Split dataframe based on versus in "Versus"-column#
#SplitVers <- split(joinedDF_GO_FULL, joinedDF_GO_FULL$Versus)    # Split data frame in list

#okay, PLOT AWAY!

#Change name to what you want the plot to be named
pdf("VersPlotsALL.pdf", width = 20, height =20)

for (el in uniqVers_ALL) 
{
  meg <- subset(joinedDF_GO_FULL, Versus == el)
  meg <- unique(meg$SYMBOL)
  identifier = el
  MeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL")
  Melk <- as.data.frame(MeGo)
  Mel <- c(Melk$geneID[1])
  ForCol <- str_split(Mel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  MeDF_FULL[el] <- ForCol
  if (nrow(Melk) > 0)
  {
    print(cnetplot(MeGo, color_category='#1b9e77',
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(MeGo$ID) >1)
    {
      print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  BeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL")
  Belk <- as.data.frame(BeGo)
  Bel <- c(Belk$geneID[1])
  ForCol <- str_split(Bel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  BeDF_FULL[el] <- ForCol
  if (nrow(Belk) > 0)
  {
    print(cnetplot(BeGo, color_category='#1b9e77',
                 color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(BeGo$ID) >1)
    {
      print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  CeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL")
  Celk <- as.data.frame(CeGo)
  Cel <- c(Celk$geneID[1])
  ForCol <- str_split(Cel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  CeDF_FULL[el] <- ForCol
  if (nrow(Celk) > 0)
  {
    print(cnetplot(CeGo, color_category='#1b9e77', 
                 color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(CeGo$ID) >1)
    {
    print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  AeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "SYMBOL")
  Aelk <- as.data.frame(AeGo) #to check for hits
  if (nrow(Aelk) > 0)
  {
    print(cnetplot(AeGo, color_category='#1b9e77', 
                 color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  }
}


dev.off()

################################################################################
#'*LOOP OVER TERMS AND DEGs for DOWN*
################################################################################

#This is for down
uniqGOTermVec <- unique(goTermDF$GOALL)
uniqVers_DOWN <- unique(TippyTopGeneDF_DOWN$Versus) #Done again because not all versuses might get hits

joinedDF_GO_DOWN <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(joinedDF_GO_DOWN) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENTREZID", "ENSEMBL", "logFC", "logCPM", "F", "PValue", "FDR", "Versus")

for (goterm in uniqGOTermVec)
{
  dfGO <- goTermDF[goTermDF$GOALL==goterm,]
  for (vers in uniqVers_DOWN)
  {
    dfVERS <- TippyTopGeneDF_DOWN[TippyTopGeneDF_DOWN$Versus==vers,]
    GoVersJoined <- inner_join(dfGO, dfVERS, by = "ENTREZID")
    joinedDF_GO_DOWN <- rbind(joinedDF_GO_DOWN, GoVersJoined)
  }
}

################################################################################
#'*ADD HGNC SYMBOLS TO joinedDF_GO_DOWN*
################################################################################
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

valsFull <- joinedDF_GO_DOWN$ENSEMBL

convertsFull <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = valsFull, mart = mart)
convertsFull <- convertsFull %>% rename("ensembl_gene_id" = "ENSEMBL", "hgnc_symbol" = "SYMBOL")

joinedDF_GO_DOWN <- inner_join(joinedDF_GO_DOWN, convertsFull, by = "ENSEMBL")

################################################################################
#'*MAKE PLOTS FOR DOWN!*
################################################################################

# Create DF for genes
MeDF_DOWN <- data.frame(matrix(ncol = length(uniqVers_DOWN), nrow = 500 ))
colnames(MeDF_DOWN) <- uniqVers_DOWN
BeDF_DOWN <- data.frame(matrix(ncol = length(uniqVers_DOWN), nrow = 500 ))
colnames(BeDF_DOWN) <- uniqVers_DOWN
CeDF_DOWN <- data.frame(matrix(ncol = length(uniqVers_DOWN), nrow = 500 ))
colnames(CeDF_DOWN) <- uniqVers_DOWN

#Split dataframe based on versus in "Versus"-column#
SplitVers <- split(joinedDF_GO_DOWN, joinedDF_GO_DOWN$Versus)    # Split data frame in list

#okay, PLOT AWAY!

#Change name to what you want the plot to be named
pdf("VersPlotsDown.pdf", width = 20, height =20)

for (el in uniqVers_DOWN) 
{
  meg <- subset(joinedDF_GO_DOWN, Versus == el)
  meg <- unique(meg$SYMBOL)
  identifier = el
  MeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL")
  Melk <- as.data.frame(MeGo)
  Mel <- c(Melk$geneID[1])
  ForCol <- str_split(Mel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  MeDF_DOWN[el] <- ForCol
  if (nrow(Melk) > 0)
  {
    print(cnetplot(MeGo, color_category='#1b9e77',
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(MeGo$ID) >1)
    {
      print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  BeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL")
  Belk <- as.data.frame(BeGo)
  Bel <- c(Belk$geneID[1])
  ForCol <- str_split(Bel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  BeDF_DOWN[el] <- ForCol
  if (nrow(Belk) > 0)
  {
    print(cnetplot(BeGo, color_category='#1b9e77',
                 color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(BeGo$ID) >1)
    {
      print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    } 
  }

  CeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL")
  Celk <- as.data.frame(CeGo)
  Cel <- c(Celk$geneID[1])
  ForCol <- str_split(Cel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  CeDF_DOWN[el] <- ForCol
  if (nrow(Celk) > 0)
  {
    print(cnetplot(CeGo, color_category='#1b9e77', 
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(CeGo$ID) >1)
    {
      print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  AeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "SYMBOL")
  Aelk <- as.data.frame(AeGo) #to check for hits
  if (nrow(Aelk) > 0)
  {
    print(cnetplot(AeGo, color_category='#1b9e77', 
                 color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  }
}

dev.off()



################################################################################
#'*LOOP OVER TERMS AND DEGs for UP*
################################################################################

#This is for all
uniqGOTermVec <- unique(goTermDF$GOALL)
uniqVers_UP <- unique(TippyTopGeneDF_UP$Versus) #Done again because not all versuses might get hits

joinedDF_GO_UP <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(joinedDF_GO_UP) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENTREZID", "ENSEMBL", "logFC", "logCPM", "F", "PValue", "FDR", "Versus")

for (goterm in uniqGOTermVec)
{
  dfGO <- goTermDF[goTermDF$GOALL==goterm,]
  for (vers in uniqVers_UP)
  {
    dfVERS <- TippyTopGeneDF_UP[TippyTopGeneDF_UP$Versus==vers,]
    GoVersJoined <- inner_join(dfGO, dfVERS, by = "ENTREZID")
    joinedDF_GO_UP <- rbind(joinedDF_GO_UP, GoVersJoined)
  }
}

################################################################################
#'*ADD HGNC SYMBOLS TO joinedDF_GO_UP*
################################################################################
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

valsFull <- joinedDF_GO_UP$ENSEMBL

convertsFull <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = valsFull, mart = mart)
convertsFull <- convertsFull %>% rename("ensembl_gene_id" = "ENSEMBL", "hgnc_symbol" = "SYMBOL")

joinedDF_GO_UP <- inner_join(joinedDF_GO_UP, convertsFull, by = "ENSEMBL")

################################################################################
#'*MAKE PLOTS FOR UP!*
################################################################################


# Create DF for genes
MeDF_UP <- data.frame(matrix(ncol = length(uniqVers_UP), nrow = 500 ))
colnames(MeDF_UP) <- uniqVers_UP
BeDF_UP <- data.frame(matrix(ncol = length(uniqVers_UP), nrow = 500 ))
colnames(BeDF_UP) <- uniqVers_UP
CeDF_UP <- data.frame(matrix(ncol = length(uniqVers_UP), nrow = 500 ))
colnames(CeDF_UP) <- uniqVers_UP
#okay, PLOT AWAY!

#Change name to what you want the plot to be named
pdf("VersPlotsUp.pdf", width = 20, height =20)

for (el in uniqVers_UP) 
{
  meg <- subset(joinedDF_GO_UP, Versus == el)
  meg <- unique(meg$SYMBOL)
  identifier = el
  MeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL")
  Melk <- as.data.frame(MeGo)
  Mel <- c(Melk$geneID[1])
  ForCol <- str_split(Mel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  MeDF_UP[el] <- ForCol
  if (nrow(Melk) > 0)
  {
    print(cnetplot(MeGo, color_category='#1b9e77',
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(MeGo$ID) >1)
    {
    print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  BeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL")
  Belk <- as.data.frame(BeGo)
  Bel <- c(Belk$geneID[1])
  ForCol <- str_split(Bel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  BeDF_UP[el] <- ForCol
  if (nrow(Belk) > 0)
  {
    print(cnetplot(BeGo, color_category='#1b9e77',
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(BeGo$ID) >1)
    {
    print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  CeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL")
  Celk <- as.data.frame(CeGo)
  Cel <- c(Celk$geneID[1])
  ForCol <- str_split(Cel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  CeDF_UP[el] <- ForCol
  if (nrow(Celk) > 0)
  {
    print(cnetplot(CeGo, color_category='#1b9e77', 
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(CeGo$ID) >1)
    {
    print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }
  AeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "SYMBOL")
  Aelk <- as.data.frame(AeGo) #to check for hits
  if (nrow(Aelk) > 0)
  {
    print(cnetplot(AeGo, color_category='#1b9e77', 
                 color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  }
}



dev.off()

################################################################################
#'*WRITE THE DFs TO EXCEL * 
################################################################################

write.xlsx(MeDF_FULL, file = "MeDF_FULL.xlsx", row.names = FALSE)
write.xlsx(BeDF_FULL, file = "BeDF_FULL.xlsx", row.names = FALSE)
write.xlsx(CeDF_FULL, file = "CeDF_FULL.xlsx", row.names = FALSE)

write.xlsx(MeDF_DOWN, file = "MeDF_DOWN.xlsx", row.names = FALSE)
write.xlsx(BeDF_DOWN, file = "BeDF_DOWN.xlsx", row.names = FALSE)
write.xlsx(CeDF_DOWN, file = "CeDF_DOWN.xlsx", row.names = FALSE)

write.xlsx(MeDF_UP, file = "MeDF_UP.xlsx", row.names = FALSE)
write.xlsx(BeDF_UP, file = "BeDF_UP.xlsx", row.names = FALSE)
write.xlsx(CeDF_UP, file = "CeDF_UP.xlsx", row.names = FALSE)




################################################################################
#'*PATHWAY ANALYSIS FOR DOWN * 
################################################################################


################################################################################
#'*LOOP OVER EACH GO TERM * 
################################################################################

save.image("StartingGoPlots.RData")

uniGO_UP <- unique(joinedDF_GO_UP$GOALL)

MeDF_UP_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))
BeDF_UP_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))
CeDF_UP_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))



pdf("VersPlotsUpPerGO.pdf", width = 20, height =20)

for (el in uniqVers_UP) 
{
  meg2 <- subset(joinedDF_GO_UP, Versus == el)
  for (gop in uniGO_UP)
  {
  meg3 <- subset(meg2, GOALL == gop)
  meg <- unique(meg3$SYMBOL)
  #finding go term as keywords:
  wgo <- wantedGo[which(wantedGo$goID == gop), ]
  wgo <- wgo$rogo
  identifier = paste(el, "for", gop, "(", wgo,")")
  MeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL")
  Melk <- as.data.frame(MeGo)
  Mel <- c(Melk$geneID[1])
  ForCol <- str_split(Mel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  MeDF_UP_GO[identifier] <- ForCol
  if (nrow(Melk) > 0)
  {
    print(cnetplot(MeGo, color_category='#1b9e77',
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(MeGo$ID) >1)
    {
    print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  BeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL")
  Belk <- as.data.frame(BeGo)
  Bel <- c(Belk$geneID[1])
  ForCol <- str_split(Bel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  BeDF_UP_GO[identifier] <- ForCol
  if (nrow(Belk) > 0)
  {
    print(cnetplot(BeGo, color_category='#1b9e77',
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(BeGo$ID) >1)
    {
    print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  CeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL")
  Celk <- as.data.frame(CeGo)
  Cel <- c(Celk$geneID[1])
  ForCol <- str_split(Cel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  CeDF_UP_GO[identifier] <- ForCol
  if (nrow(Celk) > 0)
  {
    print(cnetplot(CeGo, color_category='#1b9e77', 
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(CeGo$ID) >1)
    {
    print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }
  AeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "SYMBOL")
  Aelk <- as.data.frame(AeGo) #to check for hits
  if (nrow(Aelk) > 0)
  {
    print(cnetplot(AeGo, color_category='#1b9e77', 
                 color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  }
}
}

dev.off()

save.image("FinishedUpGOplots.RData")

uniGO_DOWN <- unique(joinedDF_GO_DOWN$GOALL)

MeDF_DOWN_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))
BeDF_DOWN_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))
CeDF_DOWN_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))


pdf("VersPlotsDownPerGO.pdf", width = 20, height =20)

for (el in uniqVers_DOWN) 
{
  meg2 <- subset(joinedDF_GO_DOWN, Versus == el)
  for (gop in uniGO_DOWN)
  {
  meg3 <- subset(meg2, GOALL == gop)
  meg <- unique(meg3$SYMBOL)
  #finding go term as keywords:
  wgo <- wantedGo[which(wantedGo$goID == gop), ]
  wgo <- wgo$rogo
  identifier = paste(el, "for", gop, "(", wgo,")")
  MeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL")
  Melk <- as.data.frame(MeGo)
  Mel <- c(Melk$geneID[1])
  ForCol <- str_split(Mel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  MeDF_DOWN_GO[identifier] <- ForCol
  if (nrow(Melk) > 0)
  {
    print(cnetplot(MeGo, color_category='#1b9e77',
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(MeGo$ID) >1)
    {
      print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  BeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL")
  Belk <- as.data.frame(BeGo)
  Bel <- c(Belk$geneID[1])
  ForCol <- str_split(Bel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  BeDF_DOWN_GO[identifier] <- ForCol
  if (nrow(Belk) > 0)
  {
    print(cnetplot(BeGo, color_category='#1b9e77',
                 color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(BeGo$ID) >1)
    {
      print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    } 
  }

  CeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL")
  Celk <- as.data.frame(CeGo)
  Cel <- c(Celk$geneID[1])
  ForCol <- str_split(Cel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  CeDF_DOWN_GO[identifier] <- ForCol
  if (nrow(Celk) > 0)
  {
    print(cnetplot(CeGo, color_category='#1b9e77', 
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(CeGo$ID) >1)
    {
      print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  AeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "SYMBOL")
  Aelk <- as.data.frame(AeGo) #to check for hits
  if (nrow(Aelk) > 0)
  {
    print(cnetplot(AeGo, color_category='#1b9e77', 
                 color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  }
}
}

dev.off()

save.image("FinishedDownGOplots.RData")

uniGO_FULL <- unique(joinedDF_GO_FULL$GOALL)

MeDF_FULL_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))
BeDF_FULL_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))
CeDF_FULL_GO <- data.frame(matrix(ncol = 0, nrow = 500 ))


pdf("VersPlotsALLPerGO.pdf", width = 20, height =20)

for (el in uniqVers_ALL) 
{
  meg2 <- subset(joinedDF_GO_FULL, Versus == el)
  for (gop in uniGO_FULL)
  {
  meg3 <- subset(meg2, GOALL == gop)
  meg <- unique(meg3$SYMBOL)
  #finding go term as keywords:
  wgo <- wantedGo[which(wantedGo$goID == gop), ]
  wgo <- wgo$rogo
  identifier = paste(el, "for", gop, "(", wgo,")")
  MeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL")
  Melk <- as.data.frame(MeGo)
  Mel <- c(Melk$geneID[1])
  ForCol <- str_split(Mel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  MeDF_FULL_GO[identifier] <- ForCol
  if (nrow(Melk) > 0)
  {
    print(cnetplot(MeGo, color_category='#1b9e77',
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(MeGo$ID) >1)
    {
      print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  BeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL")
  Belk <- as.data.frame(BeGo)
  Bel <- c(Belk$geneID[1])
  ForCol <- str_split(Bel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  BeDF_FULL_GO[identifier] <- ForCol
  if (nrow(Belk) > 0)
  {
    print(cnetplot(BeGo, color_category='#1b9e77',
                 color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process:",identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(BeGo$ID) >1)
    {
      print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  CeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL")
  Celk <- as.data.frame(CeGo)
  Cel <- c(Celk$geneID[1])
  ForCol <- str_split(Cel, "/")
  ForCol <- unlist(ForCol)
  while (length(ForCol) < 500)
  {
    ForCol <- append(ForCol, NA)
  }
  CeDF_FULL_GO[identifier] <- ForCol
  if (nrow(Celk) > 0)
  {
    print(cnetplot(CeGo, color_category='#1b9e77', 
                 color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(CeGo$ID) >1)
    {
    print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    }
  }

  AeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "SYMBOL")
  Aelk <- as.data.frame(AeGo) #to check for hits
  if (nrow(Aelk) > 0)
  {
    print(cnetplot(AeGo, color_category='#1b9e77', 
                 color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  }
}
}

dev.off()
save.image("FinishedallGOplots.RData")

################################################################################
#'*WRITE THE DFs FOR PER GO TO FILE * 
################################################################################


write.csv(MeDF_FULL_GO, "MeDF_FULL_perGO.csv", row.names = FALSE, na="")
write.csv(BeDF_FULL_GO,"BeDF_FULL_perGO.csv", row.names = FALSE, na="")
write.csv(CeDF_FULL_GO,"CeDF_FULL_perGO.csv", row.names = FALSE, na="")

write.csv(MeDF_DOWN_GO,"MeDF_DOWN_perGO.csv", row.names = FALSE, na="")
write.csv(BeDF_DOWN_GO,"BeDF_DOWN_perGO.csv", row.names = FALSE, na="")
write.csv(CeDF_DOWN_GO,"CeDF_DOWN_perGO.csv", row.names = FALSE, na="")

write.csv(MeDF_UP_GO, "MeDF_UP_perGO.csv", row.names = FALSE, na="")
write.csv(BeDF_UP_GO,"BeDF_UP_perGO.csv", row.names = FALSE, na="")
write.csv(CeDF_UP_GO,"CeDF_UP_perGO.csv", row.names = FALSE, na="")


# write.xlsx2(MeDF_FULL_GO,"MeDF_FULL_perGO.xlsx", row.names = FALSE)
# write.xlsx2(BeDF_FULL_GO,"BeDF_FULL_perGO.xlsx", row.names = FALSE)
# write.xlsx2(CeDF_FULL_GO,"CeDF_FULL_perGO.xlsx", row.names = FALSE)

# write.xlsx2(MeDF_DOWN_GO,"MeDF_DOWN_perGO.xlsx", row.names = FALSE)
# write.xlsx2(BeDF_DOWN_GO,"BeDF_DOWN_perGO.xlsx", row.names = FALSE)
# write.xlsx2(CeDF_DOWN_GO,"CeDF_DOWN_perGO.xlsx", row.names = FALSE)

# write.xlsx2(MeDF_UP_GO,"MeDF_UP_perGO.xlsx", row.names = FALSE)
# write.xlsx2(BeDF_UP_GO,"BeDF_UP_perGO.xlsx", row.names = FALSE)
# write.xlsx2(CeDF_UP_GO,"CeDF_UP_perGO.xlsx", row.names = FALSE)

