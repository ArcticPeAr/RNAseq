######Read libraries####
library(edgeR)
library(xlsx)
library(tidyverse)
library(clusterProfiler)
library("org.Hs.eg.db")
library(biomaRt)
################################################################################
#'*Point code to where files are:*
################################################################################
setwd("/media/petear/SharedPart/RNAseq/Bam/")


######Files for analyzos ######
bambamFiles <- list.files(pattern = "\\.bam$")

#GTF is of full
gtfFile <- "/home/petear/Documents/Homo_sapiens.GRCh38.107.chr_patch_hapl_scaff.gtf"

##featCounter <- featureCounts(bambamFiles, annot.inbuilt = "hg38", isPairedEnd = TRUE, nthreads=18)
#Already done. Dont need to do it again for now.
featCounter <- readRDS("featCounter.RDS")

fCounterCount <- featCounter$counts


#Create vector for the different groupings for EdgeR (This also keeps them in order of appearance)
grVector <- c()
for (el in colnames(fCounterCount)) {
  grMem = substring(el, 1, nchar(el)-5)
  grVector[length(grVector) + 1] <- grMem
}
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


################################################################################
#'*TORMODS GO TERMS <3*
################################################################################

wantedGo <- read.xlsx2("/media/petear/SharedPart/wantedGo.xlsx", sheetIndex = 1)
wantGoVec <- c(wantedGo$goID)


goTermDF <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(goTermDF) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENTREZID", "ENSEMBL")


#Make a DF of all genes connected with certain Go terms
for (term in wantGoVec)
{
  rtrv <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=term, columns=c("ENTREZID","ENSEMBL"))
  goTermDF <- rbind(goTermDF, rtrv)
}

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

################################################################################
#'*MAKE PLOTS FOR FULL!*
################################################################################
setwd("/home/petear/")


# Create vectors for genes
MeVec_FULL <- c()
BeVec_FULL <- c()
CeVec_FULL <- c()

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
  if (nrow(Melk) > 0)
  {
    MeVec_FULL <- append(MeVec_FULL, el)
    MeVec_FULL <- append(MeVec_FULL, Mel)
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
  if (nrow(Belk) > 0)
  {
    BeVec_FULL <- append(BeVec_FULL, el)
    BeVec_FULL <- append(BeVec_FULL, Bel)
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
  if (nrow(Celk) > 0)
  {
    CeVec_FULL <- append(CeVec_FULL, el)
    CeVec_FULL <- append(CeVec_FULL, Cel)
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

# Create vectors for genes
MeVec_DOWN <- c()
BeVec_DOWN <- c()
CeVec_DOWN <- c()

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
  if (nrow(Melk) > 0)
  {
    MeVec_DOWM <- append(MeVec_DOWN, el)
    MeVec_DOWN <- append(MeVec_DOWN, Mel)
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
  if (nrow(Belk) > 0)
  {
    BeVec_DOWM <- append(BeVec_DOWN, el)
    BeVec_DOWN <- append(BeVec_DOWN, Bel)
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
  
  if (nrow(Celk) > 0)
  {
    CeVec_DOWM <- append(CeVec_DOWN, el)
    CeVec_DOWN <- append(CeVec_DOWN, Cel)
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


# Create vectors for genes
MeVec_UP <- c()
BeVec_UP <- c()
CeVec_UP <- c()


#okay, PLOT AWAY!

#Change name to what you want the plot to be named
pdf("VersPlotsUp.pdf", width = 20, height =20)

for (el in uniqVers_UP) 
{
  cat(el)
  meg <- subset(joinedDF_GO_UP, Versus == el)
  meg <- unique(meg$SYMBOL)
  identifier = el
  MeGo <- enrichGO(meg, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL")
  Melk <- as.data.frame(MeGo)
  Mel <- c(Melk$geneID[1])
  MeVec_FULL <- append(MeVec_FULL, Mel)
  if (nrow(Melk) > 0)
  {
    MeVec_UP <- append(MeVec_UP, el)
    MeVec_UP <- append(MeVec_UP, Mel)
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
  if (nrow(Belk) > 0)
  {
    BeVec_UP <- append(BeVec_UP, el)
    BeVec_UP <- append(BeVec_UP, Bel)
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
  if (nrow(Celk) > 0)
  {
    CeVec_UP <- append(CeVec_UP, el)
    CeVec_UP <- append(CeVec_UP, Cel)
    print(cnetplot(CeGo, color_category='#1b9e77', 
                   color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component:", identifier)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
    if (length(BeGo$ID) >1)
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
#'*PATHWAY ANALYSIS FOR DOWN * 
################################################################################
