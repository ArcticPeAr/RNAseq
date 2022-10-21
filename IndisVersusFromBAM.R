
######Read libraries####
library(edgeR)
library(xlsx)
library(tidyverse)
library(clusterProfiler)
library("org.Hs.eg.db")

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
  TippyTopGeneR2C <- tibble::rownames_to_column(data.frame(SpecTestTippyTopGenes), "ENTREZID")
  #assign(paste0("downReg", versus), listNegPosDFs[[versus]])
  #<- tippyTopTagsR2C[tippyTopTagsR2C$logFC < 0,]
  versVec <- rep(c(versus),times=nrow(TippyTopGeneR2C))
  TippyTopGeneR2C[["Versus"]] <- versVec
  if (nrow(SpecTestTippyTopGenes)> 0)
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
  TippyTopGeneR2C <- tibble::rownames_to_column(data.frame(SpecTestTippyTopGenes), "ENTREZID")
  #assign(paste0("downReg", versus), listNegPosDFs[[versus]])
  #<- tippyTopTagsR2C[tippyTopTagsR2C$logFC < 0,]
  versVec <- rep(c(versus),times=nrow(TippyTopGeneR2C))
  TippyTopGeneR2C[["Versus"]] <- versVec
  if (nrow(SpecTestTippyTopGenes)> 0)
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
  TippyTopGeneR2C <- tibble::rownames_to_column(data.frame(SpecTestTippyTopGenes), "ENTREZID")
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

uniqGOTermVec <- unique(goTermDF$GOALL)
uniqVers_ALL <- unique(TippyTopGeneDF_ALL$Versus) #Done again because not all versuses might get hits

joinedDF_GO_Vers <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(joinedDF_GO_Vers) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENTREZID", "ENSEMBL", "logFC", "logCPM", "F", "PValue", "FDR", "Versus")

for (goterm in uniqGOTermVec)
{
  dfGO <- goTermDF[goTermDF$GOALL==goterm,]
  for (vers in uniqVers_ALL)
  {
    dfVERS <- TippyTopGeneDF_ALL[TippyTopGeneDF_ALL$Versus==vers,]
    GoVersJoined <- inner_join(dfGO, dfVERS, by = "ENTREZID")
    joinedDF_GO_Vers <- rbind(joinedDF_GO_Vers, GoVersJoined)
  }
}

