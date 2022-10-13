
######Read libraries####
library(edgeR)
library(xlsx)
library(tidyverse)
library(clusterProfiler)
library("org.Hs.eg.db")

library(Rsamtools)  #another library needs part of this package to be unloaded. Always load Rsamtools first to avoid this
library(Rsubread)
library(DESeq2)
library(GenomicFeatures)
library("edgeR")
library("GenomicAlignments")
library("BiocParallel") #enable multicore
library(DESeq2)
library(ggplot2)
library(mygene)  #Using rentrez because of a 404 error
library(rentrez) #using mygene instead for now
#library("staplr")
library(qpdf)
library(biomaRt)


library(clusterProfiler)
library(enrichplot)

library(ggrepel)
########
#Point code to where files are:
setwd("/media/petear/SharedPart/RNAseq/Bam/")


######File for analyzos ######
bambamFiles <- list.files(pattern = "\\.bam$")

#GTF is of full
#gtfFile <- snakemake@input[["GTF"]]
gtfFile <- "/home/petear/Documents/Homo_sapiens.GRCh38.107.chr_patch_hapl_scaff.gtf"

########


#¤Note: make sure that the chromosome names of the genomic features in the annotation you use are consistent with the chromosome names of the reference used for read alignment. Otherwise, the scripts might fail to count any reads to features due to the mismatching names. For example, a common mistake is when the alignment files contain chromosome names in the style of 1 and the gene annotation in the style of chr1, or the other way around. See the seqlevelsStyle function in the GenomeInfoDb package for solutions. We can check the chromosome names (here called “seqnames”) in the alignment files like so:

#featCounter <- featureCounts(bambamFiles, annot.inbuilt = "hg38", isPairedEnd = TRUE, nthreads=18)
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




##############################################################
#######################PAIRWISE COMPARISONS####################


contrastDF <- data.frame(t(combn(unique(grVector),2)))
colnames(contrastDF) <- c("Y1", "Y2")

contrastDF <- contrastDF %>% unite("contCol", Y1:Y2, sep = "-", remove = FALSE)

contrastVec2 <- c(contrastDF$contCol)


contrasts <- makeContrasts(contrasts=contrastVec2, levels=design)

################################################################################
##For gener
#nedregulert i ‘Aβ vs TPA’       <- C7 vs C1  n1
#OG
#nedregulert i ‘LPS+Aβ vs Aβ’    <- T8 vs C7 n10
#OG
#oppregulert i ‘DHA+NS+Aβ vs Aβ’ <- T11 vs C7
################################################################################
SpecTestVec <- c("C_1-C_7","C_7-T_8","C_7-T_11")
SpecTestWantVec <- c("UP", "UP","DOWN")

SpecTestTippyTopGeneDFrr <- data.frame(matrix(ncol = 7, nrow = 0)) # create4 vector for list in loop
colnames(SpecTestTippyTopGeneDFrr) <- c("ENTREZID", "logFC", "logCPM", "F", "PValue", "FDR", "Versus")
SpecTestTippyTopGeneDFss <- SpecTestTippyTopGeneDFrr

#pdf("SpecTestEnrichPlotsALL.pdf", width = 20, height =20)


#The loop for adding all
for (versus in SpecTestVec)
{
  qlf <- glmQLFTest(fit, contrast=contrasts[,versus])
  SpecTestTippyTopGenes <- topTags(qlf, sort.by="PValue", n=Inf, p=0.05)
  SpecTestTippyTopGeneR2C <- tibble::rownames_to_column(data.frame(SpecTestTippyTopGenes), "ENTREZID")
  #assign(paste0("downReg", versus), listNegPosDFs[[versus]])
  #<- tippyTopTagsR2C[tippyTopTagsR2C$logFC < 0,]
  versVec <- rep(c(versus),times=nrow(SpecTestTippyTopGeneR2C))
  SpecTestTippyTopGeneR2C[["Versus"]] <- versVec
  if (nrow(SpecTestTippyTopGenes)> 0)
  {
    if (versus == "C_1-C-7" || versus == "C_7-T_8") #UP-REGULATED
    {
    df2<-subset(SpecTestTippyTopGeneR2C, logFC > 0)
    SpecTestTippyTopGeneDFrr <- rbind(SpecTestTippyTopGeneDFrr, df2)
    }
    if (versus == "C_7-T_11") #DOWN-REGULATED 
    {
      df2<-subset(SpecTestTippyTopGeneR2C, logFC < 0)
      SpecTestTippyTopGeneDFrr <- rbind(SpecTestTippyTopGeneDFrr, df2)
    }
  }
}

rr <- SpecTestTippyTopGeneDFrr$ENTREZID
counts <- table(rr)
multiplesrr <- names(counts)[counts > 1]
multiplesrr

#The loop for adding all INVERSE
#Can inverse be added to the full above?

for (versus in SpecTestVec)
{
  qlf <- glmQLFTest(fit, contrast=contrasts[,versus])
  SpecTestTippyTopGenes <- topTags(qlf, sort.by="PValue", n=Inf, p=0.05)
  SpecTestTippyTopGeneR2C <- tibble::rownames_to_column(data.frame(SpecTestTippyTopGenes), "ENTREZID")
  #assign(paste0("downReg", versus), listNegPosDFs[[versus]])
  #<- tippyTopTagsR2C[tippyTopTagsR2C$logFC < 0,]
  versVec <- rep(c(versus),times=nrow(SpecTestTippyTopGeneR2C))
  SpecTestTippyTopGeneR2C[["Versus"]] <- versVec
  if (nrow(SpecTestTippyTopGenes)> 0)
  {
    if (versus == "C_7-T_11") #DOWN-REGULATED
    {
      df2<-subset(SpecTestTippyTopGeneR2C, logFC > 0)
      SpecTestTippyTopGeneDFss <- rbind(SpecTestTippyTopGeneDFss, df2)
    }
    if (versus == "C_1-C-7" || versus == "C_7-T_8") #UP-REGULATED
    {
      df2<-subset(SpecTestTippyTopGeneR2C, logFC < 0)
      SpecTestTippyTopGeneDFss <- rbind(SpecTestTippyTopGeneDFss, df2)
    }
  }
}

ss <- SpecTestTippyTopGeneDFss$ENTREZID
counts <- table(ss)
multiplesss <- names(counts)[counts > 1]
multiplesss


#READ THE GENES TORMOD LIKES <3
################################################################################

wantedGo <- read.xlsx2("/media/petear/SharedPart/wantedGo.xlsx", sheetIndex = 1)
#wantedGo <- wantedGo[!(wantedGo$rogo=="GO:0044440"),] #Obsolete GO term - new term is already present. (commented out because removed while preparing data)
wantGoVec <- c(wantedGo$goID)



goTermDF <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(goTermDF) <- c("GOALL", "EVIDENCEALL", "ONTOLOGYALL", "ENTREZID")

#Make a DF of all genes connected with certain Go terms
for (term in wantGoVec)
{
  rtrv <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=term, columns="ENTREZID")
  goTermDF <- rbind(goTermDF, rtrv)
}

#
# ##############################################################################
# #Because i will sort on wanted GO terms I will set column ENTREZID as row names and tidy it on row names with tippyTopTags. But I will then need to remove duplicates (from 12976 to 5702 rows) ::
# #This SHOULD be okay because as long as one gene is present it will be participating in downstream analysis
#
# goTermDFnoDupes <- goTermDF[!duplicated(goTermDF$ENTREZID), ]
# rownames(goTermDFnoDupes) <- goTermDFnoDupes$ENTREZID
# drop <- "ENTREZID"
# goTermDFnoDupes = goTermDFnoDupes[,!(names(goTermDFnoDupes) %in% drop)]
#
#

########################################
#Find all genes


tippyTopGeneDF <- data.frame(matrix(ncol = 5, nrow = 0)) # create4 vector for list in loop
colnames(tippyTopGeneDF) <- c("logFC", "logCPM", "F", "PValue", "FDR")

geneDF <- data.frame(matrix(ncol = 10, nrow = 0)) # create4 vector for list in loop
colnames(geneDF) <- c("ENTREZID", "GOALL", "EVIDENCEALL", "ONTOLOGYALL", "logFC", "logCPM", "F","PValue", "FDR", "Versus"  )


pdf("VersaEnrichPlotsALL.pdf", width = 20, height =20)


#The loop for adding all
for (versus in contrastVec2)
{
  qlf <- glmQLFTest(fit, contrast=contrasts[,versus])
  tippyTopTags <- topTags(qlf, sort.by="PValue", n=Inf, p=0.05)
  tippyTopTagsR2C <- tibble::rownames_to_column(data.frame(tippyTopTags), "ENTREZID")
  #assign(paste0("downReg", versus), listNegPosDFs[[versus]])
  #<- tippyTopTagsR2C[tippyTopTagsR2C$logFC < 0,]
  merged <- merge(x=tippyTopTagsR2C, y=goTermDF,  by="ENTREZID", all.x=FALSE)
  versVec <- rep(c(versus),times=nrow(merged))
  merged[["Versus"]] = versVec
  #merged <- merged[order(merged$PValue),]
  geneDF <- rbind(geneDF, merged)
}


dev.off()
write.xlsx(geneDF, "/media/petear/SharedPart/TopGenesGOsorted.xlsx", rowNames = FALSE)

#################################################################
##Make plots on upreg and downreg genes
################################################################

##Find unique versuses:
uniques <- unique(geneDF$Versus)

##add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
geneDFUpDown <- geneDF
geneDFUpDown$DifferentiallyExpressed <- "NO"

geneDFUpDown$DifferentiallyExpressed[geneDFUpDown$logFC > 0.6 & geneDFUpDown$PValue < 0.05] <- "UP"
geneDFUpDown$DifferentiallyExpressed[geneDFUpDown$logFC < -0.6 & geneDFUpDown$PValue < 0.05] <- "DOWN"


## Now write down the name of genes beside the points...
## Create a new column "delabel" to dataframe, that will contain the name of genes differentially expressed (NA in case they are not)
geneDFUpDown$DFlabel <- NA
geneDFUpDown$DFlabel[geneDFUpDown$DifferentiallyExpressed != "NO"] <- geneDFUpDown$ENTREZID[geneDFUpDown$DifferentiallyExpressed != "NO"]


pdf("volcano.pdf", width = 20, height =20)
for (versa in uniques)
{
dfSubset <- geneDFUpDown[geneDFUpDown$Versus == versa,]
print(ggplot(data=dfSubset, aes(x=logFC, y=-log10(PValue), col=DifferentiallyExpressed, label=DFlabel)) +
        geom_point() +
        geom_text_repel() +
        scale_color_manual(values=c("#1b9e77", "black", "#d95f02")) +
        geom_vline(xintercept=c(-0.6, 0.6), col="#7570b3") +
        geom_hline(yintercept=-log10(0.05), col="#7570b3") +
        ggtitle("Volcanoplot of", versa))
}
dev.off()



pdf("volcano.pdf", width = 20, height =20)
for (versa in uniques)
{
dfSubset <- geneDF[geneDF$Versus == versa,]
print(ggplot(data=dfSubset, aes(x=logFC, y=PValue)) + geom_point())

}
dev.off()

###################################################################
##Make plots on gene Ontology and kegg based on gene ID only
##################################################################

##Remove ENTREZID duplicates
geneDFuniq <- geneDF[!duplicated(geneDF$ENTREZID),]



##subset geneDF based on versus (and maybe then subset the subset based on Ontology



pdf("VersaEnrichPlotsVvV.pdf", width = 20, height =20)
  
##Dont need to take Versus (or anything else) into account because its only based on the gene name itself
for (vers in uniques) 
{
  tryCatch(
  {
  versSubset <- geneDFuniq[geneDFuniq$Versus == vers,]  
  pdfName <- paste(vers,".pdf", sep="")
  pdf(pdfName, width = 20, height =20)
  print(vers)
  MeGo <- enrichGO(versSubset$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "ENTREZID")
  print(cnetplot(MeGo, color_category='#1b9e77',
         color_gene='#d95f02') + ggtitle(paste("Ontology for Molecular Function",vers)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  print(goplot(MeGo) + ggtitle(paste("Ontology for Molecular Function of", vers)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))

  BeGo <- enrichGO(versSubset$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID")
  print(cnetplot(BeGo, color_category='#1b9e77',
         color_gene='#d95f02') + ggtitle(paste("Ontology for Biological Process",vers)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  print(goplot(BeGo) + ggtitle(paste("Ontology for Biological Process", vers)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))

  CeGo <- enrichGO(versSubset$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "ENTREZID")
  print(cnetplot(CeGo, color_category='#1b9e77', 
         color_gene='#d95f02') + ggtitle(paste("Ontology for Cellular Component", vers)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  print(goplot(CeGo) + ggtitle(paste("Ontology for Cellular Component", vers)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))

  AeGo <- enrichGO(versSubset$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "ENTREZID")
  print(cnetplot(AeGo, color_category='#1b9e77', 
         color_gene='#d95f02') + ggtitle(paste("Ontology for all three GO classificiations", vers)) + theme(plot.margin=unit(c(0.0,0.4,0.0,0.4), 'cm')))
  dev.off()
  }, error=function(e){dev.off()})
}

dev.off()

pdf("KEGG.pdf")
eKEGG <- enrichKEGG(reactomeVec)
  cnetplot(eKEGG, color_category='#1b9e77', 
  color_gene='#d95f02')  + ggtitle("KEGG") + theme(plot.margin=unit(c(0.0,0.2,0.0,0.2), 'cm'))
  cnetplot(eKEGG, showCategory = 22, color_category='#1b9e77', 
         color_gene='#d95f02') + ggtitle("KEGG -forced extended connections") + theme(plot.margin=unit(c(0.0,0.2,0.0,0.2), 'cm'))
  
dev.off()

