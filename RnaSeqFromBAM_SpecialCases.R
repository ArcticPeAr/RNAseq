library(Rsubread)
library(DESeq2)
library(ggplot2)
library(openxlsx)
gtfFile <- "/media/petear/7BFE-B23C/RNAseq/Homo_sapiens.GRCh38.111.chr_patch_hapl_scaff.gtf"

# Specify the directory containing BAM files
bamFilesDirectory <- "/media/petear/7BFE-B23C/RNAseq/Bam/"

# List all .bam files with specific prefixes
samplePrefixes <- c("T_9", "T_10", "T_11")

pattern <- paste0("^(", paste(samplePrefixes, collapse = "|"), ").*\\.bam$")

bamFiles <- list.files(path = bamFilesDirectory, pattern = pattern, full.names = TRUE)

fc <- featureCounts(files = bamFiles, isGTFAnnotationFile=TRUE, annot.ext = gtfFile, isPairedEnd=TRUE, nthreads=15)
#fc <- readRDS("/home/petear/MEGA/TormodGroup/InputData/fc.rds")

counts <- fc$counts

condition <- substr(colnames(counts), 1, 4)
condition <- as.factor(condition)

sampleInfo <- data.frame(
    SampleName = colnames(counts),
    Condition = condition
)

sampleInfo$group <- ifelse(sampleInfo$Condition %in% c("T_9", "T_10"), "T_9_T_10", "T_11")

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleInfo, design = ~ group)
dds <- DESeq(dds)

res_T11_vs_T9_T10 <- results(dds, contrast = c("group", "T_11", "T_9_T_10"))

res_T11_vs_T9_T10 <- res_T11_vs_T9_T10[!(res_T11_vs_T9_T10[,1] == 0 & rowSums(is.na(res_T11_vs_T9_T10[, -1])) == (ncol(res_T11_vs_T9_T10) - 1)), ]

res_T11_vs_T9_T10 <- res_T11_vs_T9_T10[order(res_T11_vs_T9_T10$pvalue), ]



write.xlsx(res_T11_vs_T9_T10, "/home/petear/MEGA/TormodGroup/InputData/T11vsT9andT10_dseq2.xlsx")
