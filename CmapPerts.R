################################################################################
#This program accepts up and down lists of genes. Converts them to entrez IDs and submits them to Clue.io to find perturbagens.
#It is required to have an .Renvironment API key for Clue.io. This file can only consist of the API key as a single string.
################################################################################
library(cluequery)
library(org.Hs.eg.db)

UpList <- c('CD22', 'EPHA3', 'RAB4A', 'PICK1', 'ABCA7', 'ARC', 'CLIP3', 'LDLRAP1', 'VPS28', 'APPL2', 'ANKRD13B', 'ADRB2', 'LDLRAD4', 'CD4', 'CLCN3', 'CLCN4', 'ACKR2', 'GPER1', 'HSPD1', 'KCNQ1', 'KIFC1', 'PHB1', 'RAP1GAP', 'TLR4', 'DYSF', 'AP4M1', 'ZFYVE9', 'PTP4A3', 'INPP5F', 'GGA2', 'MLC1', 'ASTN2', 'BACE1', 'SNX5', 'TMEM230', 'NEURL1B', 'CLN6', 'PLEKHJ1', 'WLS', 'DNER', 'RNF157', 'PLD4', 'NAPEPLD', 'SNX30')
DownList <- c('B2M', 'CLN3', 'PTPN1', 'RAB5A', 'NUMB', 'SPHK1', 'RUBCN', 'ZFYVE16', 'RAB31', 'RAB21', 'RABGEF1', 'EHD4', 'ANKFY1', 'RAB17', 'RIN3', 'APBB2', 'APP', 'SERPINB1', 'RAPGEF1', 'HLA-A', 'HLA-B', 'HLA-DRA', 'HLA-E', 'MR1', 'LNPEP', 'MME', 'SLC11A2', 'NTRK1', 'PML', 'MAP2K1', 'PSEN1', 'RAB1A', 'RET', 'SH3GL1', 'SIAH2', 'SIGLEC1', 'SNX2', 'TLR3', 'UVRAG', 'VCAM1', 'PTP4A1', 'FZD5', 'PDLIM4', 'NRP1', 'RAB29', 'ATP6V0D1', 'HGS', 'LIPG', 'LITAF', 'TOM1', 'STX6', 'IFITM3', 'EHD1', 'RFTN1', 'FKBP15', 'STX12', 'SGK3', 'PARM1', 'RPS6KC1', 'VPS41', 'LAMP3', 'VPS4A', 'DBNL', 'CD274', 'APH1A', 'VAC14', 'GRIPAP1', 'PMEPA1', 'ACKR3', 'ARRDC3', 'WDFY1', 'WDFY4', 'GPR107', 'ZFYVE28', 'SNX27', 'MVB12B', 'SLC15A4', 'SNX20', 'DTX3L', 'RASGEF1B', 'SAMD9L', 'MARCHF8', 'FGD2', 'ATP6V0D2', 'WASHC2C', 'WASHC2A')


UpList <- mapIds(org.Hs.eg.db, UpList, "ENTREZID", "SYMBOL")
UpList <- unique(UpList)

DownList <- mapIds(org.Hs.eg.db, DownList, "ENTREZID", "SYMBOL")
DownList <- unique(DownList)

pre_gmt <- clue_gmt_from_list(UpList, DownList, "Testing")

submission_result <- clue_query_submit(
    pre_gmt[["up"]], pre_gmt[["down"]], name="test34R"
)

subID <- submission_result$result$job_id

clue_query_wait(submission_result, interval = 120, timeout = 1200)