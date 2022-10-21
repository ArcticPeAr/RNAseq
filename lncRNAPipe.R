library(rjson)
library(biomaRt)
RNAcentralJson <- fromJSON(file = "/media/petear/SharedPart/mmp_AND_so_rna_type_nameLncRNA_AND_TAXONOMY9606_AND_entry_typeSequence.json")

MMPidVec <- c()

n <- length(RNAcentralJson)

#Collecting all ids of MMP found in RNAcentral.org 
for (var in 1:n)
{
  MMPidVec <- append(MMPidVec, RNAcentralJson[[var]]$id)
}


################################################################################
#'*Convert from RNACentral to ENSEMBL*
################################################################################
marty <- useMart(biomart = "rnacentral",
                 dataset = "hsapiens_gene_ensembl",
                 host = "https://www.ensembl.org")

genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","entrezgene_id"),
               values = n1_3JoinSorted$ENSEMBL, 
               mart = marty)