library(rjson)

RNAcentralJson <- fromJSON(file = "/media/petear/SharedPart/mmp_AND_so_rna_type_nameLncRNA_AND_TAXONOMY9606_AND_entry_typeSequence.json")

MMPidVec <- c()

n <- length(RNAcentralJson)

#Collecting all ids  
for (var in 1:n)
{
  MMPidVec <- append(idVec, RNAcentralJson[[var]]$id)
}
