################################################################################
#This program accepts up and down lists of genes. Converts them to entrez IDs and submits them to Clue.io to find perturbagens.
#It is required to have an .Renvironment API key for Clue.io. This file can only consist of the API key as a single string.
################################################################################
library(cluequery)
library(org.Hs.eg.db)

#Read in dataframes
bdd <- read.csv("BeDF_DOWN_perGO.csv")
cdd <- read.csv("CeDF_DOWN_perGO.csv")
mdd <- read.csv("MeDF_DOWN_perGO.csv")

bdu <- read.csv("BeDF_UP_perGO.csv")
cdu <- read.csv("CeDF_UP_perGO.csv")
mdu <- read.csv("MeDF_UP_perGO.csv")

Samples <- c(
"T_11",
"T_8",
"T_12",
"C_7"
)
combos <- combn(Samples, 2, simplify = TRUE)

comboVec <- c()
#for column in combos create string of values and add to comboVec
for (i in 1:ncol(combos)){
  comboVec <- c(comboVec, paste(combos[,i], collapse = "."))
}

#make another vector of the same length as comboVec but with the opposite order of the samples
comboVec2 <- c()
for (i in 1:ncol(combos)){
  comboVec2 <- c(comboVec2, paste(rev(combos[,i]), collapse = "."))
}
#merge comboVec and comboVec2
comboVec <- c(comboVec, comboVec2)



#merge bdd,cdd and mdd into one dataframe
down <- rbind(bdd,cdd,mdd)
#merge bdu,cdu and mdu into one dataframe
up <- rbind(bdu,cdu,mdu)




Opptak <- c(
"phagocytosis",
"endosome",
"endocytosis",
"early endosome",
"endocytic vesicle",
"transport vesicle"
)

#Create a vector of the GO terms to be used
downOpptak <- down %>% select(contains(c(Opptak)))
upOpptak <- up %>% select(contains(c(Opptak)))

#
downOpptak <- downOpptak %>% select(contains(c(comboVec)))
upOpptakVs <-  upOpptak %>% select(contains(c(comboVec)))

Degradering <- c(
"proteolysis",
"proteasome-mediated ubiquitin-dependent protein catabolic process",
"ubiquitin-dependent protein catabolic process",
"metallopeptidase activity",
"lysosome",
"endolysosome",
"lysosomal protein catabolic process"
)

downDegradering <- down %>% select(contains(c(Degradering)))
upDegradering <- up %>% select(contains(c(Degradering)))

Inflammasjon <- c(
"inflammasome complex",
"neuroinflammatory response",
"inflammatory response"
)

downInflammasjon <- down %>% select(contains(c(Inflammasjon)))
upInflammasjon <- up %>% select(contains(c(Inflammasjon)))


UpList <- mapIds(org.Hs.eg.db, UpList, "ENTREZID", "SYMBOL")
UpList <- unique(UpList)

DownList <- mapIds(org.Hs.eg.db, DownList, "ENTREZID", "SYMBOL")
DownList <- unique(DownList)

pre_gmt <- clue_gmt_from_list(UpList, DownList, "Testing")

submission_result <- clue_query_submit(
    pre_gmt[["up"]], pre_gmt[["down"]], name="testNew"
)

subID <- submission_result$result$job_id

clue_query_wait(submission_result, interval = 120, timeout = 1200)