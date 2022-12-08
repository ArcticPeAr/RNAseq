################################################################################
#This program accepts up and down lists of genes. Converts them to entrez IDs and submits them to Clue.io to find perturbagens.
#It is required to have an .Renvironment API key for Clue.io. This file can only consist of the API key as a single string.
################################################################################
library(cluequery)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
#Read in dataframes
bdd <- read.csv("BeDF_DOWN_perGO.csv")
cdd <- read.csv("CeDF_DOWN_perGO.csv")
mdd <- read.csv("MeDF_DOWN_perGO.csv")

bdu <- read.csv("BeDF_UP_perGO.csv")
cdu <- read.csv("CeDF_UP_perGO.csv")
mdu <- read.csv("MeDF_UP_perGO.csv")

#merge bdd,cdd and mdd into one dataframe
down <- rbind(bdd,cdd,mdd)
#merge bdu,cdu and mdu into one dataframe
up <- rbind(bdu,cdu,mdu)

down[is.na(down)] <- ""
up[is.na(up)] <- ""

#Read in TippyTopGeneDF_ALL
TippyTopGeneDF_ALL <- read.csv("TippyTopGeneDF_ALL.csv")

################################################################################
#CREATE VECTORS OF VERSUSES
#ALL vs C7
#LOWEST VS HIGHEST OPPTAK
#LOWEST VS HIGHEST DEGRADETION
#LOWEST VS HIGHEST CLEARANCE
################################################################################
#Special cases to be investigated
OpptakVec <- c("T_11.T_12")
OpptakVec2 <- c("T_12.T_11")
#DegradationVec <- c("T_11","T_12")  #Same as OpptakVec
#ClearanceVec <- c("T_11","T_12")   #Same as OpptakVec


VersVec <- c(OpptakVec,OpptakVec2)

################################################################################
#CREATE VECTORS OF VERSUSES
################################################################################
#General cases to be investigated
TestVec <- c("C_7")

Samples <- c(
"T_5",
"T_8",
"T_12",
"T_11"
)

#Create a vector of all possible combinations of TestVec and Samples
RunVec <- c()
for (i in 1:length(TestVec)){
for (j in 1:length(Samples)){
RunVec <- c(RunVec, paste(TestVec[i], Samples[j], sep = "."))
}
}
#make another vector of the same length as comboVec but with the opposite order of the samples
RunVec2 <- c()
for (i in 1:length(TestVec)){
for (j in 1:length(Samples)){
RunVec2 <- c(RunVec2, paste(Samples[j], TestVec[i], sep = "."))
}
}

RunVec <- c(RunVec, RunVec2)


################################################################################
#THIS IS FOR ALL COMBINATIONS OF TESTVEC AND SAMPLE
#Not needed for this analysis
################################################################################
# combos <- combn(Samples, 2, simplify = TRUE)

# comboVec <- c()
# #for column in combos create string of values and add to comboVec
# for (i in 1:ncol(combos)){
#   comboVec <- c(comboVec, paste(combos[,i], collapse = "."))
# }

# #make another vector of the same length as comboVec but with the opposite order of the samples
# comboVec2 <- c()
# for (i in 1:ncol(combos)){
#   comboVec2 <- c(comboVec2, paste(rev(combos[,i]), collapse = "."))
# }
# #merge comboVec and comboVec2
# comboVec <- c(comboVec, comboVec2)


################################################################################
#Load GO terms
################################################################################
anaVec <- c(RunVec, VersVec)

Opptak <- c(
"phagocytosis",
"endosome",
"endocytosis",
"early.endosome",
"endocytic.vesicle",
"transport.vesicle"
)
Degradering <- c(
"proteolysis",
"proteasome-mediated.ubiquitin-dependent.protein.catabolic.process",
"ubiquitin-dependent.protein.catabolic.process",
"metallopeptidase.activity",
"lysosome",
"endolysosome",
"lysosomal.protein.catabolic.process"
)
Inflammasjon <- c(
"inflammasome.complex",
"neuroinflammatory.response",
"inflammatory.response",
"autophagy"
)

Clearance <- c(Opptak, Degradering)
#Make dataframes smaller by removing columns with no relation to the terms:

################################################################################
#OPPTAK ClueList creation 
################################################################################
downOpptak <- down %>% select(contains(c(Opptak)))
upOpptak <- up %>% select(contains(c(Opptak)))


OpptakUp <- c()
clueList_Opptak_Up <- list()

OpptakDown <- c()
clueList_Opptak_Down <- list()

# For every value in anavec, find the up genes for opptak
for (i in 1:length(anaVec)){
    OpptakUp <- c()
    OpptakDown <- c()
    downOpptakAna <- downOpptak %>% select(contains(c(anaVec[i])))
    upOpptakAna <- upOpptak %>% select(contains(c(anaVec[i])))

    for (column in colnames(downOpptakAna)){
        for (row in 1:nrow(downOpptakAna)){
            if (downOpptakAna[row, column] != ""){
                OpptakDown <- c(OpptakDown, downOpptakAna[row, column])
            }
        }
    OpptakDown <- unique(OpptakDown)
    #add "Down" to list item
    clueList_Opptak_Down[[anaVec[i]]] <- OpptakDown
    
    }

    for (column in colnames(upOpptakAna)){
        for (row in 1:nrow(upOpptakAna)){
            if (upOpptakAna[row, column] != ""){
                OpptakUp <- c(OpptakUp, upOpptakAna[row, column])
            }        
        }
    OpptakUp <- unique(OpptakUp)
    #add "Up" to list item
    clueList_Opptak_Up[[anaVec[i]]] <- OpptakUp
    }
}


################################################################################
#DEGRADERING ClueList creation 
################################################################################
downDegradering <- down %>% select(contains(c(Degradering)))
upDegradering <- up %>% select(contains(c(Degradering)))


DegraderingUp <- c()
clueList_Degradering_Up <- list()

DegraderingDown <- c()
clueList_Degradering_Down <- list()

# For every value in anavec, find the up genes for opptak
for (i in 1:length(anaVec)){
    DegraderingUp <- c()
    DegraderingDown <- c()
    downDegraderingAna <- downDegradering %>% select(contains(c(anaVec[i])))
    upDegraderingAna <- upDegradering %>% select(contains(c(anaVec[i])))

    for (column in colnames(downDegraderingAna)){
        for (row in 1:nrow(downDegraderingAna)){
            if (downDegraderingAna[row, column] != ""){
                DegraderingDown <- c(DegraderingDown, downDegraderingAna[row, column])
            }
        }
    DegraderingDown <- unique(DegraderingDown)
    #add "Down" to list item
    clueList_Degradering_Down[[anaVec[i]]] <- DegraderingDown
    
    }

    for (column in colnames(upDegraderingAna)){
        for (row in 1:nrow(upDegraderingAna)){
            if (upDegraderingAna[row, column] != ""){
                DegraderingUp <- c(DegraderingUp, upDegraderingAna[row, column])
            }        
        }
    DegraderingUp <- unique(DegraderingUp)
    #add "Up" to list item
    clueList_Degradering_Up[[anaVec[i]]] <- DegraderingUp
    }
}

################################################################################
#INFLAMMASJON ClueList creation 
################################################################################

downInflammasjon <- down %>% select(contains(c(Inflammasjon)))
upInflammasjon <- up %>% select(contains(c(Inflammasjon)))

downInflammasjon <- down %>% select(contains(c(Inflammasjon)))
upInflammasjon <- up %>% select(contains(c(Inflammasjon)))


InflammasjonUp <- c()
clueList_Inflammasjon_Up <- list()

InflammasjonDown <- c()
clueList_Inflammasjon_Down <- list()

# For every value in anavec, find the up genes for opptak
for (i in 1:length(anaVec)){
    InflammasjonUp <- c()
    InflammasjonDown <- c()
    downInflammasjonAna <- downInflammasjon %>% select(contains(c(anaVec[i])))
    upInflammasjonAna <- upInflammasjon %>% select(contains(c(anaVec[i])))

    for (column in colnames(downInflammasjonAna)){
        for (row in 1:nrow(downInflammasjonAna)){
            if (downInflammasjonAna[row, column] != ""){
                InflammasjonDown <- c(InflammasjonDown, downInflammasjonAna[row, column])
            }
        }
    InflammasjonDown <- unique(InflammasjonDown)
    #add "Down" to list item
    clueList_Inflammasjon_Down[[anaVec[i]]] <- InflammasjonDown
    
    }

    for (column in colnames(upInflammasjonAna)){
        for (row in 1:nrow(upInflammasjonAna)){
            if (upInflammasjonAna[row, column] != ""){
                InflammasjonUp <- c(InflammasjonUp, upInflammasjonAna[row, column])
            }        
        }
    InflammasjonUp <- unique(InflammasjonUp)
    #add "Up" to list item
    clueList_Inflammasjon_Up[[anaVec[i]]] <- InflammasjonUp
    
    }
}

################################################################################
#CLEARANCE ClueList creation
################################################################################
downClearance <- down %>% select(contains(c(Clearance)))
upClearance <- up %>% select(contains(c(Clearance)))

ClearanceUp <- c()
clueList_Clearance_Up <- list()

ClearanceDown <- c()
clueList_Clearance_Down <- list()

# For every value in anavec, find the up genes for opptak
for (i in 1:length(anaVec)){
    ClearanceUp <- c()
    ClearanceDown <- c()
    downClearanceAna <- downClearance %>% select(contains(c(anaVec[i])))
    upClearanceAna <- upClearance %>% select(contains(c(anaVec[i])))

    for (column in colnames(downClearanceAna)){
        for (row in 1:nrow(downClearanceAna)){
            if (downClearanceAna[row, column] != ""){
                ClearanceDown <- c(ClearanceDown, downClearanceAna[row, column])
            }
        }
    ClearanceDown <- unique(ClearanceDown)
    #add "Down" to list item
    clueList_Clearance_Down[[anaVec[i]]] <- ClearanceDown
    
    }

    for (column in colnames(upClearanceAna)){
        for (row in 1:nrow(upClearanceAna)){
            if (upClearanceAna[row, column] != ""){
                ClearanceUp <- c(ClearanceUp, upClearanceAna[row, column])
            }        
        }
    ClearanceUp <- unique(ClearanceUp)
    #add "Up" to list item
    clueList_Clearance_Up[[anaVec[i]]] <- ClearanceUp
    
    }
}

#downOpptakVS <- downOpptak %>% select(contains(c(anaVec)))
#upOpptakVS <-  upOpptak %>% select(contains(c(anaVec)))

# Loop over values in anaVec and find columns in downOpptak and upOpptak that contain the value in anaVec
# Then create vectors of the values in 


################################################################################
#GET GO'ING
################################################################################
#The lists are
#clueList_Opptak_Up
#clueList_Opptak_Down
#ClueList_Degradering_Up
#ClueList_Degradering_Down
#ClueList_Inflammasjon_Up
#ClueList_Inflammasjon_Down


################################################################################
#SENDING OPPTAK TO CLUE
################################################################################
#get list names from clueList_Opptak_Up
listNamesOpptakUP <- names(clueList_Opptak_Up)
listNamesOpptakDOWN <- names(clueList_Opptak_Down)
if (length(setdiff(listNamesOpptakUP, listNamesOpptakDOWN)>0)){
    print("OPPTAKLISTENE are not the same")
    stop()
}

for (i in 1:length(listNamesOpptakUP)){
    tryCatch({
    namestringOpptakF <- paste(listNamesOpptakUP[i], "Opptak", sep = "_")
    namestringOpptakR <- paste(listNamesOpptakUP[i], "Opptak_Reversed", sep = "_")
    nameString <- gsub(".",  "-", listNamesOpptakUP[i], fixed = TRUE)
    SpecTT <- TippyTopGeneDF_ALL %>% filter(Versus == nameString)
    SpecTT$ENTREZID <- as.character(SpecTT$ENTREZID)
    #
    UpList <- clueList_Opptak_Up[[i]]
    UpListDF <- as.data.frame(mapIds(org.Hs.eg.db, UpList, "ENTREZID", "SYMBOL"))
    colnames(UpListDF) <- "ENTREZID"
    UpListDF$Symbol <- rownames(UpListDF)
    UpJoined <- inner_join(UpListDF, SpecTT, by="ENTREZID")
    UpJoined <- UpJoined[order(UpJoined$logFC, decreasing=TRUE),]
    #if (nrow(UpJoined) > 150){
    #    UpJoined <- UpJoined[1:150,]
    #}
    DownList <- clueList_Opptak_Down[[i]]
    DownListDF <- as.data.frame(mapIds(org.Hs.eg.db, DownList, "ENTREZID", "SYMBOL"))
    colnames(DownListDF) <- "ENTREZID"
    DownListDF$Symbol <- rownames(DownListDF)
    DownJoined <- inner_join(DownListDF, SpecTT, by="ENTREZID")
    DownJoined <- DownJoined[order(DownJoined$logFC, decreasing=FALSE),]
    #if (nrow(DownJoined) > 150){
    #    DownJoined <- DownJoined[1:150,]
    #}
    pre_gmtF <- clue_gmt_from_list(UpJoined$ENTREZID, DownJoined$ENTREZID, namestringOpptakF)
    submission_result <- clue_query_submit(
        pre_gmtF[["up"]], pre_gmt[["down"]], name=namestringOpptakF
    )
    pre_gmtR <- clue_gmt_from_list(DownJoined$ENTREZID, UpJoined$ENTREZID, namestringOpptakR)
    submission_result <- clue_query_submit(
        pre_gmtR[["up"]], pre_gmtR[["down"]], name=namestringOpptakR
    )
}, error = function(e) {
    print(e)
})
}


################################################################################
#SENDING DEGRADERING TO CLUE
################################################################################
listNamesDegraderingUP <- names(clueList_Degradering_Up)
listNamesDegraderingDOWN <- names(clueList_Degradering_Down)
if (length(setdiff(listNamesDegraderingUP, listNamesDegraderingDOWN)>0)){
    print("DEGRADERINGLISTENE are not the same")
    stop()
}

for (i in 1:length(listNamesDegraderingUP)){
    tryCatch({
    namestringDegraderingF <- paste(listNamesDegraderingUP[i], "Degradering", sep = "_")
    namestringDegraderingR <- paste(listNamesDegraderingUP[i], "Degradering_Reversed", sep = "_")
    nameString <- gsub(".",  "-", listNamesDegraderingUP[i], fixed = TRUE)
    SpecTT <- TippyTopGeneDF_ALL %>% filter(Versus == nameString)
    SpecTT$ENTREZID <- as.character(SpecTT$ENTREZID)
    #
    UpList <- clueList_Degradering_Up[[i]]
    UpListDF <- as.data.frame(mapIds(org.Hs.eg.db, UpList, "ENTREZID", "SYMBOL"))
    colnames(UpListDF) <- "ENTREZID"
    UpListDF$Symbol <- rownames(UpListDF)
    UpJoined <- inner_join(UpListDF, SpecTT, by="ENTREZID")
    UpJoined <- UpJoined[order(UpJoined$logFC, decreasing=TRUE),]
    #if (nrow(UpJoined) > 150){
    #    UpJoined <- UpJoined[1:150,]
    #}
    DownList <- clueList_Degradering_Down[[i]]
    DownListDF <- as.data.frame(mapIds(org.Hs.eg.db, DownList, "ENTREZID", "SYMBOL"))
    colnames(DownListDF) <- "ENTREZID"
    DownListDF$Symbol <- rownames(DownListDF)
    DownJoined <- inner_join(DownListDF, SpecTT, by="ENTREZID")
    DownJoined <- DownJoined[order(DownJoined$logFC, decreasing=FALSE),]
    #if (nrow(DownJoined) > 150){
    #    DownJoined <- DownJoined[1:150,]
    #}
    pre_gmt <- clue_gmt_from_list(UpJoined$ENTREZID, DownJoined$ENTREZID, namestringDegraderingF)
    submission_result <- clue_query_submit(
        pre_gmt[["up"]], pre_gmt[["down"]], name=namestringDegraderingF
    )
    pre_gmtR <- clue_gmt_from_list(DownJoined$ENTREZID, UpJoined$ENTREZID, namestringDegraderingR)
    submission_result <- clue_query_submit(
        pre_gmtR[["up"]], pre_gmtR[["down"]], name=namestringDegraderingR
    )
}, error = function(e) {
    print(e)
})
}


################################################################################
#SENDING INFLAMMASJON TO CLUE
################################################################################
listNamesInflammasjonUP <- names(clueList_Inflammasjon_Up)
listNamesInflammasjonDOWN <- names(clueList_Inflammasjon_Down)
if (length(setdiff(listNamesInflammasjonUP, listNamesInflammasjonDOWN)>0)){
    print("INFLAMMASJONLISTENE are not the same")
    stop()
}

for (i in 1:length(listNamesInflammasjonUP)){
    tryCatch({
    namestringInflammasjonF <- paste(listNamesInflammasjonUP[i], "Inflammasjon", sep = "_")
    namestringInflammasjonR <- paste(listNamesInflammasjonUP[i], "Inflammasjon_Reversed", sep = "_")
    nameString <- gsub(".",  "-", listNamesInflammasjonUP[i], fixed = TRUE)
    SpecTT <- TippyTopGeneDF_ALL %>% filter(Versus == nameString)
    SpecTT$ENTREZID <- as.character(SpecTT$ENTREZID)
    #
    UpList <- clueList_Inflammasjon_Up[[i]]
    UpListDF <- as.data.frame(mapIds(org.Hs.eg.db, UpList, "ENTREZID", "SYMBOL"))
    colnames(UpListDF) <- "ENTREZID"
    UpListDF$Symbol <- rownames(UpListDF)
    UpJoined <- inner_join(UpListDF, SpecTT, by="ENTREZID")
    UpJoined <- UpJoined[order(UpJoined$logFC, decreasing=TRUE),]
    #if (nrow(UpJoined) > 150){
    #    UpJoined <- UpJoined[1:150,]
    #}
    DownList <- clueList_Inflammasjon_Down[[i]]
    DownListDF <- as.data.frame(mapIds(org.Hs.eg.db, DownList, "ENTREZID", "SYMBOL"))
    colnames(DownListDF) <- "ENTREZID"
    DownListDF$Symbol <- rownames(DownListDF)
    DownJoined <- inner_join(DownListDF, SpecTT, by="ENTREZID")
    DownJoined <- DownJoined[order(DownJoined$logFC, decreasing=FALSE),]
    #if (nrow(DownJoined) > 150){
    #    DownJoined <- DownJoined[1:150,]
    #}
    pre_gmt <- clue_gmt_from_list(UpJoined$ENTREZID, DownJoined$ENTREZID, namestringInflammasjonF)
    submission_result <- clue_query_submit(
        pre_gmt[["up"]], pre_gmt[["down"]], name=namestringInflammasjonF
    )
    pre_gmtR <- clue_gmt_from_list(DownJoined$ENTREZID, UpJoined$ENTREZID, namestringInflammasjonR)
    submission_result <- clue_query_submit(
        pre_gmtR[["up"]], pre_gmtR[["down"]], name=namestringInflammasjonR
    )
}, error = function(e) {
    print(e)
})
}

################################################################################
#SENDING CLEARANCE TO CLUE
################################################################################
listNamesClearanceUP <- names(clueList_Clearance_Up)
listNamesClearanceDOWN <- names(clueList_Clearance_Down)
if (length(setdiff(listNamesClearanceUP, listNamesClearanceDOWN)>0)){
    print("CLEARANCELISTENE are not the same")
    stop()
}

for (i in 1:length(listNamesClearanceUP)){
    tryCatch({
    namestringClearanceF <- paste(listNamesClearanceUP[i], "Clearance", sep = "_")
    namestringClearanceR <- paste(listNamesClearanceUP[i], "Clearance_Reversed", sep = "_")
    nameString <- gsub(".",  "-", listNamesClearanceUP[i], fixed = TRUE)
    SpecTT <- TippyTopGeneDF_ALL %>% filter(Versus == nameString)
    SpecTT$ENTREZID <- as.character(SpecTT$ENTREZID)
    #
    UpList <- clueList_Clearance_Up[[i]]
    UpListDF <- as.data.frame(mapIds(org.Hs.eg.db, UpList, "ENTREZID", "SYMBOL"))
    colnames(UpListDF) <- "ENTREZID"
    UpListDF$Symbol <- rownames(UpListDF)
    UpJoined <- inner_join(UpListDF, SpecTT, by="ENTREZID")
    UpJoined <- UpJoined[order(UpJoined$logFC, decreasing=TRUE),]
    #if (nrow(UpJoined) > 150){
    #    UpJoined <- UpJoined[1:150,]
    #}
    DownList <- clueList_Clearance_Down[[i]]
    DownListDF <- as.data.frame(mapIds(org.Hs.eg.db, DownList, "ENTREZID", "SYMBOL"))
    colnames(DownListDF) <- "ENTREZID"
    DownListDF$Symbol <- rownames(DownListDF)
    DownJoined <- inner_join(DownListDF, SpecTT, by="ENTREZID")
    DownJoined <- DownJoined[order(DownJoined$logFC, decreasing=FALSE),]
    #if (nrow(DownJoined) > 150){
    #    DownJoined <- DownJoined[1:150,]
    #}
    pre_gmt <- clue_gmt_from_list(UpJoined$ENTREZID, DownJoined$ENTREZID, namestringClearanceF)
    submission_result <- clue_query_submit(
        pre_gmt[["up"]], pre_gmt[["down"]], name=namestringClearanceF
    )
    pre_gmtR <- clue_gmt_from_list(DownJoined$ENTREZID, UpJoined$ENTREZID, namestringClearanceR)
    submission_result <- clue_query_submit(
        pre_gmtR[["up"]], pre_gmtR[["down"]], name=namestringClearanceR
    )
}, error = function(e) {
    print(e)
})
}



#Get the results from CLUE
clue_query_wait(submission_result, interval = 120, timeout = 1200)

