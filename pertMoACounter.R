library(readxl)
library(writexl)
library(openxlsx)
library(reshape2)
library(tidyverse)

###########################################################################ALL
##########################################################################

pertDFAll <- read_excel("/home/petear/MEGA/TormodGroup/InputData/AllAlgoSi09Aug.xlsx", col_names = FALSE, sheet = 1)
pertDFAll <- as.data.frame(pertDFAll)
#remove column 2-5 and 7-10
pertDFAll <- pertDFAll[,c(1,6)]

#remove duplicates in ...1
pertDFAll <- pertDFAll %>%
  distinct(...1, .keep_all = TRUE)

# 1. Generate the count dataframe
countDFAll <- as.data.frame(table(pertDFAll$"...6"))
colnames(countDFAll) <- c("...6", "n")

# Sort the dataframe by count in descending order
countDFAll <- countDFAll[order(-countDFAll$n), ]

# Reshape the dataframe
wideDFAll <- pertDFAll %>%
  group_by(...6) %>%
  mutate(row = row_number()) %>%
  spread(key = ...6, value = ...1) %>%
  select(-row)


###########################################################################ECYT
##########################################################################

pertDFEcyt <- read_excel("/home/petear/MEGA/TormodGroup/InputData/EndocyticAlgoSi09Aug.xlsx", col_names = FALSE, sheet = 1)
pertDFEcyt <- as.data.frame(pertDFEcyt)
#remove column 2-5 and 7-10
pertDFEcyt <- pertDFEcyt[,c(1,6)]

#remove duplicates in ...1
pertDFEcyt <- pertDFEcyt %>%
  distinct(...1, .keep_all = TRUE)

# 1. Generate the count dataframe
countDFEcyt <- as.data.frame(table(pertDFEcyt$"...6"))
colnames(countDFEcyt) <- c("...6", "n")

# Sort the dataframe by count in descending order
countDFEcyt <- countDFEcyt[order(-countDFEcyt$n), ]

# Reshape the dataframe
wideDFEcyt <- pertDFEcyt %>%
  group_by(...6) %>%
  mutate(row = row_number()) %>%
  spread(key = ...6, value = ...1) %>%
  select(-row)

###########################################################################ELYS
##########################################################################

pertDFElys <- read_excel("/home/petear/MEGA/TormodGroup/InputData/endolysAlgoSi09Aug.xlsx", col_names = FALSE, sheet = 1)
pertDFElys <- as.data.frame(pertDFElys)
#remove column 2-5 and 7-10
pertDFElys <- pertDFElys[,c(1,6)]

#remove duplicates in ...1
pertDFElys <- pertDFElys %>%
  distinct(...1, .keep_all = TRUE)

# 1. Generate the count dataframe
countDFElys <- as.data.frame(table(pertDFElys$"...6"))
colnames(countDFElys) <- c("...6", "n")

# Sort the dataframe by count in descending order
countDFElys <- countDFElys[order(-countDFElys$n), ]

# Reshape the dataframe
wideDFElys <- pertDFElys %>%
  group_by(...6) %>%
  mutate(row = row_number()) %>%
  spread(key = ...6, value = ...1) %>%
  select(-row)


###########################################################################Print
##########################################################################

write.xlsx(list("CountsAll" = countDFAll, "PertsPrMoAAll" = wideDFAll, "CountsEcyt" = countDFEcyt, "PertsPrMoAEcyt" = wideDFEcyt, "CountsElys" = countDFElys, "PertsPrMoAElys" = wideDFElys), "MoACounts.xlsx")
