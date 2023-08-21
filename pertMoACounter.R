library(readxl)
library(writexl)
library(dplyr)

pertDF <- read_excel("/home/petear/MEGA/TormodGroup/InputData/AllAlgoSi09Aug.xlsx", col_names = FALSE, sheet = 1)
pertDF <- as.data.frame(pertDF)
#remove column 2-5 and 7-10
pertDF <- pertDF[,c(1,6)]

#remove duplicates in ...1
pertDF <- pertDF %>%
  distinct(...1, .keep_all = TRUE)

#count each string in ...6
pertTup <- pertDF %>% group_by(...6) %>% summarise(n = n())
pertTup <- pertTup[order(pertTup$n, decreasing = TRUE),] 

library(writexl)
write_xlsx(pertTup, "/home/petear/MEGA/TormodGroup/InputData/MoACounted.xlsx")
##############################################################################################
library(reshape2)


# 1. Generate the count dataframe
count_df <- as.data.frame(table(pertDF$"...6"))
colnames(count_df) <- c("...6", "n")

# Sort the dataframe by count in descending order
count_df <- count_df[order(-count_df$n), ]

# Print the count dataframe
head(count_df, 10)

library(tidyverse)

# Reshape the dataframe
wide_df <- pertDF %>%
  group_by(...6) %>%
  mutate(row = row_number()) %>%
  spread(key = ...6, value = ...1) %>%
  select(-row)

head(wide_df)


library(openxlsx)

write.xlsx(list("Counts" = count_df, "Wide Data" = wide_df), "MoACounts.xlsx")