#read in excel file with several sheets
library(readxl)
library(writexl)
library(dplyr)

MoAFile <- read_excel("/home/petear/MEGA/TormodGroup/InputData/ClueHackedMoAfromMorph.xlsx", col_names = FALSE, sheet = 1)
MoAFile <- as.data.frame(MoAFile)

vennPertsPath <- "/home/petear/MEGA/TormodGroup/InputData/VennDiagramPerts.xlsx"
sheetNames <- excel_sheets(vennPertsPath)
vennPertsSheets <- lapply(sheetNames, function(sheet) {
  read_excel(vennPertsPath, sheet = sheet)
})


# Update vennPertsSheets with matching rows and remove rows with all identical values
updatedVennPertsSheets <- lapply(vennPertsSheets, function(df) {
  matchingRows <- filter(MoAFile, MoAFile[[1]] %in% df[[1]])
  combinedDf <- bind_rows(df, matchingRows)
  filteredDf <- filter(combinedDf, !apply(combinedDf, 1, function(row) length(unique(row)) == 1))
  return(filteredDf)
})

# Determine column to count based on sheet index
getColumnToCount <- function(sheetIndex, totalSheets) {
  if (sheetIndex %in% c(totalSheets, totalSheets - 1)) {
    return(10)
  } else {
    return(9)
  }
}

# Add counts to the end of each sheet and sort them in descending order
totalSheets <- length(updatedVennPertsSheets)
updatedVennPertsSheets <- lapply(seq_along(updatedVennPertsSheets), function(index) {
  df <- updatedVennPertsSheets[[index]]
  colToCount <- getColumnToCount(index, totalSheets)
  stringCountsDf <- as.data.frame(table(df[[colToCount]]))
  colnames(stringCountsDf) <- c("String", "Count")
  
  # Sort the counts in descending order
  stringCountsDf <- stringCountsDf[order(-stringCountsDf$Count), ]
  
  separator <- data.frame(matrix(ncol = ncol(df), nrow = 1))
  colnames(separator) <- colnames(df)
  combinedDf <- bind_rows(df, separator, stringCountsDf)
  return(combinedDf)
})

# Save to Excel with original sheet names
names(updatedVennPertsSheets) <- sheetNames
write_xlsx(updatedVennPertsSheets, "PertMoAVenn.xlsx")