#import xlsx files
library(readxl)
endocyt <- read_excel("/home/petear/MEGA/TormodGroup/InputData/EndocyticAlgoSi09Aug.xlsx", col_names = FALSE)
endocyt <- as.data.frame(endocyt)

endolys <- read_excel("/home/petear/MEGA/TormodGroup/InputData/endolysAlgoSi09Aug.xlsx", col_names = FALSE)
endolys <- as.data.frame(endolys)

all <- read_excel("/home/petear/MEGA/TormodGroup/InputData/AllAlgoSi09Aug.xlsx", col_names = FALSE)
all <- as.data.frame(all)

library(ggvenn)
library(ggplot2)

list_of_data <- list(Endocyt = endocyt[[1]], Endolys = endolys[[1]], All = all[[1]])

pdf("VennDiagram.pdf")
venn_plot <- ggvenn(list_of_data)
print(venn_plot)
dev.off()

#First 50 rows of each data frame
endocyt50 <- head(endocyt, 50)
endolys50 <- head(endolys, 50)
all50 <- head(all, 50)

pdf("VennDiagram50.pdf")
venn_plot50 <- ggvenn(list(Endocyt = endocyt50[[1]], Endolys = endolys50[[1]], All = all50[[1]]))
print(venn_plot50)
dev.off()

#First 100 rows of each data frame
endocyt100 <- head(endocyt, 100)
endolys100 <- head(endolys, 100)
all100 <- head(all, 100)

pdf("VennDiagram100.pdf")
venn_plot100 <- ggvenn(list(Endocyt = endocyt100[[1]], Endolys = endolys100[[1]], All = all100[[1]]))
print(venn_plot100)
dev.off()

#find the strings that are in both endocyt and endolys but not in all
EcEl_notA100 <- setdiff(intersect(endocyt100[[1]], endolys100[[1]]), all100[[1]])
EcEl_notA50 <- setdiff(intersect(endocyt50[[1]], endolys50[[1]]), all50[[1]])
EcA_notEl100 <- setdiff(intersect(endocyt100[[1]], all100[[1]]), endolys100[[1]])
EcA_notEl50 <- setdiff(intersect(endocyt50[[1]], all50[[1]]), endolys50[[1]])
ElA_notEc100 <- setdiff(intersect(endolys100[[1]], all100[[1]]), endocyt100[[1]])
ElA_notEc50 <- setdiff(intersect(endolys50[[1]], all50[[1]]), endocyt50[[1]])


#find strings that are in all three
EcElA100 <- intersect(intersect(endocyt100[[1]], endolys100[[1]]), all100[[1]])
EcElA50 <- intersect(intersect(endocyt50[[1]], endolys50[[1]]), all50[[1]])

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "EcEl_notA100")
writeData(wb, "EcEl_notA100", EcEl_notA100)

addWorksheet(wb, "EcEl_notA50")
writeData(wb, "EcEl_notA50", EcEl_notA50)

addWorksheet(wb, "EcA_notEl100")
writeData(wb, "EcA_notEl100", EcA_notEl100)

addWorksheet(wb, "EcA_notEl50")
writeData(wb, "EcA_notEl50", EcA_notEl50)

addWorksheet(wb, "ElA_notEc100")
writeData(wb, "ElA_notEc100", ElA_notEc100)

addWorksheet(wb, "ElA_notEc50")
writeData(wb, "ElA_notEc50", ElA_notEc50)

addWorksheet(wb, "EcElA100")
writeData(wb, "EcElA100", EcElA100)

addWorksheet(wb, "EcElA50")
writeData(wb, "EcElA50", EcElA50)

saveWorkbook(wb, "VennDiagramPerts.xlsx", overwrite = TRUE)









#find the strings that are in both endocyt and endolys
endocyt_endolys100 <- intersect(endocyt100[[1]], endolys100[[1]])


#make excel files of the strings that are in both endocyt and endolys
library(xlsx)

write.xlsx(endocyt_endolys50, "endocyt_endolys50.xlsx")
write.xlsx(endocyt_endolys100, "endocyt_endolys100.xlsx")
