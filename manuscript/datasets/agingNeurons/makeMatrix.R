library(Matrix)
fileList <- list.files(pattern = 'GS[[:print:]]+.txt.gz')
metaData <- read.csv('GSE129788_Supplementary_meta_data_Cell_Types_Etc.txt', sep = '\t')
metadata <- metaData[metaData$cluster %in% 'mNEUR',]
metadata <- metadata[metadata$animal_type %in% c('young','old'),]
metadata$NAME <- gsub('Aging_mouse_brain_portal_data_','X',metadata$NAME)


mNeur <- lapply(fileList, function(X){
  input <- read.csv(X,sep = '\t')
  input <- input[,colnames(input) %in% metadata$NAME]
  input <- as.matrix(input)
  input <- as(input,'dgCMatrix')
  input
})

allNeuron <- mNeur[[1]]
for(i in seq_along(mNeur)[-1]){
  allNeuron <- cbind(allNeuron, mNeur[[i]])
}

yCells <- metadata$NAME[metadata$animal_type == 'young']
oCells <- metadata$NAME[metadata$animal_type == 'old']

yCells <- allNeuron[,yCells]
oCells <- allNeuron[,oCells]

writeMM(yCells,'yMatrix.mtx')
writeLines(rownames(yCells),'yGenes.txt')
writeLines(colnames(yCells),'yBarcodes.txt')
