library(Matrix)
SERGIO <- readMM('simulationSERGIO_12x3000.mtx')

gName <- formatC(1:100, 3, format = 's')
gName <- gsub('[[:space:]]','0', gName)
gName <- paste0('g',gName)

writeLines(gName, 'geneList.txt')



sapply(1:11, function(X){
  sCells <- c(((3000 * X)+1):(3000 * (X+1)))
  sCells <- Matrix(as.matrix(SERGIO[,sCells]))
  writeMM(sCells, file = paste0('data/simulation_', X+1))
})
