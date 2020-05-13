library(Matrix)

makeSymmetric <- function(X,Y, O){
  X <- (X + t(X))/2
  Y <- (Y + t(Y))/2
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(X,Y)
  DR <- scTenifoldNet::dRegulation(MA)
  write.csv(DR,file = O)
}

X <- readMM('tensorOutput/X_10X500ANEURON_Itensor.mtx')
rownames(X) <- colnames(X) <- readLines('tensorOutput/genes_10X500ANEURON_Itensor.mtx')
Y <- readMM('tensorOutput/Y_10X500ANEURON_Itensor.mtx')
rownames(Y) <- colnames(Y) <- readLines('tensorOutput/genes_10X500ANEURON_Itensor.mtx')
O <- 'sym10X500ANEURON_Itensor_Dalignment.csv'
makeSymmetric(X,Y,O)

X <- readMM('tensorOutput/X_10X500DF_Itensor.mtx')
rownames(X) <- colnames(X) <- readLines('tensorOutput/genes_10X500DF_Itensor.mtx')
Y <- readMM('tensorOutput/Y_10X500DF_Itensor.mtx')
rownames(Y) <- colnames(Y) <- readLines('tensorOutput/genes_10X500DF_Itensor.mtx')
O <- 'sym10X500DF_Itensor_Dalignment.csv'
makeSymmetric(X,Y,O)

X <- readMM('tensorOutput/X_10X500morphineNeuron_Itensor.mtx')
rownames(X) <- colnames(X) <- readLines('tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')
Y <- readMM('tensorOutput/Y_10X500morphineNeuron_Itensor.mtx')
rownames(Y) <- colnames(Y) <- readLines('tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')
O <- 'sym10X500morphineNeuron_Itensor_Dalignment.csv'
makeSymmetric(X,Y,O)
