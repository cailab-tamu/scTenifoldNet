nGenes <- 10
nConditions <- 2
nNetworks <- 5
X <- array(dim = c(nGenes,nGenes,nConditions,nNetworks))

# Generating input
for (i in 1:nConditions){
  for(j in 1:nNetworks){
    X[,,i,j] <- matrix(runif(nGenes^2), ncol = nGenes, nrow = nGenes) 
  }
}
X

library(rTensor)
# Making tensor object
X <- as.tensor(X)


# Running CP decomposition
tensorOutput <- cp(tnsr = X, num_components = 3, max_iter = 100)

# Decomposed matrices
U1 <- tensorOutput$U[[1]]
U2 <- tensorOutput$U[[2]]
U3 <- tensorOutput$U[[3]]
if(nConditions == 2){
  U4 <- tensorOutput$U[[4]]
}

# Reconstructed Matrices
tensorOutput$est@data
