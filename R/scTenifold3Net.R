#' @export scTenifold3Net
#' @title scTenifoldNet
#' @importFrom methods as
#' @description Construct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data sets collected from different conditions based on principal component regression, tensor decomposition, and manifold alignment.
#' @param X Raw counts matrix with cells as columns and genes (symbols) as rows.
#' @param Y Raw counts matrix with cells as columns and genes (symbols) as rows.
#' @param qc_minLibSize An integer value. Defines the minimum library size required for a cell to be included in the analysis.
#' @param qc_removeOutlierCells A boolean value (TRUE/FALSE), if TRUE, the identified cells with library size greater than 1.58 IQR/sqrt(n) computed from the sample, are removed. For further details see: \code{?boxplot.stats}
#' @param qc_minPCT A decimal value between 0 and 1. Defines the minimum fraction of cells where the gene needs to be expressed to be included in the analysis.
#' @param qc_maxMTratio A decimal value between 0 and 1. Defines the maximum ratio of mitochondrial reads (mithocondrial reads / library size) present in a cell to be included in the analysis. It's computed using the symbol genes starting with 'MT-' non-case sensitive.
#' @param nc_nNet An integer value. The number of networks based on principal components regression to generate.
#' @param nc_nCells An integer value. The number of cells to subsample each time to generate a network.
#' @param nc_nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param nc_symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param nc_scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param nc_q A decimal value between 0 and 1. Defines the cut-off threshold of top q\% relationships to be returned.
#' @param td_K An integer value. Defines the number of rank-one tensors used to approximate the data using CANDECOMP/PARAFAC (CP) Tensor Decomposition. 
#' @param td_maxIter An integer value. Defines the maximum number of iterations if error stay above \code{td_maxError}.
#' @param td_maxError A decimal value between 0 and 1. Defines the relative Frobenius norm error tolerance.
#' @param ma_nDim An integer value. Defines the number of dimensions of the low-dimensional feature space to be returned from the non-linear manifold alignment.
#' @return A list with 3 slots as follows: 
#' \itemize{
#' \item{tensorNetworks:} The generated weight-averaged denoised gene regulatory networks using CANDECOMP/PARAFAC (CP) Tensor Decomposition.
#' \item{manifoldAlignment:} The generated low-dimensional features result of the non-linear manifold alignment.
#' \item{diffRegulation} The results of the differential regulation analysis.
#' }
#' @examples 
#' library(scTenifoldNet)
#' 
#' # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
#' nCells = 2000
#' nGenes = 100
#' set.seed(1)
#' X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
#' X <- round(X)
#' X <- matrix(X, ncol = nCells)
#' rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
#' 
#' # Generating a perturbed network modifying the expression of genes 10, 2 and 3
#' Y <- X
#' Y[10,] <- Y[50,]
#' Y[2,] <- Y[11,]
#' Y[3,] <- Y[5,]
#' 
#' \dontrun{
#' # scTenifoldNet
#' Output <- scTenifoldNet(X = X, Y = Y,
#'                        nc_nNet = 10, nc_nCells = 500,
#'                        td_K = 3, qc_minLibSize = 30)
#' 
#' # Structure of the output
#' str(Output)
#' 
#' # Accessing the computed weight-averaged denoised gene regulatory networks
#' 
#' # Network for sample X
#' igraph::graph_from_adjacency_matrix(adjmatrix = Output$tensorNetworks$X, weighted = TRUE)
#' # IGRAPH 15cbeea DNW- 100 2836 -- 
#' # + attr: name (v/c), weight (e/n)
#' # + edges from 15cbeea (vertex names):
#' #   [1] ng6 ->ng1 ng12->ng1 ng14->ng1 ng24->ng1 ng28->ng1
#' # [6] ng31->ng1 ng42->ng1 ng44->ng1 ng49->ng1 ng55->ng1
#' # [11] ng56->ng1 ng59->ng1 ng62->ng1 ng63->ng1 ng72->ng1
#' # [16] ng73->ng1 ng74->ng1 ng77->ng1 ng80->ng1 ng82->ng1
#' # [21] ng83->ng1 ng87->ng1 ng89->ng1 mt-1->ng1 mt-5->ng1
#' # [26] mt-7->ng1 ng27->ng3 ng28->ng3 ng31->ng3 ng32->ng3
#' # [31] ng44->ng3 ng59->ng3 ng62->ng3 ng72->ng3 ng73->ng3
#' # [36] ng74->ng3 ng77->ng3 ng82->ng3 ng87->ng3 ng89->ng3
#' # + ... omitted several edges
#' 
#' # Network for sample Y
#' igraph::graph_from_adjacency_matrix(adjmatrix = Output$tensorNetworks$Y, weighted = TRUE)
#' #IGRAPH 3ad1533 DNW- 100 725 -- 
#' # + attr: name (v/c), weight (e/n)
#' # + edges from 3ad1533 (vertex names):
#' #   [1] ng2 ->ng2 ng3 ->ng2 ng5 ->ng2 ng6 ->ng2
#' # [5] ng7 ->ng2 ng8 ->ng2 ng9 ->ng2 ng10->ng2
#' # [9] ng11->ng2 ng12->ng2 ng13->ng2 ng15->ng2
#' # [13] ng16->ng2 ng17->ng2 ng18->ng2 ng20->ng2
#' # [17] ng21->ng2 ng22->ng2 ng23->ng2 ng24->ng2
#' # [21] ng25->ng2 ng26->ng2 ng28->ng2 ng29->ng2
#' # [25] ng30->ng2 ng31->ng2 ng33->ng2 ng34->ng2
#' # [29] ng35->ng2 ng36->ng2 ng38->ng2 ng39->ng2
#' # + ... omitted several edges
#' 
#' # Accessing the manifold alignment result
#' 
#' head(Output$manifoldAlignment)
#' #            NLMA 1      NLMA 2      NLMA 3       NLMA 4        NLMA 5
#' # X_ng1  0.0068499391  0.01096706  0.03077900  0.002655469 -0.0136455614
#' # X_ng2  0.3356288575 -0.03551752 -0.18463680 -0.193353751  0.3398606363
#' # X_ng3 -0.1285177133 -0.20064344  0.20926567  0.059542294 -0.0099528441
#' # X_ng4  0.0029881645 -0.01267593  0.01195683  0.007331123  0.0003031888
#' # X_ng5 -0.1192632208 -0.18475439  0.27616148  0.112944009 -0.0281827702
#' # X_ng6  0.0005911568  0.02557475  0.07527792 -0.191180647 -0.1165095115
#' #            NLMA 6      NLMA 7       NLMA 8       NLMA 9     NLMA 10
#' # X_ng1 -0.029852128 0.007539925  0.009299591 -0.009813157 -0.01360414
#' # X_ng2 -0.313361443 0.146429589  0.006286777  0.162023788 -0.04307899
#' # X_ng3 -0.008733285 0.172084611  0.508056218  0.199322512 -0.07935797
#' # X_ng4 -0.004680652 0.005344541  0.002634755 -0.003376544 -0.01100757
#' # X_ng5 -0.126328797 0.190769152 -0.468107666  0.170278281 -0.06744795
#' # X_ng6 -0.051266264 0.063822269  0.011060924 -0.134880459 -0.02579998
#' #           NLMA 11      NLMA 12      NLMA 13      NLMA 14      NLMA 15
#' # X_ng1 -0.0199528840  0.008035130  0.004631187  0.000807797  0.011960838
#' # X_ng2 -0.0138200390 -0.002847701 -0.004404942  0.008024704  0.006040799
#' # X_ng3  0.0232384468 -0.031398116 -0.007026934  0.028956700 -0.002112626
#' # X_ng4  0.0012864539 -0.018915289  0.003835404  0.004054159 -0.002546324
#' # X_ng5  0.0232899093 -0.040974531 -0.006759459  0.025415953 -0.007518957
#' # X_ng6 -0.0001650355  0.023277338  0.006646904 -0.002683418 -0.112688129
#' #          NLMA 16      NLMA 17     NLMA 18      NLMA 19       NLMA 20
#' # X_ng1 -0.016962988 -0.016649748  0.01140020 -0.006632691 -0.0005015655
#' # X_ng2  0.007543775 -0.016188689  0.02517684  0.014814415  0.0162617154
#' # X_ng3 -0.005598267 -0.006975026  0.05218029  0.006731063  0.0183436415
#' # X_ng4  0.003207934 -0.001784120  0.01093237 -0.001192860  0.0028746990
#' # X_ng5 -0.009555879 -0.007429166  0.05206441  0.006534604  0.0170071357
#' # X_ng6 -0.065437425  0.110728870 -0.12746932  0.335610531  0.1341842827
#' #          NLMA 21      NLMA 22      NLMA 23      NLMA 24      NLMA 25
#' # X_ng1  0.003113385 -0.023311350 -0.026415944 7.085995e-04  0.053898102
#' # X_ng2  0.001390569  0.001191301 -0.015621435 2.359703e-03 -0.013418093
#' # X_ng3 -0.007483171  0.011496519  0.004164546 2.764407e-02 -0.004527981
#' # X_ng4  0.020316634 -0.002796092  0.032119363 4.203867e-05 -0.002251366
#' # X_ng5 -0.004963436  0.016525449  0.009683698 2.564700e-02  0.002286340
#' # X_ng6  0.229199525  0.340639745 -0.041216345 3.599596e-03  0.008572652
#' #          NLMA 26       NLMA 27      NLMA 28     NLMA 29      NLMA 30
#' # X_ng1  0.065832029 -0.0080248854  0.107300843 -0.02902323 -0.005337500
#' # X_ng2 -0.007982259 -0.0026295392 -0.001765851  0.01491257 -0.003546343
#' # X_ng3  0.009770602  0.0008819272  0.014564070 -0.01568192 -0.017450667
#' # X_ng4  0.015802609  0.0012975576 -0.003406675 -0.01774975 -0.003300053
#' # X_ng5  0.003401007  0.0001761177  0.013622016 -0.01224127 -0.013909178
#' # X_ng6 -0.089450710 -0.0763838722 -0.107751916 -0.05841353 -0.059217012
#' 
#' # Differential Regulation
#' head(Output$diffRegulation,n = 10)
#' # gene    distance        Z        FC      p.value      p.adj
#' # 2   ng2 0.023526702 2.762449 12.193413 0.0004795855 0.02414332
#' # 50 ng50 0.023514429 2.761550 12.180695 0.0004828665 0.02414332
#' # 11 ng11 0.022443941 2.681598 11.096894 0.0008647241 0.02882414
#' # 3   ng3 0.020263415 2.508478  9.045415 0.0026335445 0.06583861
#' # 10 ng10 0.019194561 2.417929  8.116328 0.0043868326 0.07711821
#' # 5   ng5 0.019079975 2.407977  8.019712 0.0046270923 0.07711821
#' # 31 ng31 0.013632541 1.865506  4.094085 0.0430335257 0.61476465
#' # 96 mt-6 0.011401177 1.589757  2.863536 0.0906081350 0.90977795
#' # 59 ng59 0.009835354 1.368238  2.130999 0.1443466682 0.90977795
#' # 62 ng62 0.007995812 1.067193  1.408408 0.2353209153 0.90977795
#' 
#' # Plotting
#' # Genes with FDR < 0.1 are labeled as red
#' set.seed(1)
#' qChisq <- rchisq(100,1)
#' geneColor <- rev(ifelse(Output$diffRegulation$p.adj < 0.1, 10,1))
#' qqplot(qChisq, Output$diffRegulation$FC, pch = 16, main = 'H0', col = geneColor, 
#'        xlab = expression(X^2~Quantiles), ylab = 'FC', xlim=c(0,8), ylim=c(0,13))
#' qqline(qChisq)
#' legend('bottomright', legend = c('FDR < 0.1'), pch = 16, col = 'red', bty='n', cex = 0.7)
#' }

scTenifold3Net <- function(X, Y, qc_minLibSize = 1000, qc_removeOutlierCells = TRUE,
                          qc_minPCT = 0.05, qc_maxMTratio = 0.1, nc_nNet = 10,
                          nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                          nc_q = 0.05, td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5, ma_nDim = 30){
  # Single-cell Quality Control
  X <- scQC(X, minLibSize = qc_minLibSize, removeOutlierCells = qc_removeOutlierCells, minPCT = qc_minPCT, maxMTratio = qc_maxMTratio)
  Y <- scQC(Y, minLibSize = qc_minLibSize, removeOutlierCells = qc_removeOutlierCells, minPCT = qc_minPCT, maxMTratio = qc_maxMTratio)
  
  # Counts per million (CPM) normalization
  X <- cpmNormalization(X)
  Y <- cpmNormalization(Y)

  # Comparing gene ids.
  xNames <- rownames(X)
  yNames <- rownames(Y)

  sharedGenes <- intersect(xNames, yNames)
  nGenes <- length(sharedGenes)
  
  # Filtering out non-shared genes
  X <- X[sharedGenes,]
  Y <- Y[sharedGenes,]
  
  # Construction of gene-regulatory networks based on principal component regression (pcNet) and random subsampling.
  set.seed(1)
  xList <- makeNetworks(X = X, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))
  set.seed(1)
  yList <- makeNetworks(X = Y, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))

  # CANDECOMP/PARAFRAC Tensor Decomposition
      # for(M in c('I','3d','4d')){
  set.seed(1)
  tensorOut <- tensor3Decomposition(xList, yList, K = td_K, maxIter = td_maxIter, maxError = td_maxError)
      # Matrix::writeMM(tensorOut$X,paste0('X_',id,'_',M,'tensor.mtx'))
      # Matrix::writeMM(tensorOut$Y,paste0('Y_',id,'_',M,'tensor.mtx'))
      # writeLines(sharedGenes, paste0('genes_',id,'_',M,'tensor.mtx'))
  
  # Split of tensor output
  tX <- as.matrix(tensorOut$X)
  tY <- as.matrix(tensorOut$Y)
  
  # Making it symmetric to fulfill the requirements of the MA
  tX <- (tX + t(tX))/2
  tY <- (tY + t(tY))/2
  
  # Non-linear manifold alignment
      # for(A in c('O','D','P')){
  set.seed(1)
  mA <- manifoldAlignment(tX , tY, d = ma_nDim)
  rownames(mA) <- c(paste0('X_', sharedGenes),paste0('y_', sharedGenes))
      # outFile <-paste0(id,'_',M,'tensor_',A,'alignment.csv')
      # write.csv(mA, outFile)
  
  # Differential regulation testing
  dR <- dRegulation(manifoldOutput = mA)
      # write.csv(dC, paste0('dCoex_',id,'_',M,'tensor_',A,'alignment.csv'))
      # }
      # }
  
  # Return preparation
  outputResult <- list()
  outputResult$tensorNetworks <- tensorOut
  outputResult$manifoldAlignment <- mA
  outputResult$diffRegulation <- dR

  # Return
  return(outputResult)
}
