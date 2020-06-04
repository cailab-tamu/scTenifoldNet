# setwd("E:\\GitHub\\scTenifoldNet\\MATLAB\\test\\R_Matlab_Julia\\network_construction")
X<-read.delim('X.txt',header=0)
X<-as.matrix(X)
rownames(X) <- paste0('g', 1:100)

library(scTenifoldNet)
out <- pcNet (X, nComp = 3, scaleScores = TRUE, symmetric = FALSE, q = 0, verbose = TRUE)
 
#  t(c(out[1,1],out[2,1],out[3,1]))
#  t(c(out[1,2],out[2,2],out[3,2]))

out[1:5,1:5]

#5 x 5 sparse Matrix of class "dgCMatrix"
#g1           g2           g3           g4         g5
#g1  .         -0.278870471 -0.125725713  0.487261818 0.21654237
#g2 -0.2476546  .           -0.006820856 -0.007739348 0.00330269
#g3 -0.1085018  0.031873180  .            0.181248760 0.19074912
#g4  0.3253655  0.002249899  0.215115726  .           0.21577932
#g5 -0.0237925  0.037821700  0.140198583  0.246401901 .         
