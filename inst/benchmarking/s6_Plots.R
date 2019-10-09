PCR <- read.csv('metrics/PCR.csv')
SCC <- read.csv('metrics/SCC.csv')
MI <- read.csv('metrics/MI.csv')
GENIE3 <- read.csv('metrics/GENIE3.csv')



PCR <- PCR[PCR$q == 0,]
SCC <- SCC[SCC$q == 0,]
MI <- MI[MI$q == 0,]
GENIE3 <- GENIE3[GENIE3$q == 0,]

par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.8,0.5,0))
plot(PCR$nCells, PCR$acc, ylim=c(0.5,0.8), type='b',
     ylab = 'Accuracy', xlab = 'Number of Cells', pch = 15)
arrows(x0 = PCR$nCells, x1 = PCR$nCells, y0 = PCR$accLB, y1 = PCR$accUB,
       length = 0.03, code = 3, angle = 90, col = 'black')
points(SCC$nCells, SCC$acc, col = 'red', type='b', pch = 16)
arrows(x0 = SCC$nCells, x1 = SCC$nCells, y0 = SCC$accLB, y1 = SCC$accUB,
       length = 0.03, code = 3, angle = 90, col = 'red')
points(MI$nCells, MI$acc, col = 'blue', type='b', pch = 17)
arrows(x0 = MI$nCells, x1 = MI$nCells, y0 = MI$accLB, y1 = MI$accUB,
       length = 0.03, code = 3, angle = 90, col = 'blue')
points(GENIE3$nCells, GENIE3$acc, col = 'forestgreen', type='b', pch = 18)
arrows(x0 = GENIE3$nCells, x1 = GENIE3$nCells, y0 = GENIE3$accLB, y1 = GENIE3$accUB,
       length = 0.03, code = 3, angle = 90, col = 'forestgreen')

plot(PCR$nCells, PCR$recall, ylim=c(0,1), type='b',
     ylab = 'Accuracy', xlab = 'Number of Cells', pch = 15)
arrows(x0 = PCR$nCells, x1 = PCR$nCells, y0 = PCR$recallLB, y1 = PCR$recallUB,
       length = 0.03, code = 3, angle = 90, col = 'black')
points(SCC$nCells, SCC$recall, col = 'red', type='b', pch = 16)
arrows(x0 = SCC$nCells, x1 = SCC$nCells, y0 = SCC$recallLB, y1 = SCC$recallUB,
       length = 0.03, code = 3, angle = 90, col = 'red')
points(MI$nCells, MI$recall, col = 'blue', type='b', pch = 17)
arrows(x0 = MI$nCells, x1 = MI$nCells, y0 = MI$recallLB, y1 = MI$recallUB,
       length = 0.03, code = 3, angle = 90, col = 'blue')
points(GENIE3$nCells, GENIE3$recall, col = 'forestgreen', type='b', pch = 18)
arrows(x0 = GENIE3$nCells, x1 = GENIE3$nCells, y0 = GENIE3$recallLB, y1 = GENIE3$accUB,
       length = 0.03, code = 3, angle = 90, col = 'forestgreen')
