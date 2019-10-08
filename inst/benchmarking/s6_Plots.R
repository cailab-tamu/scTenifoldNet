PCR <- read.csv('metrics/PCR.csv')
SCC <- read.csv('metrics/SCC.csv')
MI <- read.csv('metrics/MI.csv')



PCR <- PCR[PCR$q == 0,]
SCC <- SCC[SCC$q == 0,]
MI <- MI[MI$q == 0,]

plot(PCR$nCells, PCR$acc, ylim=c(0,1))
points(SCC$nCells, SCC$acc, col = 'red')
points(MI$nCells, MI$acc, col = 'blue')



plot(PCR$nCells, PCR$recall, ylim=c(0,1))
points(SCC$nCells, SCC$recall, col = 'red')


PC