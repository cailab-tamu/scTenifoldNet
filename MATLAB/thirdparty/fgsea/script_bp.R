library(data.table)
library(fgsea)
library(ggplot2)

# data(examplePathways)
# data(exampleRanks)


rnk.file<-"input.txt"
# gmt.file<-"scsig.all.v1.0.symbols.gmt"
gmt.file<-"c5.bp.v7.0.symbols.gmt"
# gmt.file<-"c5.mf.v7.0.symbols.gmt"

ranks <- read.table(rnk.file,
                    header=TRUE, colClasses = c("character", "numeric"))
exampleRanks <- setNames(ranks$t, ranks$ID)
# str(exampleRanks)

examplePathways <- gmtPathways(gmt.file)
# str(head(pathways))


fgseaRes <- fgsea(pathways = examplePathways, 
                  stats = exampleRanks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)

topPathwaysUp <- fgseaRes[ES > 0, ][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0, ][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
#              gseaParam = 0.5)

library(data.table)
fwrite(fgseaRes, file="output.txt", sep="\t", sep2=c("", " ", ""))
write.csv(topPathways,file="output_top20.txt")
