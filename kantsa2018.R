## Load Kantsa et al 2018 pollinator network and trait data
## John Powers
## April 27 2018

#libraries
library(heatmaply)
library(bipartite)
data(mosquin1967)
library(RColorBrewer)
pal <- brewer.pal(3,"Set2")

#load network
net <- read.delim("scentnetwork/Kantsa_etal_Supp_Data1.csv", skip = 2)
rownames(net) <- net[,1]
net[,1] <- NULL
net <- as.matrix(t(net))
snet <- sortweb(net) #sort by row & column totals

#heatmap
heatmaply(net, scale="row")

#interaction matrix plot
visweb(net, type="nested", circles=T, boxes=F, circle.max=1.2)

#interaction web plot
plotweb(net, text.rot=90, method="cca", col.high=pal[2], col.low=pal[1], col.interaction=pal[3])

#network topography
networklevel(net)

