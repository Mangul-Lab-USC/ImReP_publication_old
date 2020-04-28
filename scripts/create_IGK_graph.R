library("igraph")
library("plyr")
dataSet <- read.csv("../raw_data/betaTissuePairs_SI_t1t2_median.csv")
ncdr <- read.csv("../raw_data/clonoCounts.csv")
keeps <- c("tissue1", "tissue2", "IGK")
nKeep <- c("tissue", "IGK")
ncdr <- ncdr[nKeep]
names(ncdr) <- c("tissue", "clonoCount")
dataSet <- dataSet[keeps]
names(dataSet) <- c("tissue1", "tissue2", "weight")
dataSet[, "weight"] <- 1.0 - dataSet[, "weight"]
row_keep <- apply(dataSet[c(3)], 1, function(z) any(z>=0.001))
nonZeroData <- dataSet[row_keep,]
nonZeroData <- nonZeroData[order(nonZeroData$tissue1),]
ncdr <- ncdr[order(ncdr$tissue),]
gD <- simplify(graph.data.frame(nonZeroData, directed = FALSE))
cCount <- ncdr[ncdr$tissue %in% nonZeroData$tissue1, ]
cCountApp <- ncdr[ncdr$tissue %in% nonZeroData$tissue2 & !(ncdr$tissue %in% nonZeroData$tissue1), ]
cCount <- rbind(cCount, cCountApp)
cCount <- cCount[order(cCount$tissue),]
gD <- set.vertex.attribute(gD, "ncdr", index = order(V(gD)$name), value = cCount[,"clonoCount"])
summary(gD)
library("rgexf")
nodes_df <- data.frame(ID = c(1:vcount(gD)), name = V(gD)$name)
edges_df <- as.data.frame(get.edges(gD, c(1:ecount(gD))))
nodes_att <- data.frame(NCDR = V(gD)$ncdr)
edges_att <- data.frame(WGH = E(gD)$weight)
write.gexf(nodes = nodes_df, edges = edges_df, nodesAtt = nodes_att, edgesAtt = edges_att, 
           defaultedgetype = "undirected", output = "../summary_data/IGK_pairs_final.gexf")

colnames(edges_df) <- c('Tissue1', 'Tissue2')
node_info = cbind(nodes_df, nodes_att)
edge_info = cbind(edges_df, edges_att)
summ_df_tmp = merge(node_info, edge_info, by.x = 'ID', by.y = 'Tissue1')
summ_df = merge(node_info, summ_df_tmp, by.x = 'ID', by.y = 'Tissue2')
colnames(summ_df) = c('tissue1_ID', 'tissue1_name', 'tissue1_nCDR3', 'tissue2_ID', 
                      'tissue2_name', 'tissue2_nCDR3', 'edge_weight')
write.csv(summ_df, '../summary_data/FigureS10_data.csv', row.names=FALSE)
