library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)


library(CellChat)
library(patchwork)

data.dir <- '/Users/tylertran/Desktop/spatial r/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/comparison'
dir.create(data.dir)
setwd(data.dir)

#6 groups (intercrypt, subcrypt combined), 3 databases
cellchat.d0 <- readRDS('6groupsDB3samename8umRegionsD0cellchat_trim001.200.20.rds')
cellchat.d3 <- readRDS('6groupsDB3samename8umRegionsD3cellchat_trim001.200.20.rds')
cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order = c("Crypt stroma", 'Crypt', 'Junction', 'Villus stroma', "Villus epithelium", "Tip"))
cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order = c("Crypt stroma", 'Crypt', 'Junction', 'Villus stroma', "Villus epithelium", "Tip"))

#6 groups secreted only db
cellchat.d0 <- readRDS('6groupsDBsec.samename8umRegionsD0cellchat_trim001.200.20.rds')
cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order = c("Crypt stroma", 'Crypt', 'Junction', 'Villus stroma', "Villus epithelium", "Tip"))
cellchat.d3 <- readRDS('6groupsDBsec.samename8umRegionsD3cellchat_trim001.200.20.rds')
cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order = c("Crypt stroma", 'Crypt', 'Junction', 'Villus stroma', "Villus epithelium", "Tip"))

#6 groups non protein DB
cellchat.d0 <- readRDS('6groupsDBnonproteinsamename8umRegionsD0cellchat_trim001.200.20.rds')
cellchat.d3 <- readRDS('6groupsDBnonproteinsamename8umRegionsD3cellchat_trim001.200.20.rds')
cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order = c("Crypt stroma", 'Crypt', 'Junction', 'Villus stroma', "Villus epithelium", "Tip"))
cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order = c("Crypt stroma", 'Crypt', 'Junction', 'Villus stroma', "Villus epithelium", "Tip"))



cellchat.d0 <- filterCommunication(cellchat.d0, min.cells = 10) #minimum number of cells required in each cell group for cell-cell communication
cellchat.d3 <- filterCommunication(cellchat.d3, min.cells = 10)
cellchat.d0 <- computeCommunProbPathway(cellchat.d0)
cellchat.d0 <- aggregateNet(cellchat.d0)
cellchat.d0 <- netAnalysis_computeCentrality(cellchat.d0, slot.name = 'netP')
cellchat.d3 <- computeCommunProbPathway(cellchat.d3)
cellchat.d3 <- aggregateNet(cellchat.d3)
cellchat.d3 <- netAnalysis_computeCentrality(cellchat.d3, slot.name = 'netP')


#cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order = c("Crypt", "Villus epithelium", "Tip", 'Villus stroma', 'Junction', 'Intercrypt', 'Subcrypt'))
#cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order =c('Junction', 'Intercrypt', 'Subcrypt', "Crypt", "Villus epithelium", "Tip", 'Villus stroma'))
#cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order =c("Crypt", "Villus epithelium", 'Junction', 'Intercrypt', 'Subcrypt', "Tip", 'Villus stroma'))
#cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order =c("Villus epithelium", "Crypt", 'Junction', 'Intercrypt', 'Subcrypt', "Tip", 'Villus stroma'))
#cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order =c("Villus stroma", "Crypt", 'Junction', 'Intercrypt', 'Subcrypt', "Tip", 'Villus epithelium'))
#cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order = c("Crypt", "Villus epithelium", "Tip", 'Villus stroma', 'Junction', 'Intercrypt', 'Subcrypt'))
#cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order = c("Villus epithelium", "Tip", 'Villus stroma', 'Junction', 'Intercrypt', 'Subcrypt', "Crypt"))
#cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order = c("Intercrypt", "Tip", 'Villus stroma', 'Villus epithelium', 'Junction', 'Subcrypt', "Crypt"))
#cellchat.d0 <- computeCommunProbPathway(cellchat.d0)
#cellchat.d0 <- aggregateNet(cellchat.d0)
#cellchat.d0 <- netAnalysis_computeCentrality(cellchat.d0, slot.name = 'netP')

#cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order = c("Crypt", "Villus epithelium", "Tip", 'Villus stroma', 'Junction', 'Intercrypt', 'Subcrypt'))
#cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order = c('Junction', 'Intercrypt', 'Subcrypt', "Crypt", "Villus epithelium", "Tip", 'Villus stroma'))
#cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order = c("Tip", 'Villus stroma', 'Junction', 'Intercrypt', 'Subcrypt', "Crypt", "Villus epithelium"))
#cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order =c("Crypt", "Villus epithelium", 'Junction', 'Intercrypt', 'Subcrypt', "Tip", 'Villus stroma'))
#cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order =c("Villus stroma", "Villus epithelium", 'Junction', 'Intercrypt', 'Subcrypt', "Tip", 'Crypt'))
#cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order =c("Villus epithelium", "Crypt", 'Junction', 'Intercrypt', 'Subcrypt', "Tip", 'Villus stroma'))
#cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order = c("Crypt", "Tip", 'Villus stroma', 'Junction', 'Intercrypt', 'Subcrypt', "Villus epithelium"))
#cellchat.d3 <- computeCommunProbPathway(cellchat.d3)
#cellchat.d3 <- aggregateNet(cellchat.d3)
#cellchat.d3 <- netAnalysis_computeCentrality(cellchat.d3, slot.name = 'netP')

object.list <- list(D0 = cellchat.d0, D3 = cellchat.d3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#cellchat.d0 <- netAnalysis_computeCentrality(cellchat.d0, slot.name = 'netP')
#cellchat.d3 <- netAnalysis_computeCentrality(cellchat.d3, slot.name = 'netP')
#object.list <- list(D0 = cellchat.d0, D3 = cellchat.d3)
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#cellchat

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) +  scale_y_continuous(limits = c(0,0.15)) + scale_x_continuous(limits = c(0,0.18))
}
patchwork::wrap_plots(plots = gg)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) +  scale_y_continuous(limits = c(0,1.6)) + scale_x_continuous(limits = c(0,1.9))
}
patchwork::wrap_plots(plots = gg)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], do.label = T, title = names(object.list)[i], weight.MinMax = weight.MinMax) +  scale_y_continuous(limits = c(0,0.08)) + scale_x_continuous(limits = c(0,0.07))
}
patchwork::wrap_plots(plots = gg)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) +  scale_y_continuous(limits = c(0,0.6)) + scale_x_continuous(limits = c(0,0.7))
}
patchwork::wrap_plots(plots = gg)



num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, do.label = F) +  scale_y_continuous(limits = c(0,0.6)) + scale_x_continuous(limits = c(0,0.7))
}
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Intercrypt")
#> Visualizing differential outgoing and incoming signaling changes from D0 to D3
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Subcrypt")
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tip", dot.size = 7, color.use = c('grey10', 'grey10', 'grey10'), point.shape = c(21, 21, 21, 21), do.label=F)
#> Visualizing differential outgoing and incoming signaling changes from D0 to D3
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Villus epithelium", dot.size = 7, color.use = c('grey10', 'grey10', 'grey10'), point.shape = c(21, 21, 21, 21), do.label=T)
patchwork::wrap_plots(plots = list(gg1,gg2))


gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Crypt", dot.size = 7, color.use = c('grey10', 'grey10', 'grey10'), point.shape = c(21, 21, 21, 21), do.label=T)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Crypt stroma", dot.size = 7, color.use = c('grey10', 'grey10', 'grey10'), point.shape = c(21, 21, 21, 21), do.label=T))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Villus stroma", dot.size = 7, color.use = c('grey10', 'grey10', 'grey10'), point.shape = c(21, 21, 21, 21), do.label=T)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Junction", dot.size = 7, color.use = c('grey10', 'grey10', 'grey10'), point.shape = c(21, 21, 21, 21), do.label=T)
patchwork::wrap_plots(plots = list(gg1,gg2))

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tip")

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inner EN")
#> Visualizing differential outgoing and incoming signaling changes from D0 to D3
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Outer EN")
patchwork::wrap_plots(plots = list(gg1,gg2))

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Muscle")


netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Crypt stroma", signaling.exclude = c('COLLAGEN', 'LAMININ', 'APP'))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Junction", signaling.exclude = c("GUCA", 'GALECTIN'))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Crypt Stroma", signaling.exclude = c("LAMININ", 'COLLAGEN', 'APP'))

  set.seed(1234)
library(uwot)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)



rankSimilarity(cellchat, type = "functional")

rankSimilarity(cellchat, type = "structural")


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", measure = 'weight', stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = 'weight', stacked = F, do.stat = TRUE)
gg1 + gg2

rankNet(cellchat, mode = "comparison", measure = 'weight', stacked = F, do.stat = TRUE, show.raw = T)


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pathway.union <- c('TGFb', 'GALECTIN', 'GUCA', 'COMPLEMENT', 'MIF', 'IL1', 'ANNEXIN', 'SOMATOSTATIN', 'WNT', 'BMP')
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 45, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 45, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



netVisual_bubble(cellchat, sources.use = 2, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, signaling = c('COMPLEMENT'), comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in D0", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in D0", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2




# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "D0"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "D0",ligand.logFC = 0.1, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "D3",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 2, targets.use = c(3:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 2, targets.use = c(3:7), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 2, targets.use = c(3:7), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 2, targets.use = c(3:7), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#install.packages("wordcloud")
computeEnrichmentScore(net.up, species = 'mouse')
computeEnrichmentScore(net.down, species = 'mouse')

computeEnrichmentScore(net.up, measure = 'signaling', species = 'mouse')
computeEnrichmentScore(net.down, measure = 'signaling', species = 'mouse')



cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("D0", "D3")) # set factor level
plotGeneExpression(cellchat, signaling = "COMPLEMENT", split.by = "datasets", colors.ggplot = T)

plotGeneExpression(cellchat.d0, signaling = "COMPLEMENT", colors.ggplot = T)
plotGeneExpression(cellchat.d3, signaling = "COMPLEMENT", split.by = "datasets", colors.ggplot = T)

pathways.show = 'VEGF'
netAnalysis_contribution(cellchat.d0, signaling = pathways.show, title = 'Each L-R pair contribution D0')
netAnalysis_contribution(cellchat.d3, signaling = pathways.show, title = 'Each L-R pair contribution D3')

# Circle plot
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(cellchat.d0, signaling = pathways.show, layout = "circle", title = paste(pathways.show, "D0"))
netVisual_aggregate(cellchat.d3, signaling = pathways.show, layout = "circle")

netAnalysis_signalingRole_network(cellchat.d0, signaling = pathways.show, width = 10, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.d3, signaling = pathways.show, width = 10, height = 2.5, font.size = 10)

netVisual_heatmap(cellchat, signaling = pathways.show, measure = "weight")
netVisual_heatmap(cellchat.d0, signaling = pathways.show, measure = "weight", color.heatmap = "Blues")
netVisual_heatmap(cellchat.d3, signaling = pathways.show, measure = "weight", color.heatmap = "Blues")


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], signaling = pathways.show) +  scale_y_continuous() + scale_x_continuous()
}
patchwork::wrap_plots(plots = gg)



pathways1 <- cellchat.d0@netP$pathways
pathways2 <- cellchat.d3@netP$pathways
all_pathways <- unique(c(pathways1, pathways2))
# Exclude 'COLLAGEN'
filtered_pathways <- setdiff(all_pathways, "COLLAGEN")

# Check result
print(filtered_pathways)
