pathways.show = 'GUCA'
netAnalysis_contribution(cellchat.d0, signaling = pathways.show, title = 'Each L-R pair contribution D0')
netAnalysis_contribution(cellchat.d3, signaling = pathways.show, title = 'Each L-R pair contribution D3')

# Circle plot
par(mfrow=c(1,2), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(cellchat.d0, signaling = pathways.show, layout = "circle", edge.weight.max = 0.00213, weight.scale = T, arrow.size = 0.7)
netVisual_aggregate(cellchat.d3, signaling = pathways.show, layout = "circle", edge.weight.max = 0.00213, weight.scale = T, arrow.size = 0.7)

netAnalysis_signalingRole_network(cellchat.d0, signaling = pathways.show, width = 10, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.d3, signaling = pathways.show, width = 10, height = 2.5, font.size = 10)


netVisual_heatmap(cellchat.d0, signaling = pathways.show, measure = "weight", color.heatmap = "Blues")
netVisual_heatmap(cellchat.d3, signaling = pathways.show, measure = "weight", color.heatmap = "Blues")

netVisual_heatmap(cellchat.d3, measure = "weight", color.heatmap = "Blues")


i=1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathways.show, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathways.show, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
library(ComplexHeatmap)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], signaling = pathways.show) +  scale_y_continuous() + scale_x_continuous()
}
patchwork::wrap_plots(plots = gg)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], signaling = pathways.show, weight.MinMax = c(0,100), title = names(object.list)[i]) +  scale_y_continuous(limits = c(0,0.003)) + scale_x_continuous(limits = c(0,0.005))
}
patchwork::wrap_plots(plots = gg)


pathways0 <- cellchat.d0@netP$pathways
pathways3 <- cellchat.d3@netP$pathways

unique_to_cellchat0 <- setdiff(pathways0, pathways3)
unique_to_cellchat3 <- setdiff(pathways3, pathways0)
common_pathways <- intersect(pathways0, pathways3)

for (i in 1:length(common_pathways)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = common_pathways[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg1 <- netAnalysis_contribution(cellchat.d0, signaling = common_pathways[i])
  gg2 <- netAnalysis_contribution(cellchat.d3, signaling = common_pathways[i])
  ggsave(filename=paste0(common_pathways[i], "D0_L-R_contribution.pdf"), plot=gg1, width = 6, height = 4, units = 'in', dpi = 450)
  ggsave(filename=paste0(common_pathways[i], "D3_L-R_contribution.pdf"), plot=gg2, width = 6, height = 4, units = 'in', dpi = 450)
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(common_pathways)) {
  gg1 <- netAnalysis_signalingRole_scatter(cellchat.d0, title = 'D0', signaling = common_pathways[i]) +  scale_y_continuous() + scale_x_continuous()
  gg2 <- netAnalysis_signalingRole_scatter(cellchat.d3, title = 'D3', signaling = common_pathways[i]) +  scale_y_continuous() + scale_x_continuous()
  ggsave(filename=paste0(common_pathways[i], "D0_L-R_contribution.pdf"), plot=gg1, width = 6, height = 4, units = 'in', dpi = 450)
  ggsave(filename=paste0(common_pathways[i], "D3_L-R_contribution.pdf"), plot=gg2, width = 6, height = 4, units = 'in', dpi = 450)
}

#flow = sum of communication probability for all cell groups
gg1 <- rankNet(cellchat, measure = 'weight', mode = 'comparison', do.stat=TRUE, stacked = TRUE)
gg2 <- rankNet(cellchat, measure = 'weight', mode = 'comparison', do.stat=TRUE, stacked = F)
gg1+gg2

rankNet(cellchat, measure = 'weight', mode = 'comparison', do.stat=TRUE, show.raw = TRUE)

rankNet(cellchat, measure = 'weight', mode = 'comparison', do.stat=TRUE, show.raw = TRUE, return.data = TRUE)



ranknetresults <- rankNet(cellchat, measure = 'weight', mode = 'comparison', do.stat=TRUE, show.raw = TRUE, return.data = TRUE)
#write.csv(ranknetresults$signaling.contribution, file = "Rep2.12groups.3DB.ranknet.csv", row.names = FALSE)

write.csv(CellChatDB.mouse$interaction, file = '/CellChatDB.mouse.csv')

rankNet(cellchat, measure = 'weight', mode = 'comparison', sources.use = c('Crypt', 'Intercrypt', 'Junction', 'Subcrypt', 'Villus stroma'), targets.use = c('Crypt', 'Intercrypt', 'Junction', 'Subcrypt', 'Villus stroma'), signaling = c('GRN', 'Complement'), do.stat=TRUE, show.raw = FALSE)
rankNet(cellchat, measure = 'weight', mode = 'comparison', signaling = c('GRN', 'Complement'), do.stat=TRUE, show.raw = FALSE)



netVisual_bubble(cellchat, signaling = c('TGFb'), comparison = c(1, 2), angle.x = 45)

netVisual_bubble(cellchat, signaling = c('TGFb'), sources.use = c(1), comparison = c(1, 2), angle.x = 45)




seurat_obj <- LoadSeuratRds('/rep 2/Jej_1B5_abc15_keep/8um/load10xspatial/8umcreatedobjnorm.rds')

metadata <- read.csv('/rep 2/Jej_1B5_abc15_keep/8um/25.4.8 villus epi and stroma D0 and D3.csv', row.names = 1)
metadata <- metadata[rownames(metadata) %in% colnames(seurat_obj), , drop = FALSE]
seurat_obj <- AddMetaData(seurat_obj, metadata)
Idents(seurat_obj) <- seurat_obj$Stromal.regions
regions_to_keep <- c('D0', 'D3')
seurat_subset <- subset(seurat_obj, subset = FullSamples.NoC2 %in% regions_to_keep)
regions_to_keep <- c("Day 0 crypt", "Day 0 tip C4", "Day 0 Mature epi C1", "Day 0 Young villus C8", "Day 0 junction C7", "D0 Villus stroma", "D0 Intercrypt", "D0 Subcrypt", "Day 3 Crypt C5", "Day 3 intercrypt", "Day 3 Junction C10", "Day 3 subcrypt", "Day 3 tip", "day 3 villus", "Day 3 Villus epi C13")
seurat_subset <- subset(seurat_subset, subset = Stromal.regions %in% regions_to_keep)
table(seurat_subset$Stromal.regions)

seurat_subset <- subset(seurat_subset, subset = nCount_Spatial > 40)

fts <- c('Tgfb1','Acvr1b', 'Tgfbr2','Tgfbr1', 'Acvr1','Acvr1c', 'Tgfb2','Tgfb3')

fts <- c('Col4a1', 'Itga3', 'Itgb1', 'Sdc1', 'Col1a1', 'Col1a2', 'Sdc4', 
         'Itga1', 'Itga2', 'Itga9', 'Cd44', 'Col6a1', 'Itga11', 'Col6a3', 
         'Col6a2', 'Col4a2', 'Col6a4', 'Col4a5', 'Col4a6', 'Col6a5', 
         'Itgav', 'Itgb8', 'Col4a3', 'Col9a3', 'Col4a4', 'Col2a1')

fts <- c("Aldh1a1", "Aldh1a2", "Aldh1a3", "Rara", "Rxra", "Crabp2", "Rxrb", "Rarb", "Rarg", "Rorb")

fts <- c("PGE2", "Ptges3", "Ptger4", "Ptges2", "Ptger2", "Ptger1", "Ptges", "Ptger3", 
           "PGF2a", "Prxl2b", "Ptgfr", "PGI2", "Ptgis", "Ptgir", "TXA2", "Tbxas1", "Tbxa2r")

fts <- c('Il6', 'Il6ra', 'Il6st')

main <- c('Tgfb1', 'Mif', 'Nrg1', 'C3', 'Areg', 'Vegfa', 'Angptl2', 'Col4a1', 'Tnf', 'Lgals9', 'Il6', 'Igf1', 'Aldh1a1', 'Ptges2', 'Ptges3')

fts <- c("Nrg1", "Itga6", "Itgb4", "Erbb3", "Erbb2", "Nrg2", "Itgav", "Itgb3")

fts <- c("Igf2", "Itga6", "Itgb4", "Igf1", "Igf2r", "Itgav", "Itgb3", "Igf1r")


DotPlot(seurat_subset, features = fts, group.by = 'FullSamples.NoC2', scale = F)
DotPlot(seurat_subset, features = main, group.by = 'FullSamples.NoC2', scale = F, dot.scale = 12)

VlnPlot(seurat_subset, features = main, group.by = 'FullSamples.NoC2', stack = F)



DoHeatmap(seurat_subset, features=main, group.by='FullSamples.NoC2', slot = 'data')
