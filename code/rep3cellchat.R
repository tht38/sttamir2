library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(bsicons)

#seurat_obj <- Load10X_Spatial(data.dir = '/rep 3')
#seurat_obj <- subset(seurat_obj, subset = nCount_Spatial > 40)

#seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


metadata <- read.csv('/rep 3/Times.csv', row.names = 1)
metadata <- metadata[rownames(metadata) %in% colnames(seurat_obj), , drop = FALSE]
seurat_obj <- AddMetaData(seurat_obj, metadata)

metadata <- read.csv('/rep 3/Rep3stroma and epi regions.csv', row.names = 1)
metadata <- metadata[rownames(metadata) %in% colnames(seurat_obj), , drop = FALSE]
seurat_obj <- AddMetaData(seurat_obj, metadata)

#saveRDS(seurat_obj, file = '/rep 3/8umnorm40countmin.rds')


seurat_obj <- subset(seurat_obj, subset = Times != "")
seurat_obj <- subset(seurat_obj, subset = Graph.based != "")

seurat_obj <- LoadSeuratRds('/rep 3/8umnorm40countmin.rds')


#seurat_subset <- subset(seurat_obj, subset = FullSamples.NoC2 == "D0")
#saveRDS(seurat_subset, file = '/rep 2/Jej_1B5_abc15_keep/8um/load10xspatial/8umcreatedobjnormD0.rds')
#seurat_subset <- LoadSeuratRds('/rep 2/Jej_1B5_abc15_keep/8um/load10xspatial/8umcreatedobjnormD0.rds')

table(seurat_obj@meta.data$regions)
seurat_subset <- subset(seurat_obj, subset = regions %in% c("D0 crypt", 'D0 crypt stroma', 'D0 junction', 'D0 tip', 'D0 villus bottom', 'D0 villus epi', 'D0 villus stroma'))
table(seurat_subset@meta.data$regions)

seurat_subset <- subset(seurat_obj, subset = regions %in% c("D3 crypt", 'D3 crypt stroma', 'D3 junction', 'D3 Tip', 'D3 villus epi', 'D3 villus stroma'))
table(seurat_subset@meta.data$regions)


seurat_subset <- subset(seurat_obj, subset = regions %in% c("D0 crypt", 'D0 crypt stroma', 'D0 junction', 'D0 tip', 'D0 villus bottom', 'D0 villus epi', 'D0 villus stroma', "D3 crypt", 'D3 crypt stroma', 'D3 junction', 'D3 Tip', 'D3 villus epi', 'D3 villus stroma'))


#regions_to_keep <- c("Day 0 crypt", "Day 0 tip C4", "Day 0 Mature epi C1", "Day 0 Young villus C8", "Day 0 junction C7", "D0 Villus stroma", "D0 Intercrypt", "D0 Subcrypt")
#seurat_subset <- subset(seurat_obj, subset = Stromal.regions %in% regions_to_keep)
#table(seurat_subset$Stromal.regions)


seurat_subset$regions[seurat_subset$regions == "D0 crypt"] <- "Crypt"
seurat_subset$regions[seurat_subset$regions == "D0 crypt stroma"] <- "Crypt Stroma"
seurat_subset$regions[seurat_subset$regions == "D0 junction"] <- "Junction"
seurat_subset$regions[seurat_subset$regions == "D0 tip"] <- "Tip"
seurat_subset$regions[seurat_subset$regions == "D0 villus bottom"] <- "Villus epithelium"
seurat_subset$regions[seurat_subset$regions == "D0 villus epi"] <- "Villus epithelium"
seurat_subset$regions[seurat_subset$regions == "D0 villus stroma"] <- "Villus stroma"
table(seurat_subset$regions)
Idents(seurat_subset) <- seurat_subset$regions


seurat_subset$regions[seurat_subset$regions == "D3 crypt"] <- "Crypt"
seurat_subset$regions[seurat_subset$regions == "D3 crypt stroma"] <- "Crypt Stroma"
seurat_subset$regions[seurat_subset$regions == "D3 junction"] <- "Junction"
seurat_subset$regions[seurat_subset$regions == "D3 Tip"] <- "Tip"
seurat_subset$regions[seurat_subset$regions == "D3 villus epi"] <- "Villus epithelium"
seurat_subset$regions[seurat_subset$regions == "D3 villus stroma"] <- "Villus stroma"
table(seurat_subset$regions)
Idents(seurat_subset) <- seurat_subset$regions


seurat_subset$regions[seurat_subset$regions == "D0 crypt"] <- "D0 Crypt"
seurat_subset$regions[seurat_subset$regions == "D0 crypt stroma"] <- "D0 Crypt stroma"
seurat_subset$regions[seurat_subset$regions == "D0 junction"] <- "D0 Junction"
seurat_subset$regions[seurat_subset$regions == "D0 tip"] <- "D0 Tip"
seurat_subset$regions[seurat_subset$regions == "D0 villus bottom"] <- "D0 Villus epithelium"
seurat_subset$regions[seurat_subset$regions == "D0 villus epi"] <- "D0 Villus epithelium"
seurat_subset$regions[seurat_subset$regions == "D0 villus stroma"] <- "D0 Villus stroma"
seurat_subset$regions[seurat_subset$regions == "D3 crypt"] <- "D3 Crypt"
seurat_subset$regions[seurat_subset$regions == "D3 crypt stroma"] <- "D3 Crypt Stroma"
seurat_subset$regions[seurat_subset$regions == "D3 junction"] <- "D3 Junction"
seurat_subset$regions[seurat_subset$regions == "D3 Tip"] <- "D3 Tip"
seurat_subset$regions[seurat_subset$regions == "D3 villus epi"] <- "D3 Villus epithelium"
seurat_subset$regions[seurat_subset$regions == "D3 villus stroma"] <- "D3 Villus stroma"
table(seurat_subset$regions)
Idents(seurat_subset) <- seurat_subset$regions


library(Seurat)
library(DESeq2)
library(BiocNeighbors)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(reticulate)
library(ggalluvial)
library(NMF)
library(presto)


library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

color.use <- scPalette(nlevels(seurat_subset)); names(color.use) <- levels(seurat_subset)
#Seurat::SpatialDimPlot(seurat_obj, label = T, label.size = 3, cols = color.use)
Seurat::SpatialDimPlot(seurat_subset, label = T, label.size = 3, cols = color.use)


data.input = Seurat::GetAssayData(seurat_subset, layer = "data", assay = "Spatial") # normalized data matrix

labels = Seurat::Idents(seurat_subset)
meta = data.frame(labels = Seurat::Idents(seurat_subset), samples = "sample1", row.names = names(Seurat::Idents(seurat_subset))) # manually create a dataframe consisting of the cell labels
meta$samples <- factor(meta$samples)
unique(meta$labels) # check the cell labels
unique(meta$samples) # check the sample labels

spatial.locs = Seurat::GetTissueCoordinates(seurat_subset, scale = NULL, cols = c("imagerow", "imagecol"))
#spatial.locs <- spatial.locs[, !colnames(spatial.locs) %in% "cell"]
spatial.locs <- as.matrix(spatial.locs[, c("x", "y")])
rownames(spatial.locs) <- rownames(spatial.locs)
str(spatial.locs)

scalefactors = jsonlite::fromJSON('/rep 3/spatial/scalefactors_json.json')
spot.size = 8 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

#d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
#min(d.spatial[d.spatial!=0])
#batch_size <- 500
#num_batches <- ceiling(nrow(spatial.locs) / batch_size)
#d.spatial_list <- lapply(1:num_batches, function(i) {
#  start <- (i - 1) * batch_size + 1
#  end <- min(i * batch_size, nrow(spatial.locs))
#  computeCellDistance(coordinates = spatial.locs[start:end, ], 
#                      ratio = spatial.factors$ratio, 
#                      tol = spatial.factors$tol)
#})
#min_distance <- min(sapply(d.spatial_list, function(mat) min(mat[mat != 0], na.rm = TRUE)), na.rm = TRUE)




cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
#                           datatype = "RNA", coordinates = spatial.locs, spatial.factors = spatial.factors)
print(cellchat)






CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)


# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling) that can be only estimated from gene expression data. 
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", 'ECM-Receptor', 'Cell-Cell Contact'), key = "annotation") # use Secreted Signaling

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 7) # do parallel
future::plan("multicore", workers = 7)
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
#cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0, thresh.fc = 0, thresh.p = 0.05)
options(future.globals.maxSize = 13000 * 1024^2)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
#library(BiocParallel)
#register(SnowParam(8))  # Set number of parallel workers

# project gene expression data onto PPI network (optional)
#cellchat <- smoothData(cellchat, adj = PPI.mouse)

#queryknn error
#rownames(cellchat@meta) <- as.character(rownames(cellchat@meta))
#cellchat@meta[] <- lapply(cellchat@meta, function(x) {
#  if (is.factor(x)) as.character(x) else x
#})

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.001,
                              distance.use = TRUE, interaction.range = 200, scale.distance = 0.344,
                              contact.dependent = TRUE, contact.range = 20)


computeAveExpr(cellchat, features = c("Wnt2"), type =  "truncatedMean", trim = 0.1)
computeAveExpr(cellchat, features = c("Bmp2"), type =  "truncatedMean", trim = 0.05)
computeAveExpr(cellchat, features = c("Bmp2"), type =  "truncatedMean", trim = 0.01)
computeAveExpr(cellchat, features = c("Bmp2"), type =  "truncatedMean", trim = 0.001)
computeAveExpr(cellchat, features = c("Bmp2"), type =  "truncatedMean", trim = 0.0001)
computeAveExpr(cellchat, features = c("Bmp2"), type =  "truncatedMean", trim = 0.00001)


cellchat <- filterCommunication(cellchat, min.cells = 10) #minimum number of cells required in each cell group for cell-cell communication

length(cellchat@netP)
head(cellchat@netP$prob)
cellchat@net
cellchat@netP$pathways

#cellchat <- updateClusterLabels(cellchat, new.order = c('Isg15Epi', 'Isg15negEpi', 'Isg15Stroma', 'Isg15negStroma'))


# infers cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

cellchat@netP$pathways

pathways.show <- c("BMP") 
# Circle plot
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 4, vertex.size.max = 0.1, alpha.image = 0.1, vertex.label.cex = 3.5, point.size = 0.000008)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 10, height = 2.5, font.size = 10)


spatialFeaturePlot(cellchat, pairLR.use = "BMP2_BMPR1A_BMPR2", point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
# Take an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(cellchat, pairLR.use = "BMP2_BMPR1A_BMPR2", point.size = 1, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    xlim = c(500, 1000), ylim = c(2900, 3000))  # Adjust to your desired zoom area


# shows cell-cell communication network
mat <- cellchat@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = rowSums(cellchat@net$weight), edge.weight.max = max(mat), weight.scale = T, title.name = rownames(mat)[i], arrow.size = 0.01)
}

# displays pathways in whihc cells have significant communication

pathways <- cellchat@netP$pathways
for (i in pathways) {
  pathways.show <- c(i) 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  
}

#saveRDS(cellchat, file = '/rep 2/Jej_1B5_abc15_keep/8um/load10xspatial/8umD0cellchat_trim01.200.20.rds')

# choosing a pthway to look at
pathways.show <- c("EGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
gg <- netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
gg <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
ggsave(filename=("test_figure.pdf"), plot=gg, width=5,height=5, units='in', dpi=300)
pathways.show <- c("PROS")

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow = c(1,1), xpd=TRUE)
cairo_pdf(file="circle_egf.pdf", width=10, height=10, fallback_resolution = 300)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()



# Access all the signaling pathways showing significant communications, this saves all the graphs made so far
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,15)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 4, units = 'in', dpi = 450)
}



#
#https://theislab.github.io/interaction-tools/14-CellChat.html
group_size <- as.numeric(table(labels))


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#global patterns
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 5)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 5)

#install.packages("uwot")
library(pander)
#functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional", thresh = 0.25)
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "functional", k = 4)
str(cellchat@netP$similarity$functional$matrix)
pander(cellchat@netP$similarity$functional$matrix$single[1:5,1:5])
pander(cellchat@netP$similarity$functional$dr$single[1:5, ])
cellchat@netP$similarity$functional$group
#structural similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural", thresh = 0.25)
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "structural", k = 4)
pander(cellchat@netP$similarity$structural$matrix$single[1:5,1:5])
pander(cellchat@netP$similarity$structural$dr$single[1:5, ])
cellchat@netP$similarity$structural$group


netVisual_aggregate(cellchat, "TGFb", vertex.receiver = 1:3, layout = 'hierarchy')
netAnalysis_river(cellchat, pattern = "outgoing")


netVisual_embedding(cellchat, type = "structural", pathway.remove.show = FALSE,
                    label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


#base cellchat tutorial
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

mat <- cellchat@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, title.name = rownames(mat)[i])
}

# Chord diagram
par(mfrow=c(1,1), xpd = TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

pathways.show <- c("GUCA")
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 'Isg15Stroma', targets.use = c('Isg15Epi', 'Isg15negEpi', 'Isg15Stroma','Isg15negStroma'), remove.isolate = FALSE)
#> Comparing communications on a single object

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("BMP", "WNT"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


selectK(cellchat, pattern = "outgoing")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")




par(mfrow = c(1,2), xpd=TRUE)
rankNet(cellchat, measure  = 'count', mode = "single")
rankNet(cellchat, measure  = 'weight', mode = "single")
rankNet(cellchat, mode = "single")


ranked_pathways <- rankNet(cellchat, measure  = 'count', mode = "single")
ranked_pathways <- rankNet(cellchat, mode = "single", return.data = TRUE)

df <- ranked_pathways$signaling.contribution
head(df)
top10 <- df[order(-df$contribution.scaled), ][1:10, ]
print(top10)
ggplot(top10, aes(x = reorder(name, contribution.scaled), y = contribution.scaled)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Pathway", y = "Scaled Contribution", title = "Top 10 Pathways") +
  theme_minimal()


runCellChatApp(cellchat)

#saveRDS(cellchat, file = '/rep 3/cellchatd0d3/groups3samename8umRegionsD0cellchat_trim001.200.20.rds')
#saveRDS(cellchat, file = '/rep 3/cellchatd0d3/groups3samename8umRegionsD3cellchat_trim001.200.20.rds')

#saveRDS(cellchat, file = '/rep 3/cellchatd0d3/12groups3dbcellchat_trim001.200.20.rds')
cellchat <- readRDS('/rep 3/cellchatd0d3/12groups3dbcellchat_trim001.200.20.rds')
cellchat <- updateClusterLabels(cellchat, new.order = c('D0 Crypt stroma', 'D0 Crypt', 'D0 Junction', 'D0 Villus stroma', 'D0 Villus epithelium', 'D0 Tip', 'D3 Crypt Stroma', 'D3 Crypt', 'D3 Junction', 'D3 Villus stroma', 'D3 Villus epithelium', 'D3 Tip'))
cellchat <- filterCommunication(cellchat, min.cells = 10) #minimum number of cells required in each cell group for cell-cell communication
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')

cellchat.d0 <- readRDS('/rep 3/cellchatd0d3/groups3samename8umRegionsD0cellchat_trim001.200.20.rds')
cellchat.d3 <- readRDS('/rep 3/cellchatd0d3/groups3samename8umRegionsD3cellchat_trim001.200.20.rds')
cellchat.d0 <- updateClusterLabels(cellchat.d0, new.order = c("Crypt Stroma", 'Crypt', 'Junction', 'Villus stroma', "Villus epithelium", "Tip"))
cellchat.d3 <- updateClusterLabels(cellchat.d3, new.order = c("Crypt Stroma", 'Crypt', 'Junction', 'Villus stroma', "Villus epithelium", "Tip"))



cellchat.d0 <- filterCommunication(cellchat.d0, min.cells = 10) #minimum number of cells required in each cell group for cell-cell communication
cellchat.d3 <- filterCommunication(cellchat.d3, min.cells = 10)
cellchat.d0 <- computeCommunProbPathway(cellchat.d0)
cellchat.d0 <- aggregateNet(cellchat.d0)
cellchat.d0 <- netAnalysis_computeCentrality(cellchat.d0, slot.name = 'netP')
cellchat.d3 <- computeCommunProbPathway(cellchat.d3)
cellchat.d3 <- aggregateNet(cellchat.d3)
cellchat.d3 <- netAnalysis_computeCentrality(cellchat.d3, slot.name = 'netP')

#write.csv(cellchat.d0@netP$prob, file = "rep3d0netp.csv", row.names = TRUE)
#write.csv(cellchat.d3@netP$prob, file = "rep3d3netp.csv", row.names = TRUE)


