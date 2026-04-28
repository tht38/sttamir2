library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(bsicons)

seurat_obj <- Load10X_Spatial(data.dir = '/rep 2/Jej_1B5_abc15_keep/8um')

metadata <- read.csv('/rep 2/Jej_1B5_abc15_keep/8um/FullSamples.NoC2 (1).csv', row.names = 1)
metadata <- metadata[rownames(metadata) %in% colnames(seurat_obj), , drop = FALSE]
seurat_obj <- AddMetaData(seurat_obj, metadata)

metadata <- read.csv('/rep 2/Jej_1B5_abc15_keep/8um/Graph-Based (1).csv', row.names = 1)
metadata <- metadata[rownames(metadata) %in% colnames(seurat_obj), , drop = FALSE]
seurat_obj <- AddMetaData(seurat_obj, metadata)

seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", 
                            scale.factor = 10000)

seurat_obj <- subset(seurat_obj, subset = FullSamples.NoC2 != "")
seurat_obj <- subset(seurat_obj, subset = Graph.based != "")

#saveRDS(seurat_obj, file = '/rep 2/Jej_1B5_abc15_keep/8um/load10xspatial/8umcreatedobjnorm.rds')
seurat_obj <- LoadSeuratRds('/rep 2/Jej_1B5_abc15_keep/8um/load10xspatial/8umcreatedobjnorm.rds')

metadata <- read.csv('rep 2/Jej_1B5_abc15_keep/8um/25.4.8 villus epi and stroma D0 and D3.csv', row.names = 1)
metadata <- metadata[rownames(metadata) %in% colnames(seurat_obj), , drop = FALSE]
seurat_obj <- AddMetaData(seurat_obj, metadata)
Idents(seurat_obj) <- seurat_obj$Stromal.regions



regions_to_keep <- c("Day 0 crypt", "Day 0 tip C4", "Day 0 Mature epi C1", "Day 0 Young villus C8", "Day 0 junction C7", "D0 Villus stroma", "D0 Intercrypt", "D0 Subcrypt", "Day 3 Crypt C5", "Day 3 intercrypt", "Day 3 Junction C10", "Day 3 subcrypt", "Day 3 tip", "day 3 villus", "Day 3 Villus epi C13")
seurat_subset <- subset(seurat_obj, subset = Stromal.regions %in% regions_to_keep)
table(seurat_subset$Stromal.regions)


seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 crypt"] <- "D0 Crypt"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 tip C4"] <- "D0 Tip"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 Mature epi C1"] <- "D0 Villus epithelium"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 Young villus C8"] <- "D0 Villus epithelium"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 junction C7"] <- "D0 Junction"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "D0 Villus stroma"] <- "D0 Villus stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "D0 Intercrypt"] <- "D0 Intercrypt"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "D0 Subcrypt"] <- "D0 Subcrypt"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 Villus stroma"] <- "D3 Villus stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "day 3 villus"] <- "D3 Villus stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 Crypt C5"] <- "D3 Crypt"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 intercrypt"] <- "D3 Intercrypt"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 Junction C10"] <- "D3 Junction"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 subcrypt"] <- "D3 Subcrypt"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 tip"] <- "D3 Tip"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 Villus epi C13"] <- "D3 Villus epithelium"
table(seurat_subset$Stromal.regions)
Idents(seurat_subset) <- seurat_subset$Stromal.regions


seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 crypt"] <- "D0 Crypt"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 tip C4"] <- "D0 Tip"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 Mature epi C1"] <- "D0 Villus epithelium"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 Young villus C8"] <- "D0 Villus epithelium"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 0 junction C7"] <- "D0 Junction"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "D0 Villus stroma"] <- "D0 Villus stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "D0 Intercrypt"] <- "D0 Crypt stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "D0 Subcrypt"] <- "D0 Crypt stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 Villus stroma"] <- "D3 Villus stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "day 3 villus"] <- "D3 Villus stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 Crypt C5"] <- "D3 Crypt"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 intercrypt"] <- "D3 Crypt stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 Junction C10"] <- "D3 Junction"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 subcrypt"] <- "D3 Crypt stroma"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 tip"] <- "D3 Tip"
seurat_subset$Stromal.regions[seurat_subset$Stromal.regions == "Day 3 Villus epi C13"] <- "D3 Villus epithelium"
table(seurat_subset$Stromal.regions)
Idents(seurat_subset) <- seurat_subset$Stromal.regions



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

scalefactors = jsonlite::fromJSON(txt = file.path('/rep 2/Jej_1B5_abc15_keep/8um/spatial', 'scalefactors_json.json'))
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
CellChatDB.use <- subsetDB(CellChatDB, search = "Non-protein Signaling", key = "annotation") # use Secreted Signaling

CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", 'ECM-Receptor', 'Cell-Cell Contact', 'Non-protein Signaling'), key = "annotation") # use Secreted Signaling

# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling) that can be only estimated from gene expression data. 
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", 'ECM-Receptor', 'Cell-Cell Contact'), key = "annotation") # use Secreted Signaling

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 7) # do parallel
future::plan("multicore", workers = 7)
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

#saveRDS(cellchat, file = '/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/npDBcombined12groups8umRegionscellchat_trim001.200.20.rds')
cellchat <- readRDS('/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/npDBcombined12groups8umRegionscellchat_trim001.200.20.rds')

#saveRDS(cellchat, file = '/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/combined15groups3samename8umRegionsD0cellchat_trim001.200.20.rds')
cellchat <- readRDS('/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/combined15groups3samename8umRegionsD0cellchat_trim001.200.20.rds')

#saveRDS(cellchat, file = '/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/combined14groups3samename8umRegionsD0cellchat_trim001.200.20.rds')
cellchat <- readRDS('/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/combined14groups3samename8umRegionsD0cellchat_trim001.200.20.rds')



#saveRDS(cellchat, file = '/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/4dbcombined12groups3samename8umRegionscellchat_trim001.200.20.rds')
cellchat <- readRDS('/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/4dbcombined12groups3samename8umRegionscellchat_trim001.200.20.rds')


#12 groups
#saveRDS(cellchat, file = '/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/combined12groups3samename8umRegionsD0cellchat_trim001.200.20.rds')
cellchat <- readRDS('/rep 2/Jej_1B5_abc15_keep/8um/cellchat.420/combined12groups3samename8umRegionsD0cellchat_trim001.200.20.rds')



# infers cell-cell communication at a signaling pathway level
cellchat <- filterCommunication(cellchat, min.cells = 10) #minimum number of cells required in each cell group for cell-cell communication
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')

cellchat <- updateClusterLabels(cellchat, old.cluster.name = c("Day 0 crypt", "Day 0 tip C4", "Day 0 Mature epi C1", "Day 0 Young villus C8", "Day 0 junction C7", "D0 Villus stroma", "D0 Intercrypt", "D0 Subcrypt", "Day 3 Crypt C5", "Day 3 intercrypt", "Day 3 Junction C10", "Day 3 subcrypt", "Day 3 tip", "day 3 villus", "Day 3 Villus epi C13"),
                                new.cluster.name = c('D0 Crypt', 'D0 Tip', 'D0 Mature epi', 'D0 Young epi', 'D0 Junction', 'D0 Villus stroma', 'D0 Intercrypt', 'D0 Subcrypt', 'D3 Crypt', 'D3 Intercrypt', 'D3 Junction', 'D3 Subcrypt', 'D3 Tip', 'D3 Villus stroma', 'D3 Villus epi'),
                                new.order = c('D0 Subcrypt', 'D3 Subcrypt', 'D0 Intercrypt', 'D3 Intercrypt', 'D0 Crypt', 'D3 Crypt', 'D0 Junction', 'D3 Junction', 'D0 Villus stroma', 'D3 Villus stroma', 'D0 Mature epi', 'D0 Young epi', 'D3 Villus epi', 'D0 Tip', 'D3 Tip'))

cellchat <- updateClusterLabels(cellchat, old.cluster.name = c("Day 0 crypt", "Day 0 tip C4", "Day 0 Mature epi C1", "Day 0 Young villus C8", "Day 0 junction C7", "D0 Villus stroma", "D0 Intercrypt", "D0 Subcrypt", "Day 3 Crypt C5", "Day 3 intercrypt", "Day 3 Junction C10", "Day 3 subcrypt", "Day 3 tip", "day 3 villus", "Day 3 Villus epi C13"),
                                new.cluster.name = c('D0 Crypt', 'D0 Tip', 'D0 Mature epi', 'D0 Young epi', 'D0 Junction', 'D0 Villus stroma', 'D0 Intercrypt', 'D0 Subcrypt', 'D3 Crypt', 'D3 Intercrypt', 'D3 Junction', 'D3 Subcrypt', 'D3 Tip', 'D3 Villus stroma', 'D3 Villus epi'),
                                new.order = c('D0 Subcrypt', 'D0 Intercrypt', 'D0 Crypt', 'D0 Junction', 'D0 Villus stroma', 'D0 Mature epi', 'D0 Young epi', 'D0 Tip', 'D3 Subcrypt', 'D3 Intercrypt', 'D3 Crypt', 'D3 Junction', 'D3 Villus stroma', 'D3 Villus epi', 'D3 Tip'))

cellchat <- updateClusterLabels(cellchat, new.order = c('D0 Subcrypt', 'D3 Subcrypt', 'D0 Intercrypt', 'D3 Intercrypt', 'D0 Crypt', 'D3 Crypt', 'D0 Junction', 'D3 Junction', 'D0 Villus stroma', 'D3 Villus stroma', 'D0 Villus epithelium', 'D3 Villus epithelium', 'D0 Tip', 'D3 Tip'))

cellchat <- updateClusterLabels(cellchat, new.order = c('D0 Subcrypt', 'D0 Intercrypt', 'D0 Crypt', 'D0 Junction', 'D0 Villus stroma', 'D0 Villus epithelium', 'D0 Tip', 'D3 Subcrypt', 'D3 Intercrypt', 'D3 Crypt', 'D3 Junction', 'D3 Villus stroma', 'D3 Villus epithelium', 'D3 Tip'))

#12 groups
cellchat <- updateClusterLabels(cellchat, new.order = c('D0 Crypt stroma', 'D0 Crypt', 'D0 Junction', 'D0 Villus stroma', 'D0 Villus epithelium', 'D0 Tip', 'D3 Crypt stroma', 'D3 Crypt', 'D3 Junction', 'D3 Villus stroma', 'D3 Villus epithelium', 'D3 Tip'))
celltype_colors <- c(
  "D0 crypt stroma" = "#d62728", "D0 Crypt" = "#aec7e8", "D0 Junction" = "#2ca02c", "D0 Villus stroma" = "#c5b0d5", "D0 villus epi" = "#F29401", "D0 Tip" = "#e377c2",
  "D3 crypt stroma" = "#d62728","D3 Crypt" = "#aec7e8","D3 Junction" = "#2ca02c","D3 Villus stroma" = "#c5b0d5","D3 Villus epi" = "#F29401","D3 Tip" = "#e377c2"
)

colormatch <- c(
  "D0 Subcrypt" = "#BF40BF",
  "D3 Subcrypt" = "#BF40BF",
  "D0 Intercrypt" = "#2ca02c",
  "D3 Intercrypt" = "#2ca02c",
  "D0 Crypt" = "#d62728",
  "D3 Crypt" = "#d62728",
  "D0 Junction" = "#e377c2",
  "D3 Junction" = "#e377c2",
  "D0 Villus stroma" = "#bcbd22",
  "D3 Villus stroma" = "#bcbd22",
  "D0 Mature epi" = "#aec7e8",
  "D0 Young epi" = "#aec7e8",
  "D3 Villus epi" = "#aec7e8",
  "D0 Tip" = "#c5b0d5",
  "D3 Tip" = "#c5b0d5"
)

colormatch <- c(
  "D0 Subcrypt" = "#bc9dcc",
  "D3 Subcrypt" = "#bc9dcc",
  "D0 Intercrypt" = "#f781bf",
  "D3 Intercrypt" = "#f781bf",
  "D0 Crypt" = "#e3221c",
  "D3 Crypt" = "#e3221c",
  "D0 Junction" = "#f29401",
  "D3 Junction" = "#f29401",
  "D0 Villus stroma" = "#984ea4",
  "D3 Villus stroma" = "#984ea4",
  "D0 Villus epithelium" = "#377eb8",
  "D3 Villus epithelium" = "#377eb8",
  "D0 Tip" = "#4eaf4a",
  "D3 Tip" = "#4eaf4a"
)

celltype_colors <- c(
  "D0 Subcrypt" = "#BF40BF",
  "D0 Intercrypt" = "#2ca02c",
  "D0 Crypt" = "#d62728",
  "D0 Junction" = "#e377c2",
  "D0 Villus stroma" = "#bcbd22",
  "D0 Mature epi" = "#aec7e8",
  "D0 Young epi" = "#aec7e8",
  "D0 Tip" = "#c5b0d5",
  "D3 Subcrypt" = "#BF40BF",
  "D3 Intercrypt" = "#2ca02c",
  "D3 Crypt" = "#d62728",
  "D3 Junction" = "#e377c2",
  "D3 Villus stroma" = "#bcbd22",
  "D3 Villus epi" = "#aec7e8",
  "D3 Tip" = "#c5b0d5"
)

timecolors <- c(
  "D0 Subcrypt" = "#bc9dcc",
  "D0 Intercrypt" = "#f781bf",
  "D0 Crypt" = "#e3221c",
  "D0 Junction" = "#f29401",
  "D0 Villus stroma" = "#984ea4",
  "D0 Villus epithelium" = "#377eb8",
  "D0 Tip" = "#4eaf4a",
  "D3 Subcrypt" = "#bc9dcc",
  "D3 Intercrypt" = "#f781bf",
  "D3 Crypt" = "#e3221c",
  "D3 Junction" = "#f29401",
  "D3 Villus stroma" = "#984ea4",
  "D3 Villus epithelium" = "#377eb8",
  "D3 Tip" = "#4eaf4a"
)

pathways.show = 'GUCA'
pathways.show = c('TGFb', 'MIF', 'NRG', 'COMPLEMENT', 'FGF', "EGF", "VEGF", "ANGPTL", "LAMININ", "CEACAM", "COLLAGEN", "TNF", 'GALECTIN', 'IL6', 'IGF')
pathways.show = c('TGFb', 'MIF', 'NRG', 'COMPLEMENT', "EGF", "VEGF", "ANGPTL", "LAMININ", "COLLAGEN", "TNF", 'GALECTIN', 'IL6', 'IGF')

netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", signaling = pathways.show, color.use = celltype_colors, width = 5, height = 6, color.heatmap = "OrRd")

netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = pathways.show, color.use = celltype_colors, width = 5, height = 6, color.heatmap = "OrRd")

pathways.show = c('TGFb', 'MIF', 'NRG', 'COMPLEMENT', "EGF", "VEGF", "ANGPTL", "TNF", 'GALECTIN', 'IL6', 'IGF', "COLLAGEN", 'TENASCIN', 'Prostaglandin', 'RA')
pathways.show = c('Prostaglandin', 'RA')
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = pathways.show, color.use = celltype_colors, width = 5, height = 6, color.heatmap = "OrRd")


pathways.show = c('TGFb', 'MIF', 'NRG', 'COMPLEMENT', "EGF", "VEGF", "ANGPTL", "TNF", 'GALECTIN', 'GUCA', 'IGF', "COLLAGEN", 'TENASCIN', 'Prostaglandin', 'RA')
