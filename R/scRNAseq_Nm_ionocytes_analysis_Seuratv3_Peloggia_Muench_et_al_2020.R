library(Seurat)
library(dplyr)
library(ggplot2)
library(WriteXLS)
library(stringr)
library(patchwork)
library(RColorBrewer)

# Configuring data and figure file paths and reading in gene information ====
dataPath <- function(filename){paste0("./data/", filename)}
figurePath <- function(filename){paste0("./figures/", filename)}

# Zebrafish gene info
gene_info <- read.delim(paste0("./Danio_Features_unique_Ens98_v1.tsv"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# Pre-processing and Quality Control ====
# Reading in the filtered matrix
filtered_mat <- Read10X("./zipped/")

# Initialize the Seurat object with the raw (non-normalized data).
seurat_obj <- CreateSeuratObject(counts = filtered_mat,
                                 project = "polar cell project", min.cells = 5)

# Adding percentage of mitochondrial genes to the metadata
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt")

# Checking distribution of genes, UMIs and mt genes for determining QC parameters
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# Print number of cells before applying QC cutoffs
n_cells1 <- dim(seurat_obj)[2]

# Remove cells with greater than 30 percent mitochondrial genes,
# more than 30.000 UMIs and less than 250 genes
seurat_obj <- subset(seurat_obj, subset = percent.mt < 30)
seurat_obj <- subset(seurat_obj, subset = nCount_RNA < 30000)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 250)

# Checking new data distribution after quality control
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# How many cells were removed with these parameters?
n_cells2 <- dim(seurat_obj)[2]
print(n_cells1 - n_cells2)

# Data normalization ====
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# Scaling data for all genes ====
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# Dimensional reduction ====

# Running principal components analysis
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Elbow plot for determining the dimensionality of the dataset. Calculating 50 dims
ElbowPlot(seurat_obj, ndims = 50)

# Clustering and UMAP
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:26)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:26)

# Looking at the UMAP for the first time
DimPlot(seurat_obj, reduction = "umap",  pt.size = 0.5, label = TRUE, label.size = 6)

# Saving file
saveRDS(seurat_obj, file = "./data/seurat_obj_new.rds")

# Finding cluster markers and exporting excel file with gene information ====

# Parallel solution for FindAllMarkers (from ddiaz - 1-15-20)
n_clust <- 1:length(levels(Idents(seurat_obj)))
clust <- levels(Idents(seurat_obj))

mcFindMarkers <- function(i) {
  ident1 <- clust[i]
  ident2 <- clust[clust != clust[i]]
  table <- FindMarkers(seurat_obj,
                       ident.1 = ident1, ident.2 = ident2, only.pos = FALSE)
  table$Gene.name.uniq <- rownames(table)
  table$cluster <- rep(clust[i], nrow(table))
  return(table)
}

marker_results <- list()[n_clust]
ptm <- proc.time()
marker_results <- parallel::mclapply(n_clust, mcFindMarkers, mc.cores = 8)
time_diff <- proc.time() - ptm
time_diff
markers <- marker_results

if(TRUE) {
  saveRDS(markers, dataPath("clustermarkers.rds"))
  markers <- readRDS(dataPath("clustermarkers.rds"))
}

markers <- dplyr::bind_rows(markers)
str(markers)
table(markers$cluster)

marker_subset <- markers[markers$p_val < 0.05,]
marker_subset$pct.ratio <- 1:nrow(marker_subset)
marker_subset$pct.ratio <- marker_subset$pct.1 / marker_subset$pct.2
marker_subset <- marker_subset[(order(marker_subset$cluster,
                                      -1 * (marker_subset$pct.ratio))),]


dim(marker_subset)
table(marker_subset$cluster)
marker_table <- inner_join(marker_subset, gene_info, by = "Gene.name.uniq")



WriteXLS::WriteXLS(marker_table,
                   figurePath("markers.xlsx"),
                   row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)

# Plotting mCherry+ cells per cluster ====

# Establishing a ggplot theme for the plot with desired characteristics
my_theme <- theme(axis.line.y = element_line(color = "black", size = 0.1),
                     axis.line.x = element_line(color = "black", size = 0.1),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA),
                     axis.title.x = element_text(size = 20, margin = margin(5,0,0,0)),
                     axis.title.y = element_text(size = 20, margin = margin(0,10,0,0)),
                     axis.text = element_text(size = 16, color = "black"),
                     axis.text.x = element_text(margin = margin(t=5)),
                     plot.title = element_text(size = 32, hjust = 0.5),
                     legend.position = c(0.9, 0.85),
                     legend.key.size = unit(1, "cm"),
                     legend.text = element_text(size = 20),
                     legend.text.align = 0,
                     legend.title = element_blank(),
                     legend.background = element_blank()
                     
)

# Subsets only mCherry+ cells
cherry <- subset(seurat_obj, subset = ENSmCherry > 0)

# This creates a data frame which contains the cluster IDs
# and the number of mCherry+ cells on each one of them
ncells <- as.data.frame(table(Idents(cherry)))

# Making a ggplot2 bar plot
ggplot(data=ncells, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill = "lightgrey", colour = "black", width = 0.8)  +
  # geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x = "Cluster Number", y ="Number of cells" ) +
  my_theme +
  geom_hline(yintercept=2, linetype="dotted", col = "red")
  
# Subsetting ionocytes and skin clusters ====

# This subsets NCC, HR, NaR and KS ionocyte clusters as well as tp63+ clusters
# into a new seurat object called "data"
data <- subset(seurat_obj, idents = c("0", "1", "4", "7", "10", "12"))

# Normalizing and scaling data again, as recommended on the Seurat tutorial following subsetting
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes, model.use = "negbinom")

# Re-running dimensionality reduction and clustering
data <- RunPCA(data, features = VariableFeatures(object = data),npcs = 50)
ElbowPlot(data, ndims = 50)
data <- FindNeighbors(data, dims = 1:25,force.recalc=T)
data <- FindClusters(data, resolution = 1,save.SNN = T,temp.file.location=".",force.recalc=T, dims=1:25)
data <- RunUMAP(data, dims = 1:25,force.recalc=T)

# New UMAP from the subset
DimPlot(data, reduction = "umap",label=T)

saveRDS(data, "./data/subset.rds")

# Renaming clusters ====
new.cluster.ids <- c("HR ionocytes", "NCC ionocytes", "NCC ionocytes", "skin", "NCC ionocytes",
                     "NaR ionocytes", "NCC ionocytes", "KS ionocytes",
                     "skin", "HR ionocytes", "NaR ionocytes", "ionocyte progenitors",
                     "skin", "skin", "skin", "skin", "KS ionocytes", "KS ionocytes", "ionocyte progenitors", "ionocyte progenitors")

names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)

saveRDS(data, "./data/subset_renamed-IDs.rds")

# Supervised clustering of notch1b+ and mCherry+ cells within the ionocyte clusters ====

ionocytes <- subset(data, idents = "skin", invert = TRUE)

subset <- subset(ionocytes, subset = notch1b > 0 | ENSmCherry > 0)

nmc <- WhichCells(subset)

Idents(data, cells = nmc) <- "notch1b+ cells"

clust <- Idents(data)

data <- AddMetaData(data, clust, col.name = "clusters")

ids <- c("notch1b+ cells", "HR ionocytes", "NaR ionocytes", "NCC ionocytes", "KS ionocytes", "ionocyte progenitors", "skin")

data@meta.data$clusters <- factor(
  data@meta.data$clusters, ordered = TRUE, levels = ids)


levels(data) <- c("notch1b+ cells", "HR ionocytes", "NaR ionocytes", "NCC ionocytes", "KS ionocytes", "ionocyte progenitors", "skin")

saveRDS(data, "./data/subset_renamed_with-notch1b_cells.rds")

# Finding markers for new model and exporting xlxs ====
n_clust <- 1:length(levels(Idents(data)))
clust <- levels(Idents(data))

mcFindMarkers <- function(i) {
  ident1 <- clust[i]
  ident2 <- clust[clust != clust[i]]
  table <- FindMarkers(data,
                       ident.1 = ident1, ident.2 = ident2, only.pos = FALSE)
  table$Gene.name.uniq <- rownames(table)
  table$cluster <- rep(clust[i], nrow(table))
  return(table)
}

marker_results <- list()[n_clust]
ptm <- proc.time()
marker_results <- parallel::mclapply(n_clust, mcFindMarkers, mc.cores = 2)
time_diff <- proc.time() - ptm
time_diff
markers <- marker_results

if(TRUE) {
  saveRDS(markers, dataPath("clustermarkers-subset.rds"))
  markers <- readRDS(dataPath("clustermarkers-subset.rds"))
}

markers <- dplyr::bind_rows(markers)
str(markers)
table(markers$cluster)

marker_subset <- markers[markers$p_val < 0.05,]
marker_subset$pct.ratio <- 1:nrow(marker_subset)
marker_subset$pct.ratio <- marker_subset$pct.1 / marker_subset$pct.2
marker_subset <- marker_subset[(order(marker_subset$cluster,
                                      -1 * (marker_subset$pct.ratio))),]


dim(marker_subset)
table(marker_subset$cluster)
marker_table <- inner_join(marker_subset, gene_info, by = "Gene.name.uniq")

WriteXLS::WriteXLS(marker_table,
                   dataPath("markers_subset.xlsx"),
                   row.names = FALSE, col.names = TRUE, AdjWidth = TRUE)



# Making paper dot plot and feature plots ====

# Features to plot on dot plot
to_plot <- c("ENSmCherry",
             "notch1b",
             "foxi3a",
             "foxi3b",
             "gcm2",
             "ca2",
             "ca15a",
             "rhcgb",
             "slc9a3.2",
             "atp1a1a.5",
             "trpv6",
             "atp2b2",
             "slc8a1b",
             "atp1a1a.2",
             "slc12a10.2",
             "clcn2c",
             "slc4a4b",
             "kcnj1a.3",
             "kcnj1a.5",
             "kcnj1a.6",
             "kcnj1a.1",
             "atp1b1b",
             "atp6v1aa",
             "atp6v1ab",
             "klf4",
             "krtt1c19e",
             "tp63",
             "krt4",
             "krt5"
)

# Dot plot for Figure 1P
DotPlot(data, features = to_plot, group.by = "clusters", cols = c("blue", "red")) + coord_flip() & RotatedAxis()

# Feature plots
FeaturePlot(data, features = "trpv6", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
FeaturePlot(data, features = "gcm2", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
FeaturePlot(data, features = "foxi3a", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
FeaturePlot(data, features = "slc12a.10.2", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
FeaturePlot(data, features = "CR936482.1", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
FeaturePlot(data, features = "atp1a1a.1", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
FeaturePlot(data, features = "foxi3b", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
FeaturePlot(data, features = "krtt1c19e", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
FeaturePlot(data, features = "tp63", cols = c("lightgrey", "midnightblue"), min.cutoff = "q10", max.cutoff = "q90") + NoAxes()

# Averaging expression for pseudo-bulk heatmap ====

cluster.averages <- AverageExpression(data)
head(cluster.averages[["RNA"]][, 1:5])

orig.levels <- levels(data)
Idents(data) <- gsub(pattern = " ", replacement = "_", x = Idents(data))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(data) <- orig.levels
cluster.averages <- AverageExpression(data, return.seurat = TRUE)
cluster.averages


cluster.markers <- FindAllMarkers(data, min.pct = 0.25)

clust <- Idents(cluster.averages)
cluster.averages <- AddMetaData(cluster.averages, clust, col.name = "clusters")

ids <- c("notch1b+ cells", "HR ionocytes", "NaR ionocytes", "NCC ionocytes", "KS ionocytes", "ionocyte progenitors", "skin")

cluster.averages@meta.data$clusters <- factor(
  cluster.averages@meta.data$clusters, ordered = TRUE, levels = ids)

levels(cluster.averages) <- c("notch1b+ cells", "HR ionocytes", "NaR ionocytes", "NCC ionocytes", "KS ionocytes", "ionocyte progenitors", "skin")

genes <- scan(paste0("genes.txt"), what = character(), sep = "\n")

DoHeatmap(cluster.averages, features = genes) + scale_fill_gradient2(
  low = rev(c('#D1E5F0','#67A9CF','#2166AC')),
  mid = "white", high = rev(c('#B2182B','#EF8A62','#FDDBC7')),
  midpoint = 0, guide = "colourbar", aesthetics = "fill")






