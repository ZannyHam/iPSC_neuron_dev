#Setup
TASKDIR = "/net/bbi/vol1/data/aham7/nobackup/containers/analysis/"
setwd("/net/bbi/vol1/data/aham7/nobackup/containers/analysis/")

#load some packages
library(ggplot2)
library(tidyverse)
library(viridis)
library(scales)
library(Matrix)
library(Seurat)
library(dplyr)
library(presto)
library(MAST)

#####load in the CRISPR screen file and filter

crisprscreen <- read.table("neuro_CRISPR_scores.csv", sep=",", header=T)
crisprscreen$X <- NULL #drop column, there is nothing there

## make volcano plot
ggplot(crisprscreen, aes(x=epsilon, y=-log10(pvalue)))+
  geom_point(shape=21)+
  theme_light()

#only the unclear genes that did not significantly decrease cell #s
filt.crisprscreen <- filter(crisprscreen, gene != "NTC", epsilon >= -2) 

## Extract the gene names from each Seurat object
filt.crisprscreen_genes <- filt.crisprscreen$index
head(filt.crisprscreen_genes)

##### read in and plot flo+bt data

#Load Flo's visualization Seurat object - Barbara's group's data + 5K of Flo's
flobtneurons <- readRDS("processed_seurat_obj_with_gRNA_data.rds")

#remove flo's 5000 cells & tiny clusters & final time point
missing_clusters <- is.na(flobtneurons@meta.data$seurat_clusters)
btneurons <- subset(flobtneurons, cells = colnames(flobtneurons)[!missing_clusters])

# change the default assay
Assays(btneurons)
DefaultAssay(btneurons) <- "RNA"
dim(btneurons)

Idents(btneurons) <- "seurat_clusters"
Idents(btneurons)
btneurons <- subset(btneurons, idents = setdiff(levels(Idents(btneurons)), c("0", "6", "9")))
unique(btneurons$seurat_clusters)

Idents(btneurons) <- "timepoint"
Idents(btneurons)
btneurons <- subset(btneurons, idents = setdiff(levels(Idents(btneurons)), c("w5")))
unique(btneurons$timepoint)
Idents(btneurons) <- "seurat_clusters"

#calculate average expression for heatmap
average_expression <- AverageExpression(btneurons,
                                        group.by = "seurat_clusters", 
                                        assays = "RNA",
                                        slot = "data",
                                        normalization.method = "LogNormalize"
)
dim(average_expression[[1]])
head(average_expression[[1]])

# Rotate UMAP by 270 degrees (270 clockwise)
rotate_coords <- function(coords) {
  rotation_matrix <- matrix(c(0.7071, -0.7071,
                              0.7071, 0.7071), ncol = 2)  # 270-degree rotation matrix
  as.data.frame(as.matrix(coords) %*% rotation_matrix)
}

# Extract UMAP coordinates and rotate
umap_coords <- as.data.frame(Embeddings(btneurons, reduction = "umap"))
rotated_umap <- rotate_coords(umap_coords[, 1:2])
colnames(rotated_umap) <- c("UMAP_1", "UMAP_2")

# Save the rotated coordinates back to the Seurat object
btneurons@reductions$umap@cell.embeddings <- as.matrix(rotated_umap)

Idents(btneurons) <- "timepoint"
dimplot <- DimPlot(btneurons, reduction = "umap", group.by = "timepoint") +
  scale_color_manual(values = c("d1" = "purple", "d2" = "turquoise", "d5" = "gold", "w2" = "orange", "w4" = "red"),
                                        labels = c("d0", "d1", "d5", "d14", "d28"))
print(dimplot)

DimPlot(btneurons, reduction = "umap", group.by = c("seurat_clusters")) 
 
Idents(btneurons) <- "seurat_clusters" 

#change seurat_clusters from an integer to a factor
class(btneurons$seurat_clusters)
btneurons$seurat_clusters <- as.factor(btneurons$seurat_clusters)
class(btneurons$seurat_clusters)
Layers(btneurons)

FeaturePlot(btneurons, features=c("NR2F1"))

#find the marker genes in the Treutlein data, aka gene expression that is specific to clusters and not general
bt_allmarkers.clust.wilcox <- FindAllMarkers(btneurons,
                                             assay = 'RNA',
                                             group.by = "seurat_clusters", #this directs which column the function will use to group the expression data
                                             min.pct=0.1,
                                             logfc.threshold=0.25, # could go up to 0.5 if you want a more stringent differentiation threshold.
                                             #test.use="MAST", #default is wilcox, can try both MAST or wilcox.
                                             #latent.vars=c("batch") #use for MAST when there are technical or biological replicates that could confound results.
)

#make 1gene=1row and have columns of p values for each cluster
testdat <- bt_allmarkers.clust.wilcox[,c('gene', 'cluster', 'p_val_adj')]
dim(testdat)
head(bt_allmarkers.clust.wilcox)
testdat <- testdat %>%
  filter(p_val_adj <= 0.05)
dim(testdat)
head(testdat)
btneurons.markers <- spread(testdat, cluster, p_val_adj)
dim(btneurons.markers)
head(btneurons.markers)
length(btneurons.markers$gene) #should be <=total gene count in original object

#cross btneurons and the CRISPR knockout gene list
btneuronsxcrispr <- intersect(btneurons.markers$gene, filt.crisprscreen_genes)
head(btneuronsxcrispr)
length(btneuronsxcrispr)

#rename columns and filter expression
master <- data.frame(average_expression)
master$gene <- rownames(master)
master_filtered <- master[master$gene %in% unlist(btneuronsxcrispr), ]
colnames(master_filtered) <- c("1", "2", "3", "4", "5", "7", "8", "10", "11", "12", "13", "14", "15", "16", "17", "gene")
head(master_filtered)
master_filtered$num_clusters_gt2 <- rowSums(master_filtered[, 1:(ncol(master_filtered) - 1)] > 2)

unique(master_filtered$num_clusters_gt2)
length(master_filtered$num_clusters_gt2)
head(master_filtered)

filtered_genes <- master_filtered %>%
  filter(num_clusters_gt2 == 1) %>%  # changeable 
  pull(gene)

length(filtered_genes)
head(filtered_genes)
Idents(btneurons) <- factor(as.numeric(as.character(Idents(btneurons))))
DoHeatmap(btneurons, features = filtered_genes, slot = "data", assay='RNA')
DotPlot(btneurons, features=filtered_genes, assay='RNA')

#save gene list - change name when making new lists
write.csv(filtered_genes, file = "exp>2clust=1_n35", row.names = FALSE)
