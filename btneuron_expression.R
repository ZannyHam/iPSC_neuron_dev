#Setup
TASKDIR = "/net/bbi/vol1/data/aham7/nobackup/containers/analysis/"
setwd("/net/bbi/vol1/data/aham7/nobackup/containers/analysis/")

#load some packages
library(monocle3)
library(ggplot2)
library(tidyverse)
library(viridis)
library(data.table)
library(scuttle)
library(cowplot)
library(scales)
library(stringr)
library(Matrix)
library(Seurat)
library(dplyr)
library(presto)
library(devtools)
library(MAST)
library(reshape2)

### read in and plot flo+bt data

#Load Flo's visualization Seurat object - Barbara's group's data + 5K of Flo's
flobtneurons <- readRDS("processed_seurat_obj_with_gRNA_data.rds")
flobtneurons

#Plotting UMAP with Seurat
DimPlot(flobtneurons, reduction = "umap", group.by = c("timepoint")) 
DimPlot(flobtneurons, reduction = "umap", group.by = c("batch"))
DimPlot(flobtneurons, reduction = "umap", group.by = c("seurat_clusters"))

#remove flo's 5000 cells as they are not assigned to clusters
missing_clusters <- is.na(flobtneurons@meta.data$seurat_clusters)
sum(missing_clusters)  # Count cells with missing cluster information
btneurons <- subset(flobtneurons, cells = colnames(flobtneurons)[!missing_clusters])
average_expression <- AggregateExpression(btneurons, group.by = "seurat_clusters", assays = "RNA")
dim(average_expression[[1]])
class(average_expression[[1]])
head(average_expression[[1]])
head(btneurons)

# change the default assay
Assays(btneurons)
DefaultAssay(btneurons) <- "RNA"

#sce <- as.SingleCellExperiment(btneurons)
#sce
#names(btneurons@reductions)
#umap_coords <- Embeddings(btneurons, reduction = "umap")
#umap_coords <- umap_coords[colnames(sce),]
#sce$UMAP1 <- umap_coords[,1]
#sce$UMAP2 <- umap_coords[,2]
#ggplot(data.frame(colData(sce)), aes(x=UMAP1, y=UMAP2, color=timepoint)) +
  #geom_point(size=0.5, alpha=0.3) +
  #theme_light()
#ggplot(data.frame(colData(sce)), aes(x=UMAP1, y=UMAP2, color=seurat_clusters)) +
  #geom_point(size=0.5, alpha=0.3) +
  #theme_light()

#take a look
DimPlot(btneurons, reduction = "umap", group.by = c("timepoint")) 
DimPlot(btneurons, reduction = "umap", group.by = c("batch"))
DimPlot(btneurons, reduction = "umap", group.by = c("seurat_clusters"))

#change seurat_clusters from an integer to a factor
class(btneurons$seurat_clusters)
btneurons$seurat_clusters <- as.factor(btneurons$seurat_clusters)
class(btneurons$seurat_clusters)

Idents(btneurons) <- "seurat_clusters"
Idents(btneurons)

#find the marker genes in the Treutlein data, aka gene expression that is specific to clusters and not general
bt_allmarkers.clust.wilcox <- FindAllMarkers(btneurons,
          assay = 'RNA',
          group.by = "seurat_clusters", #this directs which column the function will use to group the expression data
          min.pct=0.1,
          logfc.threshold=0.25, # could go up to 0.5 if you want a more stringent differentiation threshold.
          #test.use="wilcox", #default is wilcox, can try both MAST or wilcox.
          #latent.vars=c("batch") #use for MAST when there are technical or biological replicates that could confound results.
)

head(bt_allmarkers.clust.wilcox)
dim(bt_allmarkers.clust.wilcox) #the test split genes out in every cluster

#make 1gene=1row and have columns of p values for each cluster
testdat <- bt_allmarkers.clust.wilcox[,c('gene', 'cluster', 'p_val_adj')]
head(testdat)
btneurons.markers <- spread(testdat, cluster, p_val_adj)
dim(btneurons.markers)
head(btneurons.markers)
length(btneurons.markers$gene) #should be <=total gene count in original object

#### Make a doc with the marker gene expression by "seurat_clusters" for btneurons
write.csv(btneurons.markers, file = "bt_marker_gene_expression_by_cluster_minpct0_1_thresh0_25_wilcox_seurat.csv", row.names = FALSE)

###get averages of expression data for heatmaps and dot plots
average_expression <- AverageExpression(btneurons,
                      group.by = "seurat_clusters", 
                      assays = "RNA",
                      slot = "data",
                      normalization.method = "LogNormalize"
                      )

#master <- data.frame(average_expression)
#master$gene <- rownames(master)
#master1 <- left_join(master, bt_allmarkers.clust.wilcox)
#master <- left_join(master, btneurons.markers)


######## Make a list of genes involved with neurological disorders from Shawn's CRISPR knockout screen #####

###load in the CRISPR screen file and filter
crisprscreen <- read.table("neuro_CRISPR_scores.csv", sep=",", header=T)
table(is.na(crisprscreen$X))
crisprscreen$X <- NULL #drop column, there is nothing there
table(crisprscreen$index)
filt.crisprscreen <- filter(crisprscreen, gene != "NTC") %>%
  filter(epsilon >= -2 & epsilon <= 2) #only the unclear genes that did not significantly increase or decrease cell #s

### Extract the gene names from each Seurat object
filt.crisprscreen_genes <- filt.crisprscreen$index

#identify genes involved with neuronal disease phenotypes - data from GenCC
all_neuro_dis.all <- read.csv("gencc_allneurological.csv", header=TRUE) #pre-filtered for definitive-strong connections to neurological diseases
all_neuro_dis.all <- all_neuro_dis.all$gene_symbol #just get a list of gene names
all_neuro_dis <- unique(all_neuro_dis.all) #remove all the duplicates

#load in essential genes lists (from Shawn's papers)
hesc_ess <- read.csv("hesc_essentialome.csv", header=TRUE)
hesc_ess <- hesc_ess$Gene.symbol
hap_ess <- read.csv("hap_ess.csv", header=TRUE)
hap_ess <- hap_ess$GENE_SYMBOL

#make a list of only relevant neurological disease genes in Shawn's CRISPR knockout list
filt.all_neuro_crispr <- intersect(filt.crisprscreen_genes, all_neuro_dis) 

#filter out essential genes
filt.all_neuro_crispr <- setdiff(filt.all_neuro_crispr, hesc_ess)
filt.all_neuro_crispr <- setdiff(filt.all_neuro_crispr, hap_ess)
length(filt.all_neuro_crispr) #goal is ~250 genes

#save the list!
write.csv(filt.all_neuro_crispr, file = "final_gene_list_neuron.csv", row.names = FALSE)

##### Filter btneurons marker genes to only have genes in my list #####

# Get the gene names in the btneurons Seurat object
genes_in_btneurons <- rownames(btneurons@assays$RNA)
length(genes_in_btneurons)

# Identify genes that are in both btneurons and filt.all_neuro_crispr
matching_genes <- intersect(genes_in_btneurons, filt.all_neuro_crispr)
length(matching_genes)
head(matching_genes)

#filter marker genes to only have genes in this list
matching_markers <- intersect(btneurons.markers$gene, matching_genes)
length(matching_markers)
head(matching_markers)

#make plot of only marker genes
DotPlot(btneurons, features=matching_markers) #this is all the genes so we need to narrow it down
#DoHeatmap(btneurons, features=matching_markers)


##### need to check everything below here

#filter markers to only matching ones in the statistical test and av expression
markers <- filter(bt_allmarkers.clust.wilcox, gene %in% matching_markers) 
markers.avg <- average_expression[[1]][matching_markers,]
length(unique(markers$gene)) #266 unique genes

#compare lists of marker genes and the overall object
markergenes <- bt_allmarkers.clust.wilcox$gene 
genesinseurat <- rownames(btneurons) 
length(markergenes)
length(genesinseurat)
length(intersect(markergenes, genesinseurat))
setdiff(genesinseurat, markergenes)

####

#look at how the marker genes are expressed
hist(markers$avg_log2FC)
hist(-log10(markers$p_val_adj))

#filter the markers to be statistically sig and have a minimal threshold for fold change
markers <- filter(markers, p_val_adj < 0.05 & avg_log2FC > 0.025)
markers <- arrange(markers, desc(avg_log2FC))

#make the average expression into a data frame and organize by gene name
avgexp <- data.frame(average_expression[[1]])
avgexp$gene <- rownames(avgexp)

#reshape the table to focus on genes, clusters, and expression values
long <- melt(avgexp, id.vars=c("gene"))
long$variable <- gsub("g", "", long$variable)
table(long$variable)
table(markers$cluster)

#add average expression column, look at data and organize
markers$avg_exp <- long[match(paste(markers$gene, markers$cluster), 
                      paste(long$gene, long$variable)), 'value']
hist(markers$avg_exp)
hist(markers$pct.1)
markers <- arrange(markers, desc(avg_exp))

ggplot(markers, aes(x=pct.1, y=log2(avg_exp +1), color=p_val_adj)) +
  geom_point() +
  scale_color_viridis()

arrange(markers, p_val_adj) %>% head()

VlnPlot(btneurons, features = c("FTL"))
FeaturePlot(btneurons, features=c("NEFL"))
FeaturePlot(btneurons, features=c("VAMP2"))

# Dimensions of the raw counts matrix
dim(btneurons@assays$RNA@counts)
btneurons@assays$RNA@counts["NEFL",]

# Dimensions of the log-normalized data
dim(btneurons@assays$RNA@data)
btneurons@assays$RNA@data["NEFL",]

dim(btneurons@assays$integrated@data)
btneurons@assays$integrated@data["NEFL",]

# Dimensions of the scaled data
dim(btneurons@assays$integrated@scale.data)
btneurons@assays$integrated@scale.data["NEFL",]


ggplot() +
  geom_point(aes(x=btneurons@assays$RNA@counts["NEFL",], y=btneurons@assays$RNA@data["NEFL",]
))

ggplot() +
  geom_point(aes(x=btneurons@assays$RNA@data["NEFL",], btneurons@assays$integrated@data["NEFL",]
  ))
