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
library(scran)
library(SeuratDisk)
library(SingleCellExperiment)

###load in the CRISPR screen file and filter
crisprscreen <- read.table("neuro_CRISPR_scores.csv", sep=",", header=T)
table(is.na(crisprscreen$X))
crisprscreen$X <- NULL #drop column, there is nothing there
table(crisprscreen$index)

## make volcano plot
ggplot(crisprscreen, aes(x=epsilon, y=-log10(pvalue)))+
  geom_point(shape=21)+
  theme_light()

filt.crisprscreen <- filter(crisprscreen, gene != "NTC") %>%
  filter(epsilon >= -2) #only the unclear genes that did not significantly decrease cell #s

## Extract the gene names from each Seurat object
filt.crisprscreen_genes <- filt.crisprscreen$index
head(filt.crisprscreen_genes)

### read in and plot flo+bt data

#Load Flo's visualization Seurat object - Barbara's group's data + 5K of Flo's
flobtneurons <- readRDS("processed_seurat_obj_with_gRNA_data.rds")
flobtneurons

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

DimPlot(btneurons, reduction = "umap", group.by = c("seurat_clusters"))
DimPlot(btneurons, reduction = "umap", group.by = c("timepoint"))

#change seurat_clusters from an integer to a factor
class(btneurons$seurat_clusters)
btneurons$seurat_clusters <- as.factor(btneurons$seurat_clusters)
class(btneurons$seurat_clusters)

Idents(btneurons) <- "seurat_clusters"
Idents(btneurons)

Layers(btneurons)

#find highly variable genes based on gene expression and variance (deconvolution)
dec.neurons <- modelGeneVar(sce)

# Visualizing the fit:
fit.neurons <- metadata(dec.neurons)
plot(fit.neurons$mean, fit.neurons$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.neurons$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# Ordering by most interesting genes for inspection.
dec.neurons[order(dec.neurons$bio, decreasing=TRUE),] 
head(dec.neurons)

set.seed(0010101)
dec.pois.neurons <- modelGeneVarByPoisson(sce)
dec.pois.neurons <- dec.pois.neurons[order(dec.pois.neurons$bio, decreasing=TRUE),]
head(dec.pois.neurons)

plot(dec.pois.neurons$mean, dec.pois.neurons$total, pch=16, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(metadata(dec.pois.neurons)$trend(x), col="dodgerblue", add=TRUE)

dec.cv2.sce <- modelGeneCV2(sce)
fit.cv2.sce <- metadata(dec.cv2.sce)
plot(fit.cv2.sce$mean, fit.cv2.sce$cv2, log="xy")
curve(fit.cv2.sce$trend(x), col="dodgerblue", add=TRUE, lwd=2)
dec.cv2.sce[order(dec.cv2.sce$ratio, decreasing=TRUE),]
hvg.sce.cv2 <- getTopHVGs(dec.cv2.sce, var.field="ratio", n=1000)
str(hvg.sce.cv2)
dec.sce <- modelGeneVar(sce) v 
chosen <- getTopHVGs(dec.sce, prop=0.2)
str(chosen)

crispr.hvg <- intersect(filt.crisprscreen_genes, chosen)
head(crispr.hvg)
length(crispr.hvg)

#find the marker genes in the Treutlein data, aka gene expression that is specific to clusters and not general
bt_allmarkers.clust.wilcox <- FindAllMarkers(btneurons,
                                             assay = 'RNA',
                                             group.by = "seurat_clusters", #this directs which column the function will use to group the expression data
                                             min.pct=0.1,
                                             logfc.threshold=0.25, # could go up to 0.5 if you want a more stringent differentiation threshold.
                                             test.use="MAST", #default is wilcox, can try both MAST or wilcox.
                                             latent.vars=c("batch") #use for MAST when there are technical or biological replicates that could confound results.
)

# Get cell and feature names, and total numbers We show multiple ways to get the same output cell names


head(bt_allmarkers.clust.wilcox)
dim(bt_allmarkers.clust.wilcox) #the test split genes out in every cluster

#make 1gene=1row and have columns of p values for each cluster
testdat <- bt_allmarkers.clust.wilcox[,c('gene', 'cluster', 'p_val_adj')]
head(testdat)
btneurons.markers <- spread(testdat, cluster, p_val_adj)
dim(btneurons.markers)
head(btneurons.markers)
length(btneurons.markers$gene) #should be <=total gene count in original object


btneuronsxcrispr <- intersect(btneurons.markers$gene, filt.crisprscreen_genes)
head(btneuronsxcrispr)
length(btneuronsxcrispr)
hvg.crispr.btneurons <- intersect(btneurons.markers$gene, crispr.hvg)
head(hvg.crispr.btneurons)
length(hvg.crispr.btneurons)

# Filter btneurons.markers to include only genes in btneuronsxcrispr
filtered_markers <- btneurons.markers %>%
  filter(gene %in% btneuronsxcrispr)

filt.hvg.crispr.btneurons <- btneurons.markers %>%
  filter(gene %in% hvg.crispr.btneurons)

# View the filtered table
head(filtered_markers)
dim(filtered_markers) #should match dim of btneuronsxcrispr

head(filt.hvg.crispr.btneurons)
dim(filt.hvg.crispr.btneurons)

DoHeatmap(btneurons, features=filt.hvg.crispr.btneurons$gene, slot='data', assay='RNA')

#scale down list to only genes that are expressed in one cluster
# Count the number of clusters where each gene is expressed
filtered_unique <- filtered_markers %>%
  rowwise() %>%                                    # Operate on rows
  mutate(num_clusters = sum(!is.na(c_across(-gene)))) %>%  # Count non-NA clusters (exclude the gene column)
  filter(num_clusters == 1) %>%                    # Keep genes expressed in exactly one cluster
  select(-num_clusters)                            # Drop the helper column

filtered_unique.hvg <- filt.hvg.crispr.btneurons %>%
  rowwise() %>%                                    # Operate on rows
  mutate(num_clusters = sum(!is.na(c_across(-gene)))) %>%  # Count non-NA clusters (exclude the gene column)
  filter(num_clusters == 1) %>%               # Keep genes expressed in specified # of clusters
  select(-num_clusters)                            # Drop the helper column

# View the filtered table
head(filtered_unique)
dim(filtered_unique)

head(filtered_unique.hvg)
dim(filtered_unique.hvg)

DoHeatmap(btneurons, features=filtered_unique$gene, slot='data', assay='RNA')
DoHeatmap(btneurons, features=filtered_unique.hvg$gene, slot='data', assay='RNA', disp.min=0, disp.max=2)

##### try filtering btneurons based on highly expressed RNA #####

# Calculate total expression for each gene across all cells
total_expression <- Matrix::rowSums(btneurons@assays$RNA@counts)

# Calculate average expression for each gene
average_expression <- Matrix::rowMeans(btneurons@assays$RNA@counts)
# Define the number of top genes to retain
top_n <- 1000

# Get the names of the top genes based on total expression
top_genes <- names(sort(total_expression, decreasing = TRUE))[1:top_n]
top_genes <- names(sort(average_expression, decreasing = TRUE))[1:top_n]
# Subset the Seurat object to include only the top genes
btneurons_filtered <- subset(btneurons, features = top_genes)
dim(btneurons_filtered@assays$RNA@counts)  # Should show top_n rows (genes)

# Create a heatmap for the top genes
DoHeatmap(btneurons_filtered, features = top_genes, slot = "data")

#bt_topmarkers.clust.wilcox <- FindAllMarkers(btneurons_filtered,
                                            # assay = 'RNA',
                                            # group.by = "seurat_clusters", #this directs which column the function will use to group the expression data
                                            # min.pct=0.1,
                                           #  logfc.threshold=0.25, # could go up to 0.5 if you want a more stringent differentiation threshold.
                                             #test.use="wilcox", #default is wilcox, can try both MAST or wilcox.
                                             #latent.vars=c("batch") #use for MAST when there are technical or biological replicates that could confound results.
#)

shorter <- bt_allmarkers.clust.wilcox[,c('gene', 'cluster', 'p_val_adj')]
head(shorter)
btneurons.topmarkers <- spread(shorter, cluster, p_val_adj)
dim(btneurons.topmarkers)
head(btneurons.topmarkers)
length(btneurons.topmarkers$gene) #should be <=total gene count in original object

btneuronstopxcrispr <- intersect(btneurons.topmarkers$gene, filt.crisprscreen_genes)
head(btneuronstopxcrispr)
length(btneuronstopxcrispr)

# Filter btneurons.markers to include only genes in btneuronsxcrispr
filtered_topmarkers <- btneurons.topmarkers %>%
  filter(gene %in% btneuronstopxcrispr)

# View the filtered table
head(filtered_topmarkers)
dim(filtered_topmarkers) #should match dim of btneuronsxcrispr

DoHeatmap(btneurons, features = filtered_topmarkers$gene, slot = "data")

filtered_unique.hvg <- filtered_topmarkers %>%
  rowwise() %>%                                    # Operate on rows
  mutate(num_clusters = sum(!is.na(c_across(-gene)))) %>%  # Count non-NA clusters (exclude the gene column)
  filter(num_clusters < 3) %>%               # Keep genes expressed in specified # of clusters
  select(-num_clusters)                            # Drop the helper column

dim(filtered_unique.hvg)
DoHeatmap(btneurons, features = filtered_unique.hvg$gene, slot = "data")

#make a list of only relevant neurological disease genes in Shawn's CRISPR knockout list
#filt.all_neuro_crispr <- intersect(filt.crisprscreen_genes, all_neuro_dis) 

#filter out essential genes
#filt.all_neuro_crispr <- setdiff(filt.all_neuro_crispr, hesc_ess)
#filt.all_neuro_crispr <- setdiff(filt.all_neuro_crispr, hap_ess)
#length(filt.all_neuro_crispr) #goal is ~250 genes

#save the list!
#write.csv(filt.all_neuro_crispr, file = "final_gene_list_neuron.csv", row.names = FALSE)



DoHeatmap(btneurons, features=btneuronsxcrispr, slot="data")

#load in CRISPR hits
crisprhits<- read.csv("d21_to_iPSC_all_fdr0.05_product0.6000000000000001_hits.csv")
ggplot(crisprhits, aes(x=epsilon, y=-log10(pvalue)))+
  geom_point(shape=21)+
  theme_light()

crisprscreen$hit <- ifelse(crisprscreen$gene == "NTC", "NTC",
                           ifelse(crisprscreen$gene %in% crisprhits$gene & crisprscreen$epsilon < 0, "Negative hit",
                           ifelse(crisprscreen$gene %in% crisprhits$gene & crisprscreen$epsilon > 0, "Positive hit", "other")))
table(crisprscreen$hit)

ggplot(crisprscreen, aes(x=epsilon, y=-log10(pvalue), color=hit))+
  geom_point(size=1, shape=21) +
  theme_light() +
  scale_color_manual(values=c("darkblue", "grey", "salmon", "darkred")) +
  geom_vline(xintercept=-2, linetype="dashed", color="darkgreen") 





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

library(readr)
gencc_allneurological <- read_csv("gencc_allneurological.csv")
View(gencc_allneurological)
