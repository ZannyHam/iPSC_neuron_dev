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

### read in and plot flo+bt data

#Load Flo's visualization Seurat object - Barbara's group's data + 5K of Flo's
flobtneurons <- readRDS("processed_seurat_obj_with_gRNA_data.rds")
flobtneurons

#Plotting UMAP with Seurat
DimPlot(flobtneurons, reduction = "umap", group.by = c("timepoint")) 
DimPlot(flobtneurons, reduction = "umap", group.by = c("batch"))

### read in and fine-tune Shawn's CRISPR scores

#Read Shawn's results file
shawn <- read.table("neuro_CRISPR_scores.csv", sep=",", header=T)
table(is.na(shawn$X))
shawn$X <- NULL #drop column, there is nothing there

head(shawn)
#which(grepl("ZBTB18", Features(flobtneurons)))
counts <- LayerData(flobtneurons)
#counts[which(grepl("ZBTB18", Features(flobtneurons))), 1:10]      

#Plotting gene expression on UMAP in Seurat
FeaturePlot(object = flobtneurons, features = head(shawn$gene, n=6))

### do the same in Monocle3 for my analysis ###

#Read in our analysis of Flo's data
filteredflozh <- readRDS("20241021flo_filtered_umi1k><25k_gene750><7k_mito12.RDS")
filteredflozh

#rowData(filteredflozh)[(grepl("ZBTB18", rowData(filteredflozh)$V2)), ]

names(rowData(filteredflozh))[1] <- "gene_short_name"

#Plot gene expression on UMAP in Monocle3
#plot_cells(filteredflozh, genes = head(shawn$gene, n=6))
#plot_cells(filteredflozh, genes = head(shawn$gene, n=6), scale_to_range=F)

head(colData(filteredflozh))

#How to filter Shawn's file
head(shawn)
table(shawn$index)

filtershawn <- filter(shawn, gene != "NTC") %>%
  filter(epsilon >= -2 & epsilon <= 2)

#orig <- data.frame(fread("iPSC_d0_d21_pDNA_merged.csv", sep=",", header=T))
#head(orig)

#### Make a doc with the expression matrix by cluster for my filtered flo neurons

# Access the expression matrix
flozh_expression_matrix <- exprs(filteredflozh)

# Extract cluster assignments
#flozh_clusters <- filteredflozh@clusters$UMAP$clusters  # Adjust based on your clustering method

flozh_clusters <- clusters(filteredflozh) #same thing as line 74 - quicker

# Convert expression matrix to a data frame
flozh_expression_data <- as.data.frame(as.matrix(flozh_expression_matrix))
flozh_expression_data$cluster <- flozh_clusters[rownames(flozh_expression_data)]

# Reshape data for easier manipulation
flozh_expression_data_long <- flozh_expression_data %>%
  pivot_longer(cols = -cluster, names_to = "gene", values_to = "expression")

# Calculate mean expression per gene within each cluster
flozh_expression_summary <- flozh_expression_data_long %>%
  group_by(cluster, gene) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop")

write.csv(flozh_expression_summary, file = "flozh_gene_expression_by_cluster_monocle.csv", row.names = FALSE)

### find marker genes in flobtneurons seurat S4 object and output as a .csv

#remove flo's 5000 cells as they are not assigned to clusters
missing_clusters <- is.na(flobtneurons@meta.data$seurat_clusters)
sum(missing_clusters)  # Count cells with missing cluster information
btneurons <- subset(flobtneurons, cells = colnames(flobtneurons)[!missing_clusters])
average_expression <- AggregateExpression(btneurons, group.by = "seurat_clusters", assays = "RNA")

flobt_allmarkers <- FindAllMarkers(btneurons,
              min.pct=0.1,
              logfc.threshold=0.25, # could go up to 0.5 if you want a more stringent differentiation threshold.
              test.use="MAST", #default is wilcox and is also suitable, can try both.
              latent.vars=c("batch") #use for MAST when there are technical or biological replicates that could confound results.
              )

#to do: find out what column names there are in the seurat object and make sure all these tests are
#being done by cluster. 
#merge across this data with shawn's data and filter for thresholds
#look at other papers Shawn sent to filter out potential genes. 

#### Make a doc with the expression matrix by cluster for flobtneurons

write.csv(expression_summary, file = "bt_gene_expression_by_cluster_monocle.csv", row.names = FALSE)

