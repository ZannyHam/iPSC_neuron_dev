#Setup
TASKDIR = "~/Documents/UW/starita_lab/"
setwd("~/Documents/UW/starita_lab/")

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

flocds <- load_mm_data(mat_path = "neuron_data_filtered_features/matrix.mtx.gz", 
                    feature_anno_path = "neuron_data_filtered_features/features.tsv.gz", 
                    cell_anno_path = "neuron_data_filtered_features/barcodes.tsv.gz")

flocds$sanitycheck <- colSums((counts(flocds)))
head(colData(flocds))

flocds <- detect_genes(flocds)
flocds$lane <- as.factor(sapply(strsplit(as.character(flocds$cell), "-"), `[`, 2))
table(flocds$lane)
levels(flocds$lane)

####### let's plot UMIs per sample
ggplot(data.frame(colData(flocds)), aes(y=log10(n.umi), color=lane)) +
  geom_density() +
  coord_flip() +
  theme_light()

#summary stat function
flo_n_fun <- function(y){return(data.frame(y=max(y), label = paste0(length(y))))}
# my preferred visualization
flo_umi <- ggplot(colData(flocds), aes(x=lane, y=n.umi)) +
  #facet_wrap(~plate, ncol=1, drop=FALSE, scales="free") +
  geom_violin(aes(fill=lane)) +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 16),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size = 16)) +
  scale_y_log10() +
#this is where I could change the intercept parameters
  geom_hline(yintercept=25000, linetype="dotted", color="red") +
  geom_hline(yintercept=1000, linetype="dotted", color="red") +
  geom_hline(yintercept=3000, linetype="dotted", color="red") +
  xlab("Lane") +
  ylab("UMIs") +
  stat_summary(fun.data = flo_n_fun, geom = "text", hjust = 0.5, vjust = -0.3) +
  theme(legend.position="none")
flo_umi

data.frame(colData(flocds)) %>% 
  group_by(lane) %>%
  summarise(min=min(n.umi),
            q1=quantile(n.umi, probs=c(0.25)),
            med=median(n.umi),
            mean=mean(n.umi),
            q3=quantile(n.umi, probs=c(0.75)),
            max=max(n.umi)) %>%
  data.frame()


# What thresholds should we use? 
# The number of UMIs are not normally distributed so we use the log. Perhaps
# we exclude anything > 3 MADS and lower < 1 MAD?
flo_thresh_umi <- isOutlier(flocds$n.umi, 
                      log=TRUE,
                      nmads=3, 
                      type=c("both"), 
                      batch=paste(flocds$lane))
attr(flo_thresh_umi, "thresholds")
#3 MADS ~1K, ~25K #probably what I would use
#5 MADS ~50-70K


#### Number of genes
flo_genes <- ggplot(data.frame(colData(flocds)), aes(x=lane, y=num_genes_expressed)) +
  geom_violin(aes(fill=lane)) +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light() +
  geom_hline(yintercept = 7000, linetype="dotted", color="red") + #change to thresholds
  geom_hline(yintercept = 750, linetype="dotted", color="red") + #change to thresholds
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 16),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size = 16)) +
  scale_y_log10() +
  xlab("lane") +
  ylab("Genes") +
  stat_summary(fun.data = flo_n_fun, geom = "text", hjust = 0.5, vjust = -0.3)
flo_genes

data.frame(colData(flocds)) %>% 
  group_by(lane) %>%
  summarise(min=min(num_genes_expressed),
            q1=quantile(num_genes_expressed, probs=c(0.25)),
            med=median(num_genes_expressed),
            mean=mean(num_genes_expressed),
            q3=quantile(num_genes_expressed, probs=c(0.75)),
            max=max(num_genes_expressed)) %>%
  data.frame()


flo_thresh_genes <- isOutlier(flocds$num_genes_expressed, 
                        log=TRUE,
                        nmads=3, 
                        type=c("both"), 
                        batch=paste(flocds$lane))
attr(flo_thresh_genes, "thresholds")
#3 MADS 750, 7000

# Calculate mitochondrial and ribosomal content
#calculate mito content
fData(flocds)$MT <- grepl("^MT-", rowData(flocds)$V2)

pData(flocds)$MT_reads <- Matrix::colSums(exprs(flocds)[fData(flocds)$MT,])
pData(flocds)$perc_mitochondrial_umis <- flocds$MT_reads/flocds$n.umi*100

#calculate ribo content
fData(flocds)$ribo <- grepl("^RPL|^RPS", rowData(flocds)$V2)
pData(flocds)$ribo_reads <- Matrix::colSums(exprs(flocds)[fData(flocds)$ribo,])
pData(flocds)$perc_ribo_umis <- flocds$ribo_reads/flocds$n.umi*100

#### Mitochondrial reads
flo_mito <- ggplot(data.frame(colData(flocds)), aes(x=lane, y=perc_mitochondrial_umis)) +
  #facet_wrap(~plate, ncol=1, drop=FALSE, scales="free") +
  geom_violin(aes(fill=lane)) +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light() +
  geom_hline(yintercept = 5, linetype="dotted", ) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 16),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size = 16)) +
  scale_y_continuous(limits=c(0,100)) +
  geom_hline(yintercept=6, linetype="dotted", color="red")+
  xlab("lane") +
  ylab("% Mitochondrial") +
  stat_summary(fun.data = flo_n_fun, geom = "text", hjust = 0.5, vjust = -0.3)
flo_mito

data.frame(colData(flocds)) %>% 
  group_by(lane) %>%
  summarise(min=min(perc_mitochondrial_umis),
            q1=quantile(perc_mitochondrial_umis, probs=c(0.25)),
            med=median(perc_mitochondrial_umis),
            mean=mean(perc_mitochondrial_umis),
            q3=quantile(perc_mitochondrial_umis, probs=c(0.75)),
            max=max(perc_mitochondrial_umis)) %>%
  data.frame()

flo_thresh_mito <- isOutlier(flocds$perc_mitochondrial_umis, 
                       nmads=3, type=c("higher"), 
                       batch=paste(flocds$lane))
attr(flo_thresh_mito, "thresholds")
#Higher threshold ~12 seems good
#could eliminate mitochondria data and see how that affects Flo's UMAP

# All the QC together
ggplot(data.frame(colData(flocds)), aes(x=n.umi, y=num_genes_expressed)) +
  geom_point(aes(alpha=0.3, color=perc_mitochondrial_umis)) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 16),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size = 16),
        aspect.ratio = 1) +
  xlab("UMIs") +
  ylab("Number of Genes Captured") +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope=1, color="grey") #+
#scale_color_viridis(option="B") +
#geom_hline(yintercept = 1000, linetype="dotted", color="grey") + #change to thresholds
#geom_vline(xintercept = 300, linetype="dotted", color="grey") #+ #change to thresholds
#geom_vline(xintercept = 50000, linetype="dotted", color="grey")  #change to thresholds


#filter the data
flofiltered <- flocds[, flocds$n.umi >= 1000 & flocds$n.umi <= 25000 &
flocds$num_genes_expressed >= 750 & flocds$num_genes_expressed <= 7000
flocds$perc_mitochondrial_umis <= 12] #Zanny - change thresholds, add mitochondrial

#flofiltered <- flocds #zero filtering
#flofiltered <- flofiltered[rowData(flofiltered)$num_cells_expressed > 0,]#keep only expressed genes

set.seed(1000)
flofiltered <- preprocess_cds(flofiltered)
plot_pc_variance_explained(flofiltered)

# reduce dimensions
set.seed(1000)
flofiltered <- reduce_dimension(flofiltered)

# add UMAP coordinates to the colData for easy plotting outside of monocle3
flofiltered$UMAP1 <- reducedDim(flofiltered, "UMAP")[,1]
flofiltered$UMAP2 <- reducedDim(flofiltered, "UMAP")[,2]

#extract the meta info into a data frame
flo_meta <- data.frame(colData(flofiltered))

class(flo_meta)
str(flo_meta)

# group cells into clusters
# you can spend a lot of time tweaking clustering (and should!)
k=ceiling(sqrt(dim(flofiltered)[2])*0.25)
k
flofiltered <- cluster_cells(flofiltered, k=50)
table(clusters(flofiltered))
table(partitions(flofiltered))

plot_cells(flofiltered, color_cells_by = "partition")
plot_cells(flofiltered)

plot_cells(flofiltered, color_cells_by = "lane")
#plot_cells(flofiltered, color_cells_by = "sample") #didnt work, check columns?
# add cluster information to colData for each cell 
colData(flofiltered)$clusters = clusters(flofiltered)

# check number of clusters and partitions
head(clusters(flofiltered, reduction_method = "UMAP"))
head(partitions(flofiltered, reduction_method = "UMAP"))

# plot per lane
ggplot(colData(flofiltered), aes(x=UMAP1, y=UMAP2, color=clusters)) +
  facet_wrap(~lane, ncol=4) +
  geom_point(alpha=0.5, size=0.1) +
  theme_light() +
  theme(legend.position="bottom", aspect.ratio = 1) +
  guides(color = guide_legend(override.aes = list(size=6, alpha=1))) 

ggplot(colData(flofiltered), aes(x=UMAP1, y=UMAP2, color=perc_mitochondrial_umis)) +
  #facet_wrap(~lane, ncol=4) +
  geom_point(alpha=0.5, size=0.1) +
  theme_light() +
  scale_color_viridis() + 
  theme(legend.position="bottom", aspect.ratio = 1) +
  guides(color = guide_legend(override.aes = list(size=6, alpha=1))) 




####### Read in Shawn's file
ncs <- read.table("neuro_CRISPR_scores.csv", header=T, sep=",", stringsAsFactors=FALSE)
table(ncs$X)
table(is.na(ncs$X))
ncs$X <- NULL


ncs2 <- data.frame(fread("neuro_CRISPR_scores.csv"))
ncs2$V7 <- NULL

identical(ncs2, ncs)

ncs$treutlein <- ifelse(ncs$gene %in% rownames(rowData(flocds)), TRUE, FALSE)


rowData(flocds)$epsilon <- ncs[match(rownames(rowData(flocds)), ncs$gene), 'epsilon']

test <- data.frame(rowData(flocds)) %>% filter(!is.na(epsilon)) %>% arrange(desc(epsilon)) %>% head() 
rownames(test)


rowData(flofiltered)$gene_short_name <- rownames(rowData(flofiltered))
plot_cells(flofiltered, genes=rownames(test))
plot_cells(flofiltered, genes=rownames(test), scale_to_range = FALSE)
