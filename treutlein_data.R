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

#add all the data into the cds object
btcds <- load_mm_data(mat_path = "treutlein_data/matrix.mtx.gz", 
                    feature_anno_path = "treutlein_data/features.tsv.gz", 
                    cell_anno_path = "treutlein_data/barcodes.tsv.gz")

btcds$sanitycheck <- colSums((counts(btcds)))
head(colData(btcds))

#find genes in the cds object
btcds <- detect_genes(btcds)

btcds$lane <- sapply(strsplit(as.character(btcds$cell), "-"), `[`, 2) #splits into the lanes by separating at the hyphen
table(btcds$lane) #counts the number of cells in each lane

btcds$lane <- as.factor(sapply(strsplit(as.character(btcds$cell), "-"), `[`, 2))
table(btcds$lane)
levels(btcds$lane)
btcds$lane <- gsub(" ", "", btcds$lane)
btcds$lane2 <- as.factor(btcds$lane, levels=c("2", "3", "4", "6", "7", "8", "9", "10", "11", "12", "13"))
table(btcds$lane)

####### let's plot UMIs per sample
#a density plot of UMIs per sample
ggplot(data.frame(colData(btcds)), aes(y=log10(n.umi), color=lane)) +
  geom_density() +
  coord_flip() +
  theme_light()

#summary stat function
#displays the number of cells as a label at the max value of y
bt_n_fun <- function(y){return(data.frame(y=max(y), label = paste0(length(y))))}

#violin plot combined with box plot visualization
bt_umi <- ggplot(colData(btcds), aes(x=lane, y=n.umi)) +
  #facet_wrap(~plate, ncol=1, drop=FALSE, scales="free") +
  geom_violin(aes(fill=lane)) +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 16),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size = 16)) +
  scale_y_log10() +
#this is where I could change the intercept parameters for UMI counts
  geom_hline(yintercept=35000, linetype="dotted", color="red") +
  geom_hline(yintercept=1500, linetype="dotted", color="red") +
  geom_hline(yintercept=3000, linetype="dotted", color="red") +
  xlab("Lane") +
  ylab("UMIs") +
  stat_summary(fun.data = bt_n_fun, geom = "text", hjust = 0.5, vjust = -0.3) +
  theme(legend.position="none")
bt_umi 

#summarize the UMI distribution per lane (like means and max etc)
data.frame(colData(btcds)) %>% 
  group_by(lane) %>%
  summarise(min=min(n.umi),
            q1=quantile(n.umi, probs=c(0.25)),
            med=median(n.umi),
            mean=mean(n.umi),
            q3=quantile(n.umi, probs=c(0.75)),
            max=max(n.umi)) %>%
  data.frame()

#detect outliers based on UMI count (MAD = median absolute deviations)
# What thresholds should we use? 
# The number of UMIs are not normally distributed so we use the log. Perhaps
# we exclude anything > 3 MADS and lower < 1 MAD?
bt_thresh_umi <- isOutlier(btcds$n.umi, 
                      log=TRUE,
                      nmads=3, 
                      type=c("both"), #the mad will be applied to the high and low ends. Could do "lower" or "higher" instead
                      batch=paste(btcds$lane)) #done for each lane because of potential batch differences
attr(bt_thresh_umi, "thresholds")
#output has low threshold on top and high threshold on the bottom
#3 MADS ~1K, ~25K #probably what I would use
#5 MADS ~50-70K


#### plot Number of genes per lane in a violin and box plot
bt_genes <- ggplot(data.frame(colData(btcds)), aes(x=lane, y=num_genes_expressed)) +
  geom_violin(aes(fill=lane)) +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light() +
  geom_hline(yintercept = 6000, linetype="dotted", color="red") + #change to thresholds
  geom_hline(yintercept = 750, linetype="dotted", color="red") + #change to thresholds
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 16),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size = 16)) +
  scale_y_log10() +
  xlab("lane") +
  ylab("Genes") +
  stat_summary(fun.data = bt_n_fun, geom = "text", hjust = 0.5, vjust = -0.3)
bt_genes

data.frame(colData(btcds)) %>% 
  group_by(lane) %>%
  summarise(min=min(num_genes_expressed),
            q1=quantile(num_genes_expressed, probs=c(0.25)),
            med=median(num_genes_expressed),
            mean=mean(num_genes_expressed),
            q3=quantile(num_genes_expressed, probs=c(0.75)),
            max=max(num_genes_expressed)) %>%
  data.frame()


bt_thresh_genes <- isOutlier(btcds$num_genes_expressed, 
                        log=TRUE,
                        nmads=3, 
                        type=c("both"), 
                        batch=paste(btcds$lane))
attr(bt_thresh_genes, "thresholds")
#3 MADS 750, 7000

####Calculate mitochondrial and ribosomal content

#calculate mito content (MT = Mito genes) and % of mito reads
fData(btcds)$MT <- grepl("^MT-", rownames(rowData(btcds))) #need to modify for Treutlein data

pData(btcds)$MT_reads <- Matrix::colSums(exprs(btcds)[fData(btcds)$MT,])
pData(btcds)$perc_mitochondrial_umis <- btcds$MT_reads/btcds$n.umi*100

#calculate ribo content (ribo = ribo genes) and % of ribo reads
fData(btcds)$ribo <- grepl("^RPL|^RPS", rownames(rowData(btcds)))
pData(btcds)$ribo_reads <- Matrix::colSums(exprs(btcds)[fData(btcds)$ribo,])
pData(btcds)$perc_ribo_umis <- btcds$ribo_reads/btcds$n.umi*100

#### Plot mitochondrial reads as a violin and box plot
bt_mito <- ggplot(data.frame(colData(btcds)), aes(x=lane, y=perc_mitochondrial_umis)) +
  #facet_wrap(~plate, ncol=1, drop=FALSE, scales="free") +
  geom_violin(aes(fill=lane)) +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light() +
  geom_hline(yintercept = 5, linetype="dotted", ) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 16),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size = 16)) +
  scale_y_continuous(limits=c(0,100)) +
  geom_hline(yintercept=15, linetype="dotted", color="red")+
  xlab("lane") +
  ylab("% Mitochondrial") +
  stat_summary(fun.data = bt_n_fun, geom = "text", hjust = 0.5, vjust = -0.3)
bt_mito

###Scatter plot showing relationship between UMI counts and number of genes expressed
#Colored by mito content
data.frame(colData(btcds)) %>% 
  group_by(lane) %>%
  summarise(min=min(perc_mitochondrial_umis),
            q1=quantile(perc_mitochondrial_umis, probs=c(0.25)),
            med=median(perc_mitochondrial_umis),
            mean=mean(perc_mitochondrial_umis),
            q3=quantile(perc_mitochondrial_umis, probs=c(0.75)),
            max=max(perc_mitochondrial_umis)) %>%
  data.frame()

bt_thresh_mito <- isOutlier(btcds$perc_mitochondrial_umis, 
                       nmads=1, type=c("higher"), 
                       batch=paste(btcds$lane))
attr(bt_thresh_mito, "thresholds")
#What threshold should we use?

# All the QC together for scatter plot
ggplot(data.frame(colData(btcds)), aes(x=n.umi, y=num_genes_expressed)) +
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


#filter out cells based on UMI thresholds, gene expression, mito content
#btfiltered <- btcds[, btcds$n.umi >= 1000 &
#btcds$num_genes_expressed >= 700 &
#btcds$perc_mitochondrial_umis <= 5] #Zanny - change thresholds, add mitochondrial

btfiltered <- btcds #zero filtering
btfiltered <- btfiltered[rowData(btfiltered)$num_cells_expressed > 0,]#keep only expressed genes

head(rowData(btfiltered))


set.seed(1000) #starts with same cell to begin the calculations 
btfiltered <- preprocess_cds(btfiltered)
plot_pc_variance_explained(btfiltered) #makes a knee plot

# reduce dimensions
set.seed(1000)
btfiltered <- reduce_dimension(btfiltered)

# add UMAP coordinates to the colData for easy plotting outside of monocle3
btfiltered$UMAP1 <- reducedDim(btfiltered, "UMAP")[,1]
btfiltered$UMAP2 <- reducedDim(btfiltered, "UMAP")[,2]

#extract the meta info into a data frame
btmeta <- data.frame(colData(btfiltered))

class(btmeta)
str(btmeta)

# group cells into clusters
# you can spend a lot of time tweaking clusering (and should!)
k=ceiling(sqrt(dim(btfiltered)[2])*0.25)
k
btfiltered <- cluster_cells(btfiltered, k=50)
table(clusters(btfiltered))
table(partitions(btfiltered))

plot_cells(btfiltered, color_cells_by = "partition")
plot_cells(btfiltered)

plot_cells(btfiltered, color_cells_by = "sample")

# add cluster information to colData for each cell 
colData(btfiltered)$clusters = clusters(btfiltered)

# check number of clusters and partitions
head(clusters(btfiltered, reduction_method = "UMAP"))
head(partitions(btfiltered, reduction_method = "UMAP"))


# plot!
ggplot(colData(btfiltered), aes(x=UMAP1, y=UMAP2, color=clusters)) +
  facet_wrap(~lane, ncol=4) +
  geom_point(alpha=0.5, size=0.1) +
  theme_light() +
  theme(legend.position="bottom", aspect.ratio = 1) +
  guides(color = guide_legend(override.aes = list(size=6, alpha=1))) 

ggplot(colData(btfiltered), aes(x=UMAP1, y=UMAP2, color=perc_mitochondrial_umis)) +
  #facet_wrap(~lane, ncol=4) +
  geom_point(alpha=0.5, size=0.1) +
  theme_light() +
  scale_color_viridis() + 
  theme(legend.position="bottom", aspect.ratio = 1) +
  guides(color = guide_legend(override.aes = list(size=6, alpha=1))) 
