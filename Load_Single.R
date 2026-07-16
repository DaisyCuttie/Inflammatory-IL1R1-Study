library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)
library(harmony)
# Load the HTO data
D4neg = read.csv("D4negative_barcode/D4neg_hashtag_results.csv", row.names = 1)
typeof(D4neg)
rownames(D4neg)
rownames(D4neg) = paste(rownames(D4neg), "-1", sep = "")
D4neg$Idents.sample.hashtag. = gsub("TGATGGCCTATTGGG", "media_neg", D4neg$Idents.sample.hashtag)
D4neg$Idents.sample.hashtag. = gsub("TTCCGCCTCTCTTTG", "IL1B_neg", D4neg$Idents.sample.hashtag)
D4neg$Idents.sample.hashtag. = gsub("GGTTGCCAGATGTCA", "CD3_neg", D4neg$Idents.sample.hashtag)
D4neg$Idents.sample.hashtag.1 = c("Donor4")
rownames(D4neg)

D4pos = read.csv("D4positive_barcode/D4pos_hashtag_results.csv", row.names = 1)
typeof(D4pos)
rownames(D4pos)
rownames(D4pos) = paste(rownames(D4pos), "-2", sep = "")
D4pos$Idents.sample.hashtag. = gsub("TGATGGCCTATTGGG", "media_pos", D4pos$Idents.sample.hashtag)
D4pos$Idents.sample.hashtag. = gsub("TTCCGCCTCTCTTTG", "IL1B_pos", D4pos$Idents.sample.hashtag)
D4pos$Idents.sample.hashtag. = gsub("GGTTGCCAGATGTCA", "CD3_pos", D4pos$Idents.sample.hashtag)
D4pos$Idents.sample.hashtag.1 = c("Donor4")
rownames(D4pos)

D5neg = read.csv("D5negative_barcode/D5neg_hashtag_results.csv", row.names = 1)
typeof(D5neg)
rownames(D5neg)
rownames(D5neg) = paste(rownames(D5neg), "-3", sep = "")
D5neg$Idents.sample.hashtag. = gsub("TGATGGCCTATTGGG", "media_neg", D5neg$Idents.sample.hashtag)
D5neg$Idents.sample.hashtag. = gsub("TTCCGCCTCTCTTTG", "IL1B_neg", D5neg$Idents.sample.hashtag)
D5neg$Idents.sample.hashtag. = gsub("GGTTGCCAGATGTCA", "CD3_neg", D5neg$Idents.sample.hashtag)
D5neg$Idents.sample.hashtag.1 = c("Donor5")
rownames(D5neg)

D5pos = read.csv("D5positive_barcode/D5pos_hashtag_results.csv", row.names = 1)
typeof(D5pos)
rownames(D5pos)
rownames(D5pos) = paste(rownames(D5pos), "-4", sep = "")
D5pos$Idents.sample.hashtag. = gsub("TGATGGCCTATTGGG", "media_pos", D5pos$Idents.sample.hashtag)
D5pos$Idents.sample.hashtag. = gsub("TTCCGCCTCTCTTTG", "IL1B_pos", D5pos$Idents.sample.hashtag)
D5pos$Idents.sample.hashtag. = gsub("GGTTGCCAGATGTCA", "CD3_pos", D5pos$Idents.sample.hashtag)
D5pos$Idents.sample.hashtag.1 = c("Donor5")
rownames(D5pos)

D3neg = read.csv("../Prisila_raw/D3negative_barcode/D3neg_hashtag_results.csv", row.names = 1)
typeof(D3neg)
rownames(D3neg)
rownames(D3neg) = paste(rownames(D3neg), "-5", sep = "")
D3neg$Idents.sample.hashtag. = gsub("-TGATGGCCTATTGGG", "_neg", D3neg$Idents.sample.hashtag)
D3neg$Idents.sample.hashtag. = gsub("-TTCCGCCTCTCTTTG", "_neg", D3neg$Idents.sample.hashtag)
D3neg$Idents.sample.hashtag. = gsub("-GGTTGCCAGATGTCA", "_neg", D3neg$Idents.sample.hashtag)
rownames(D3neg)
D3neg$Idents.sample.hashtag.1 = c("Donor3")

D3pos = read.csv("../Prisila_raw/D3positive_barcode/D3pos_hashtag_results.csv", row.names = 1)
typeof(D3pos)
rownames(D3pos)
rownames(D3pos) = paste(rownames(D3pos), "-6", sep = "")
D3pos$Idents.sample.hashtag. = gsub("-TGATGGCCTATTGGG", "_pos", D3pos$Idents.sample.hashtag)
D3pos$Idents.sample.hashtag. = gsub("-TTCCGCCTCTCTTTG", "_pos", D3pos$Idents.sample.hashtag)
D3pos$Idents.sample.hashtag. = gsub("-GGTTGCCAGATGTCA", "_pos", D3pos$Idents.sample.hashtag)
D3pos$Idents.sample.hashtag.1 = c("Donor3")
rownames(D3pos)

combined_hashtag = rbind(D4neg, D4pos, D5neg, D5pos, D3neg, D3pos)

sample.data <- Read10X(data.dir = "RNA_aggr/outs/count/filtered_feature_bc_matrix/")
colnames(sample.data)
tail(colnames(sample.data))
sample <- CreateSeuratObject(counts = sample.data, project = "", min.cells = 3, min.features = 200)
common_bars = intersect(colnames(sample), rownames(combined_hashtag))
head(common_bars)
tail(common_bars)
sample = subset(sample, cells = common_bars)
combined_hashtag = combined_hashtag[common_bars,]

sample$HTO = as.character(combined_hashtag$Idents.sample.hashtag.)
sample$Donor = as.character(combined_hashtag$Idents.sample.hashtag.1)

#make sure that we have all the samples
unique(sample$HTO)
unique(sample$Donor)
rownames(sample)

#start data processing
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# check which mitochondrial genes does it have
grep ("^MT-", rownames(sample[["RNA"]]),value = T)
plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Trying without harmony
# Normalize the data
sample = SCTransform(sample, vars.to.regress = "percent.mt", verbose = TRUE)
# Perform linear dimensional reduction
sample = RunPCA(sample, verbose = TRUE)
# determine the dimensionality of the dataset
ElbowPlot(sample)
sample = RunUMAP(sample, reduction = "pca", dims = 1:12)
sample = FindNeighbors(sample, reduction = "pca", dims = 1:12)
sample = FindClusters(sample, resolution = 0.6)
DimPlot(sample, reduction = "umap", label=TRUE, group.by = "orig.ident")
DimPlot(sample, reduction = "umap", label=TRUE)
DimPlot(sample, reduction = "umap", group.by="HTO")
DimPlot(sample, reduction = "umap", group.by="Donor")

Idents(sample) = sample$seurat_clusters

VlnPlot(object=sample, features=c("nFeature_RNA","nCount_RNA"))

# Trying harmony
sample = NormalizeData(sample, verbose = FALSE)
sample = FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
sample = ScaleData(sample, verbose = FALSE)
sample = RunPCA(sample, verbose = TRUE)
# determine the dimensionality of the dataset
ElbowPlot(sample)
# run harmony to remove donor batch effect
sample = RunHarmony(sample, group.by.vars = "Donor", dims.use = 1:30)
sample <- sample %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()
# save the sample
saveRDS(sample, file = "single_cell_object_harmony.rds")
# Look at cluster IDs of the first 5 cells
#head(Idents(sample), 5)
#Idents(sample) = sample$seurat_clusters
DimPlot(sample, reduction = "umap", label=TRUE)
DimPlot(sample, reduction = "umap", group.by="HTO")
DimPlot(sample, reduction = "umap", group.by="Donor")

# Find all markers
sample.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sample.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
write.csv(sample.markers, file = "single_cell_markers.csv")

#sample.markers = read.csv("single_cell_markers.csv")
top3 <- sample.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
genes = unique(as.character(top5$gene))
DotPlot(sample, features = genes) + coord_flip()

DoHeatmap(sample, features = top3$gene) + NoLegend()

# plot stacked barplot for annotating the clusters
sample$orig.ident = sample$HTO
orig = sample$orig.ident
orig
u_orig = sort(unique(orig))
u_orig
cluster = as.character(Idents(sample))
uct = sort(unique(cluster))
z = matrix(nrow=length(uct), ncol=length(u_orig))
rownames(z) = uct
colnames(z) = u_orig
z
for (i in 1:nrow(z)){
  for (j in 1:ncol(z)){
    z[i,j] = sum(orig[cluster==uct[i]] == u_orig[j])
  }
}
# Make a plot about num cells for each cluster by orig
z = as.matrix(z)
z1 = data.frame(ncells=as.vector(z),
                type=rep(as.integer(rownames(z)), ncol(z)),
                orig=rep(colnames(z), each=nrow(z)))
ggplot(z1, aes(x=type, y=ncells, fill=orig, group=orig)) +
  geom_bar(stat="identity", width=0.7) +
  theme_classic() +
  labs(x="Cluster Number", y="Numer of Cells", fill=NULL)+
  scale_x_continuous(breaks=0:max(z1$type))

# Find average expression of all genes per cluster
cluster.averages <- AverageExpression(sample, assays = "RNA")
write.csv(cluster.averages, "ave_expr_by_cluster.csv")
