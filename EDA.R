# visualized raw reads data

plot(colMeans(gene.exp),col=annot$Group+2)
plot(rowMeans(gene.exp),col=annot$Group+2)

# after normalization the reads count
 
plot(colMeans(norm_count),col=annot$Group+2)
plot(rowMeans(norm_count),col=annot$Group+2)

# principal component analysis
# package:
library(factoextra)

pc1 = prcomp((norm_count))
#plot(pc1$rotation[,1],pc1$rotation[,2],col=annot$Group+3,pch=18,cex=1.2,
#     xlab="PC1",ylab="PC2",main="Biplot for PCA")

png("p5.png")
fviz_eig(pc1)
dev.off()
png("p6.png")
groups <- as.factor(annot$Group)
fviz_pca_ind(pc1,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             #addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
dev.off()
biplot(pc1)

dim(t(norm_count))

logcounts <- log2(norm_count+1)
pc2 <- prcomp(t(logcounts))
fviz_eig(pc2)
fviz_pca_ind(pc2,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             #addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

svd1 = svd(norm_count)
plot(svd1$v[,1],svd1$v[,2],col=annot$Group+2)

# clustering
par(mfrow=c(1,1))
plot(hclust(dist(t(logcounts))), labels=colData(dds)$groups,cex=.6)#,col=colData(dds)$groups)
plot(hclust(dist(t(norm_count))), labels=colData(dds)$groups,cex=.6)

plot(rowMeans(norm_count),rowVars(norm_count))#,labels=rownames(norm_count))

var.exp = c(1:dim(norm_count)[1])

# Clustering
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)

res <- results( dds )

# heat map
library(ComplexHeatmap)
BiocManager::install("ComplexHeatmap")
library(circlize)  

# get the top genes
sigGenes <- as.data.frame(resSort) %>% 
  rownames_to_column("GeneID") %>% 
  top_n(12, wt=-padj) %>% 
  pull("GeneID")
# filter the data for the top 200 by padj in the LRT test
plotDat <- vst(dds)[sigGenes,] %>% 
  assay()
z.mat <- t(scale(t(norm_count), center=TRUE, scale=TRUE))

# colour palette
myPalette <- c("red3", "ivory", "blue3")
myRamp = colorRamp2(c(-2, 0, 2), myPalette)

# Heatmap
Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_name = FALSE,
        cluster_columns = FALSE)

png("var.png")
plot(rowMeans(norm_count),rowVars(norm_count),col=3,pch=18)
dev.off()
length(rowMeans(norm_count))
sum(rowMeans(norm_count)>2)  
  