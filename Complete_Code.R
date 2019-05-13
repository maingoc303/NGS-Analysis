# Library/ R Packages
library(dplyr)
library(Biobase)
library(GenomicRanges)
library(SGSeq)
library(DESeq2)
####################################################################

# Import data sets

annot <- read.csv("annotation.csv",header = T, sep=",")
names(annot)<-c("Patient","Group")

gene.exp <- read.csv("data.csv",header = T,sep=",",row.names=1)
names(gene.exp)
dim(gene.exp)
head(annot)

####################################################################
#                   VISUALIZATION                                  #
####################################################################

png("p1.png")
plot(sapply(gene.exp, sum)/1000000,col=(annot$Group+2),
     xlab="subject",ylab="Sum of count (mil)",
     main="Sum of all counts by subject")
legend(34,8.8,c("Placebo","Novel"),col=c(3,2),pch=c(20,20))
abline(h=mean(sapply(gene.exp, sum)/1000000),lwd=2,col=4)
dev.off()

groups <- c(rep("Placebo", 17), rep("Novel", 24))
se <- DESeqDataSetFromMatrix(countData = gene.exp, 
                             colData = as.data.frame(groups) ,
                             design = ~ groups)
dds <- DESeqDataSetFromMatrix(countData = gene.exp, 
                              colData = as.data.frame(groups) ,
                              design = ~ groups)

#Estimate size factors
dds <- estimateSizeFactors( dds )
sizeFactors(dds)
colSums(counts(dds))

#Plot column sums according to size factor
png("p2.png")
plot(sizeFactors(dds), colSums(counts(dds))/1000000, col=4,
     xlab="scaling factor",ylab="Sum of reads count")
abline(lm(colSums(counts(dds))/1000000 ~ sizeFactors(dds) + 0))
title("Scaling factor for count data in each subject")
dev.off()

logcounts <- log2(counts(dds, normalized=TRUE)+1)
norm_count = counts(dds, normalized=TRUE)

# visualize after normalization
png("p3.png")
plot(colMeans(norm_count),col=annot$Group+2, xlab="subject",ylab="Mean of count",
     main = "Mean of reads count after normalisation")
dev.off()
png("p4.png")
plot(rowMeans(norm_count),col=annot$Group+2, xlab="gene",ylab="Mean of count",
     main = "Mean of reads count after normalisation")
legend(45000,95000,c("Placebo","Novel"),cex=0.8,col=c(3,2),pch=c(20,20))
dev.off()

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

# means versus variance of RNA-seq expression
png("var.png")
plot(rowMeans(norm_count),rowVars(norm_count),col=3,pch=18)
dev.off()

V####################################################################
#                   TESING DE GENES                                 #
####################################################################
# VOOM transformation to apply linear regression and limma
library(edgeR)

png("voom1.png")
voom1 <- voom(norm_count, mod, plot = T)
dev.off()

d0 <- DGEList(norm_count)
# Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left 

snames <- colnames(norm_count) # Sample names

mm <- model.matrix(~ groups)
png("voom.png")
par(mfrow=c(1,2))
tmp <- voom(d0, mm, plot = T)#, main="before fitering low expressed gene")
y <- voom(d, mm, plot = T)
dev.off()

# fitting linear model in limma
fit <- lmFit(y, mm)
head(coef(fit))
fit.e<-eBayes(fit)
top.table <- topTable(fit.e, sort.by = "p", n = Inf)
xtable(head(top.table, 10))
length(which(top.table$adj.P.Val < 0.05))

gene_names = select(org.Hs.eg.db, rownames(top.table)[1:7], 
                    c("SYMBOL","GENENAME","ENTREZID"),"ENSEMBL")

gene <- mapIds(org.Hs.eg.db,keys=row.names(top.table)[1:7],
               column="GENENAME", keytype="ENSEMBL", 
               multiVals="first")
#