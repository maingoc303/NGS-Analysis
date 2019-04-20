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


####################################################################

#column sums of the count by sample
sapply(gene.exp, sum)


# Normalization count data
L = dim(gene.exp)[2]

# write a function for Median Log Deviation normalization
sum_log <- function(gene_count_sample){
  S=0
  for (i in 1:length(gene_count_sample))
    S=S+sum(log(gene_count_sample[i]+1))
  return(S)
}

# applying Median Log Deviation transformation

# create new dataframe for normalized counts
norm_count <- data.frame(matrix(0, nrow = dim(gene.exp)[1], ncol = 41))
names(norm_count)=names(gene.exp)
row.names(norm_count)=row.names(gene.exp)

for (i in 1:dim(norm_count)[2])
  for (j in 1:dim(norm_count)[1])
    norm_count[j,i]=log(gene.exp[j,i]+1)-1/L*sum_log(gene.exp[,i])

# Scaling factor for each sample





groups <- c(rep("P", 17), rep("N", 24))
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
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#The argument normalized equals true, divides each column by its size factor.
logcounts <- log2(counts(dds, normalized=TRUE) + 1)
pc <- prcomp(t(logcounts))

######################## EXPLORATORY ANALYSIS ####################
library(rafalib)
#PC 1 VS PC 2
plot(pc$x[,1], pc$x[,2], 
     col=colData(dds)$groups # color is the protocol
     )# point shape by time

#hierarchical clustering
plot(hclust(dist(t(logcounts))), labels=colData(dds)$groups)

#log count of the first two samples
plot(logcounts[,1], logcounts[,2], cex=.1)

# normalization to stabilize variance (regularized logarithm)
rld <- rlog( dds ) # object of class GenomicRanges 
pc2 <- prcomp( t( assay(rld) ) )
plot(pc2$x[,1], pc2$x[,2],  col=colData(rld)$groups)
plot(hclust(dist(t(assay(rld)))), labels=colData(rld)$groups)
plot(assay(rld)[,1], assay(rld)[,2], cex=.1)

colData(dds)$groups

colData(dds)$groups <- relevel(colData(dds)$groups, "P")
levels(colData(dds)$groups)

#Change the design
design(dds) <- ~ groups

# runs the model, extract a results table for all genes
dds <- DESeq( dds )
res <- results( dds )
head(res)

#For the protocol variable I want the control over L5 SNL
head(results(dds, contrast=c("groups","P","N")))

plotMA(dds)
plotMA(res)

resBigFC <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs")
plotMA(resBigFC)#, ylim=c(-5,5))
abline(h=c(-1,1),lwd=5)

#Top genes : sort by pvalue
resSort <- res[order(res$pvalue),]
head(resSort)

#Count for the first gene: the unnormalized count
k <- counts(dds)[rownames(resSort)[1],]
#Make a stripchart by combining time and protocol
cond <- with(colData(se), factor(groups))
stripchart(log2(k + 1) ~ cond, method="jitter", vertical=TRUE, las=2)

library(org.Rn.eg.db)
BiocInstaller::biocLite("org.Rn.eg.db")

library(org.Rn.eg.db)
#The keys we can query on Ensembl
keytypes(org.Rn.eg.db)

library("AnnotationDbi")
library("org.Hs.eg.db")
res<-results(dds,alpha=.05, contrast=c("groups", "N", "P"))
res$symbol <- mapIds(org.Hs.eg.db,keys=row.names(res),
                     column="SYMBOL", keytype="ENSEMBL", 
                     multiVals="first")
res$gene <- mapIds(org.Hs.eg.db,keys=row.names(res),
                     column="GENENAME", keytype="ENSEMBL", 
                     multiVals="first")