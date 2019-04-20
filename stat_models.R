# poisson model for count data
plotDispEsts( dds)#, ylim = c(1e-6, 1e1) )

hist( res$pvalue, breaks=20, col="grey" )

# create bins using the quantile function
qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
# "cut" the genes into the bins
bins <- cut( res$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of p values less than .01 for each bin
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")

attr(res,"filterThreshold")
plot(attr(res,"filterNumRej"),type="b",
     xlab="quantiles of 'baseMean'",
     ylab="number of rejections")

sampleDists <- dist( t( assay(rld) ) )
as.matrix( sampleDists )[ 1:3, 1:3]

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$groups, 
                                     rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL   

library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

colours <- c(rgb(1:3/4,0,0),rgb(0,1:3/4,0),rgb(0,0,1:3/4),rgb(1:3/4,0,1:3/4))
plotPCA( rld, intgroup = c("sizeFactor","groups"))#, col=colours )

# gene clustering
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 15 )

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
             colData(rld)$groups ] )


library(limma)
library(edge)
library(tidyr)
mod =model.matrix(~annot$Group)
fit = lm.fit(mod,t(norm_count))

names(fit)
tidy(fit$coefficients[,1])

fit_limma1 = lmFit((norm_count),mod)
names(fit_limma)
fit_limma$coefficients[1,]

summary(fit_limma)

fit1 <- eBayes(fit_limma)
head(fit1$t)
head(fit1$p.value)

topTable <- topTable(fit1,number=nrow(norm_count),coef=2)
head(topTable)
dim(topTable)
dim(topTable[topTable[,4]<=0.05,])

length(topTable[,4])
#png("limma1.png")
#plot(t.t,-log(p.val),col=4,cex=.8,main="T-test versus LIMMA")
#points(topTable[,3],-log(topTable[,4]),col=2,cex=.7)
#legend(-2,5,legend = c("t-test","limma"),col=c(4,2),pch=c(1,1))
#dev.off()


fit_limma1 = lmFit((logcounts),mod)
names(fit_limma1)
fit_limma1$coefficients[1,]

summary(fit_limma1)

fit2 <- eBayes(fit_limma1)
head(fit2$t)
head(fit2$p.value)

topTable <- topTable(fit2,number=nrow(logcounts),coef=2)
head(topTable)
dim(topTable)
dim(topTable[topTable[,4]<=0.05,])


#### --------------------------------------------------------
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
snames
#multidimensional scaling
plotMDS(d, col = as.numeric(annot$Group)+2)

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
##############################################################
# negative binomial regression 
library(MASS)
m1 <- glm.nb(norm_count ~ mod, data = as.data.frame(norm_count))
m1 <- glm(logcounts[1,] ~ as.factor(annot$Group), data = as.data.frame(logcounts))

pval=c(1:dim(norm_count)[1])
for (i in 1:dim(norm_count)[1]){
  m1 <- glm(norm_count[i,]~as.factor(annot$Group),family = "poisson")
  pval[i]=coef(summary(m1))[,4][2]
}

holm<-mt.rawp2adjp(pval, proc=c("Holm"))
bonf<-mt.rawp2adjp(pval, proc=c("Bonferroni"))
bh<-mt.rawp2adjp(pval, proc=c("BH"))
by<-mt.rawp2adjp(pval, proc=c("BY"))
allp<-cbind(pval, bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2])
par(mfrow=c(1,1))
mt.plot(allp,plottype="pvsr", proc=c("rawp","Bonferroni","Holm","BH","BY"),lty=c(1,2,3,4,5),col=1:5,lwd=2)#,leg=c(2100,0.375))







