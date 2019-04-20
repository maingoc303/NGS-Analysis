df <- read.csv("Data1.csv", header = T)

# --------- LIBRARY ---------- #
library(Biobase)
library(multtest, verbose = F)
library(ggplot2)

# data set overview
advanced<-c(1:99)
early<-c(1:99)
for (i in 1:99){
  advanced[i]<-mean(df[df$Stage=="advanced",i+2])
  early[i]<-mean(df[df$Stage=="early",i+2])
}

png("eda.png")
plot(advanced,pch=20,col=2,main="Mean of expression levels by metabolite",
     xlab="Metabolite",ylab="Expression levels")
points(early,pch=5,col=4)
legend(0,170,pch=c(20,5),col=c(2,4),legend=c("advanced","early"))
dev.off()

png("eda2.png")
par(mfrow=c(2,3))
boxplot(df$var91~df$Stage,col="light blue")
title("metabolite 91")
boxplot(df$var50~df$Stage,col="light blue")
title("metabolite 50")
boxplot(df$var49~df$Stage,col="light blue")
title("metabolite 49")
boxplot(df$var31~df$Stage,col="light blue")
title("metabolite 31")
boxplot(df$var57~df$Stage,col="light blue")
title("metabolite 57")
boxplot(df$var45~df$Stage,col="light blue")
title("metabolite 45")
dev.off()

boxplot(log2(df$var105)~df$Stage)
boxplot(log2(df$var1)~df$Stage)

# PART 1: LARGE SCALE INFERENCE

## 1.1 t test between 2 groups of patient: advanced and early stage

dim(df)[2]-2
head(df)

attach(df)
obs <- length(df[,1])
metabolites <- dim(df)[2]-2
risk.factor<- df[,2]
table(df$Stage)
risk.label <- ifelse(risk.factor=="early",0,1)

df$class[Stage=="advanced"] <- 1
df$class[Stage=="early"] <- 0

sd.i<-e.i<-t.t<-p.val<-c(1:metabolites)

for(i in 1:metabolites){	
  #cat(i)
  x1<-df[df$class==1,i+2]
  x2<-df[df$class==0,i+2]
  t.i<-t.test(x1,x2)#, alternative="two.sided")
  t.t[i]<-t.i$statistic
  p.val[i]<-t.i$p.value
  e.i[i]<-t.i$estimate[1]-t.i$estimate[2]
  sd.i[i]<-e.i[i]/t.t[i]
}

png("signal1.png")
hist(e.i,nclass=50,col="light green",xlab="signal",
     main="Difference between two stages")
dev.off()

plot(e.i,xlab="metabolites",ylab="fold change")

png("ttest1.png")
plot(e.i,t.t,col=4,cex=.8,xlab="fold change",ylab="t statistic test",
     main="Fold change and statistic test")
abline(v=-2,lty=3,lwd=3,col="orange")
abline(v=2,lty=3,lwd=3,col="orange")
dev.off()

png("ttest2.png")
plot(t.t,sd.i,xlab="t statistic", ylab="noise",cex=1.1,col=4,
     main="Statistic test and Noise")
dev.off()
png("ttest3.png")
plot(e.i,t.t,col=4,cex=1.1,xlab="signal",ylab="test statistic",
     main="Fold change and Statistic test")
dev.off()

png("pval_ttest.png")
hist(p.val,nclass = 30,col="light blue", xlab="p-value",
     main="Histogram of raw p-value")
text(0.25,10,paste("raw p-value < 0.05:",sum(p.val<.05),"metabolites"),
     cex=1,font = 4, col="dark green")
dev.off()

sum(p.val<0.05)

plot(log(p.val),e.i)

#### adjust for multiplicity
rawp <- p.val
cbind(sort(rawp))
holm<-mt.rawp2adjp(rawp, proc=c("Holm"))
bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
bh<-mt.rawp2adjp(rawp, proc=c("BH"))
by<-mt.rawp2adjp(rawp, proc=c("BY"))
allp<-cbind(rawp, bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2])
png("adjp.png")
mt.plot(allp,plottype="pvsr", proc=c("rawp","Bonferroni","Holm","BH","BY"),
        leg=c(70,0.5),lty=1,col=1:5,lwd=2,main="Adjusted p-value for t-test")
abline(h=.05,lty=4,col=2,lwd=2)
dev.off()
mt.reject(cbind(rawp,bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2]),seq(0,1,0.05))$r


## 1.2 one-way ANOVA
p.val<-c(1:metabolites)
for(i in 1:metabolites){	
  #cat(i)
  fit.aov<-aov(df[,i+2]~as.factor(df$class))
  p.val[i]<-anova(fit.aov)[1,5]
}

rawp<-p.val
holm<-mt.rawp2adjp(rawp, proc=c("Holm"))
bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
bh<-mt.rawp2adjp(rawp, proc=c("BH"))
by<-mt.rawp2adjp(rawp, proc=c("BY"))
allp<-cbind(rawp, bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2])
png("adjp1.png")
mt.plot(allp,plottype="pvsr", proc=c("rawp","Bonferroni","Holm","BH","BY"),
        leg=c(70,0.5),lty=1,col=1:5,lwd=2,main="Adjusted p-value for ANOVA")
abline(0.05,0,col="red",lty=3)
dev.off()
mt.reject(cbind(rawp,bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2]),seq(0,1,0.05))$r

## 1.3 Analysis with LIMMA

#metabolites<-dim(df)[2]-3
sd.i<-e.i<-t.t<-p.val<-c(1:metabolites)
t.tw<-p.valw<-c(1:metabolites)

## t test ##

for(i in 1:metabolites){	
  #cat(i)
  x1<-df[df$class==1,i+2]
  x2<-df[df$class==0,i+2]
  t.i<-t.test(x1,x2, alternative="two.sided")
  t.w<-wilcox.test(x1,x2) # wilcoxon test
  t.t[i]<-t.i$statistic
  p.val[i]<-t.i$p.value
  t.tw[i]<-t.w$statistic
  p.valw[i]<-t.w$p.value # p-value from wilcoxon
  e.i[i]<-t.i$estimate[1]-t.i$estimate[2]
  sd.i[i]<-e.i[i]/t.t[i]
}

par(mfrow=c(2,1))
hist(p.valw)
hist(p.val)



par(mfrow=c(2,2))
p.val.all<-p.val
e.i.all<-e.i
plot(e.i,-log(p.val.all))
title("a")
hist(p.val.all,col=0,nclass=50,main=" ")
title("b")
plot(t.t,e.i,xlab="t-test statisttic",ylab="log(fold change)")
abline(2,0,lty=4)
abline(-2,0,lty=4)
title("c")
plot(e.i,sd.i,xlab="log(fold change)",ylab="standard error")
lines(c(2,2),c(0,6),lty=4)
lines(c(-2,-2),c(0,6),lty=4)
title("d")

### adjust for multiplicity
library(Biobase)
library(multtest, verbose = FALSE)
rawp<-p.val
holm<-mt.rawp2adjp(rawp, proc=c("Holm"))
bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
bh<-mt.rawp2adjp(rawp, proc=c("BH"))
by<-mt.rawp2adjp(rawp, proc=c("BY"))
allp<-cbind(rawp, bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2])
par(mfrow=c(1,1))
mt.plot(allp,plottype="pvsr", proc=c("rawp","Bonferroni","Holm","BH","BY"),leg=c(2100,0.375),lty=c(1,2,3,4,5),col=1:5,lwd=2)
mt.plot(allp,plottype="pvsr", proc=c("rawp","Bonferroni","Holm","BH","BY"),leg=c(-100,-100),lty=c(1,2,3,4,5),col=1:5,lwd=2)
abline(0.05,0, col="red")
text(40,0.4,"raw p values",adj=0,cex=0.75)
text(10,0.8,"BH",adj=0,cex=0.75)
text(5,0.95,"BY",adj=0,cex=0.75)
text(4,0.9,"Bonf. & Holm",adj=0,cex=0.75)

mt.reject(cbind(rawp,bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2]),seq(0,1,0.05))$r


######################### LIMMA application 
library(limma)
par(mfrow=c(1,1))

mtbl.exp<-data.frame(t(df[,c(3:101)]))
nn<-length(df$class)
nn

design <- cbind(rep(1,nn),abs(df$class-1))
design

fit <- lmFit(mtbl.exp,design)
summary(fit)

fit1 <- eBayes(fit)
head(fit1$t)
head(fit1$p.value)

topTable <- topTable(fit1,number=nrow(mtbl.exp),coef=2)
head(topTable)
dim(topTable)
dim(topTable[topTable[,4]<=0.05,])
dim(topTable[topTable[,5]<=0.05,])


help(topTable)

length(topTable[,4])
png("limma1.png")
plot(t.t,-log(p.val),col=4,cex=.8,main="T-test versus LIMMA")
points(topTable[,3],-log(topTable[,4]),col=2,cex=.7)
legend(-2,5,legend = c("t-test","limma"),col=c(4,2),pch=c(1,1))
dev.off()

y <- df$class
x <- as.matrix(df[,c(3:101)])
m1 <- glm(y~x)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#               SAM - modified test statistic                   #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(samr)
y<-df$class+1
data<-data.frame(t(df[,c(3:101)]))
dim(data)[1]
d=list(
  x=as.matrix(data),y=y,
  metaid=as.character(1:dim(data)[1]),
  metanames=paste("metabolite",
                  as.character(1:dim(data)[1]))
  ,logged2=FALSE)

samr.obj <- samr(d,resp.type="Two class unpaired")

delta.table <- samr.compute.delta.table(samr.obj)

par(mfrow=c(1,1))
delta.table

png("sam.png")
samr.plot(samr.obj,1.5)
title("SAM analysis for significant features")
dev.off()

?samr.plot

sigmeta.table<-samr.compute.siggenes.table(samr.obj,del=1.5,
                              data=d, delta.table)

sigmeta.table
t.t <- samr.obj$tt
signal <- samr.obj$foldchange
noise<-samr.obj$sd

plot(signal,t.t,xlim=c(-3,3))

p.val.sam <- pt(q = t.t,df = 148)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#             FEATURES FILTERING - BIOCONDUCTOR                 #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
library(genefilter)
  
iqr<-c(1:metabolites)
for(i in 1:metabolites){
  iqr[i]<-quantile(df[,i+2],prob=c(0.75))-
    quantile(df[,i+2],prob=c(0.25))
}
iqr1<-quantile(iqr,prob=seq(from=0,to=1,by=0.05))
plot(c(1:metabolites),iqr)
abline(0.39,0)

f2 <- function(x)(IQR(x) > 0.31)

ff <- filterfun(f2)

selected <- genefilter(t(data.frame(df[,c(2:103)])),ff)
sum(selected)
index<-c(1:ngenes)
index1<-index[selected]
length(index1)
data2<-data1[index1,] # filter in
dim(data2)
#################################################################
#                     CLASSIFICATION                            #
#################################################################

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#                       SUPERVISED GLM                          #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(glm2)

t.t<-beta.i<-p.val<-c(1:metabolites)
for(i in 1:metabolites){	
  #cat(i)
  xi<-df[,i+2]
  modi<-glm(risk.label~xi,family = "binomial"(link="logit"))
  t.t[i]<-summary(modi)$coefficients[2,3]
  beta.i[i]<-summary(modi)$coefficients[2,1]
  p.val[i]<-summary(modi)$coefficients[2,4]
}
hist(p.val,nclass=50)
sum(p.val<0.05)
sum(p.val<0.1)

plot(beta.i[beta.i<20],t.t[beta.i<20],xlab="beta",ylab="test statistic")
plot(t.t[beta.i<20],-log(p.val[beta.i<20]),xlab="test statistic",ylab="-log(P)")

plot(beta.i[beta.i<20],-log(p.val[beta.i<20]),xlab="beta",ylab="-log(P)")
p.val.sort<-sort(p.val)
beta.i.sort<-beta.i[order(p.val)]
points(beta.i.sort[1:500],-log(p.val.sort[1:500]),col=2)

which(p.val.sort<.05)

df.glm<-df[,which(p.val<.05)+2]
pz<-c(1:6)


png("glm.png")
par(mfrow=c(2,3))
for(i in 1:6){	
  #cat(i)
  xi<-df.glm[,i]
  modi<-glm(risk.label~xi,family = "binomial"(link="logit"))
  #t.t[i]<-summary(modi)$coefficients[2,3]
  #beta.i[i]<-summary(modi)$coefficients[2,1]
  #p.val[i]<-summary(modi)$coefficients[2,4]
  y<-modi$fitted.values
  plot(df.glm[,i],risk.label,
       xlab="expression level")
  lines(df.glm[,i],y,col=risk.factor,cex=.75,ylim=c(0,1))
  abline(h=.5,lty=3,col=4)
  #abline(v=mean(df.glm[,i]),lty=3,col=4)
  title(paste("metabolite",names(df.glm)[i]))
}
dev.off()

m1<-glm(risk.label~df[,24],family = "binomial"(link="logit"))
plot(df.glm$var31,risk.label)
lines(df.glm$var31,predict.glm(m1,list(xi = df.glm$var31)),type="response")
y.p<-exp(m1$coefficients[1]+m1$coefficients[2]*df.glm$var31)/(1+exp(m1$coefficients[1]+m1$coefficients[2]*df.glm$var31))
lines(df.glm$var31,y.p)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#                       SUPERVISED PCA                          #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

###### PCA for full data
pc.cr.f <- princomp(df[,c(3:101)],scores=TRUE,cor = TRUE)
names(pc.cr.f)
x<-prcomp(df[,c(3:101)])
names(x)

png("pca1.png")
plot(pc.cr.f,main="PCA with all metabolic signature",col="light blue")
dev.off()
plot(pc.cr.f$sdev^2/sum(pc.cr.f$sdev^2),#type="p",pch=20,
     col=4,cex=1,ylab="proportions of variance",xlab="Component",
     main="Proportions of variance by component")
lines(pc.cr.f$sdev^2/sum(pc.cr.f$sdev^2),lty=2,col="dark blue")
#points(pc.cr.f$sdev^2/sum(pc.cr.f$sdev^2),ltw=2,
#       col=3,pch=20)
dev.off()

###### PCA for reduced data
df.R <- df[,which(p.val<.05)+2]#df.glm

pc.cr.r <- princomp(df.R,scores=TRUE,cor = TRUE)
names(pc.cr.r)
plot(pc.cr.r,main="PCA with all metabolic signature",col="light blue")

library(mdatools)
model = pca(df[,c(3:101)], scale = TRUE, cv = 5, info = 'Simple PCA model')
model = selectCompNum(model, 1)

#
res = model$calres
summary(res)

png("pca1.png")
plot(res)
dev.off()

png("scorepca.png")
plotScores(res, cgroup = df[, index1[3]], show.labels = TRUE)
dev.off()
#

summary(model)
png("pca1.png")
plot(model, show.labels = TRUE)
dev.off()

UX1<-pc.cr.f$scores[,1]
m1.pca<-glm(risk.label~UX1,family = "binomial"(link="logit"))

y<-m1.pca$fitted.values
plot(UX1,risk.label, xlab="Fisrt component")
lines(UX1,y,col=risk.factor,cex=.75,ylim=c(0,1))
abline(h=.5,lty=3,col=4)
title(paste("metabolite",names(df.glm)[i]))

######## reduced data
df.R <- df[,which(p.val<.05)+2]
model.r = pca(df.R, scale = TRUE, cv = 10, info = 'Simple PCA model')
model.r = selectCompNum(model, 1)
#
res = model.r$calres
summary(res)

png("pca1.png")
plot(res)
dev.off()

png("scorepca.png")
plotScores(res, cgroup = df[, index1[3]], show.labels = TRUE)
dev.off()
#

summary(model.r)
png("pca2.png")
plot(model.r, show.labels = TRUE)
dev.off()

ux1<-pc.cr.r$scores[,1]
m2.pca<-glm(risk.label~ux1,family = "binomial"(link="logit"))

png("pca4.png")
y<-m2.pca$fitted.values
plot(ux1,risk.label, xlab="Fisrt component",
     main="Classification by first PCA",font=2,cex=.8)
lines(ux1,y,col=risk.factor,cex=.75,ylim=c(0,1))
abline(h=.5,lty=3,col=4)
dev.off()

## cross validation
png("cv1.png")
plotVariance(model, type ='h',main="Cross validation with full metabolite profiles")
dev.off()
png("cv2.png")
plotVariance(model.r,type='h',main="Cross validation with significant metabolite profiles")
dev.off()
##############################################

index<-c(1:metabolites)
index1<-index[order(p.val)]
p.val.sort<-sort(p.val)
kk<-index1[1]
xk<-df[,kk]
modk<-glm(risk.label~-1+xk,binomial(link = "logit"))
summary(modk)
plot(xk,risk.label)
lines(sort(xk),modk$fit[order(xk)],lwd=3,col=2)
title(paste(kk))

K<-6
p.val.beta1<-p.val.sort
cutoff<-sort(p.val)[K]
sum(p.val.beta1 <= cutoff)
min(p.val.beta1)
index<-c(1:metabolites)

index<-c(1:metabolites)
index1<-index[order(p.val)]
index2<-which(p.val<.05)+2#index1[3:(K+2)]
index2


dim(df)
data.R<-df[,index2]
dim(data.R)

png("cor.png")
pairs(data.R,col="light blue")
dev.off()

ux1<-as.factor(pc.cr.r$scores[,1])
ux2<-as.factor(pc.cr.r$scores[,2])

length(ux1)

fit<-glm2(risk.label~ux1+ux2,family = "binomial"(link="logit"))
summary(fit)
plot(ux1,risk.label,
     xlab="expression level")
lines(df.R$var91,fit$fitted.values)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#                     TREE-BASED ANALYSIS                       #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

install.packages("tree")
library(tree)
risk.label
dim(data.R)

meta91<-data.R[,1]
meta50<-data.R[,2]
meta49<-data.R[,3]
meta31<-data.R[,4]
meta57<-data.R[,5]
meta45<-data.R[,6]
#plot(g1413,g2663)

tree1.fit <- tree(as.factor(risk.label)~meta91+meta50+meta49
                  +meta45+meta57+meta31)

tree1.fit
summary(tree1.fit)
png("tree1.png")
plot(tree1.fit,col="dark blue",lwd=2)
text(tree1.fit,cex=.6)
title("Tree-based analysis")
dev.off()

partition.tree(tree1.fit, label = "yval", add = FALSE)
misclass.tree(tree1.fit, detail = FALSE)


## split training and test set
?tree()
n=dim(data.R)[1]
id<-sample(n,size = floor(n*.7),replace = F)
train<-data.R[id,]
test<-data.R[-id,]
tree1 <- tree(as.factor(risk.label[id])~train$var31+train$var45+train$var49+
                train$var50+train$var57+train$var91)
predict(tree1,test)
misclass.tree(tree1)
#prune.tree(tree1, k = NULL, best = NULL, test, 
#           method = c("deviance", "misclass"), loss, eps = 1e-3)
#  TREE with PC1 & PC2                          #
pc.cr <- princomp((data.R),scores=TRUE,cor = TRUE)
loadings(pc.cr)
plot(pc.cr)
biplot(pc.cr)
summary(pc.cr)
Ux1<-as.vector(pc.cr$scores[,1])
Ux2<-as.vector(pc.cr$scores[,2])

png("pca3.png")
plot(Ux1,Ux2,pch=20,col=3,main="Classification by first 2 PCA",
     xlab="First Component",ylab="Second Component")
points(Ux1[risk.label==1],Ux2[risk.label==1],col=2,pch=20)
legend(-5.5,3,legend = c('advanced','early'),pch=c(20,20),col=c(2,3))
dev.off()

tree2.fit <- tree(as.factor(risk.label)~Ux2+Ux1)
tree2.fit
summary(tree2.fit)

png("tree2.png")
plot(tree2.fit,col="dark blue",lwd=2)
text(tree2.fit,cex=.7)
title("Tree-based analysis by first 2 PCA")
dev.off()

partition.tree(tree2.fit)
misclass.tree(tree2.fit, detail = FALSE)
tree3.fit <- tree(as.factor(risk.label)~df.R$var91)
plot(tree3.fit)
misclass.tree(tree3.fit)

plot(Ux1,Ux2)
points(Ux1[risk.label==0],Ux2[risk.label==0],col=2)
lines(c(-1.14495,-1.14495),c(-5,5),col=2,lwd=2)

plot(Ux1,xk)
points(Ux1[risk.label==0],xk[risk.label==0],col=2)
cor(Ux1,xk)

plot(Ux1,risk.label)
points(xk,risk.label,pch="+",col=2)


z<-cbind(Ux1,Ux2)
k1<-kmeans(z, 2)
k1
summary(k1)
plot(z, col = k1$cluster)
points(k1$centers, col = c(4,5), pch ="+", cex = 3)

k2<-kmeans(df.R,2)
plot(df.R, col = k2$cluster, pch=20)

k3<-kmeans(df.R$var91,2)
plot(df.R$var91,col=k3$cluster,pch=20)
#points(k3$cluster,col=2,pch=20)

# tree-based by specific metabolite
mis<-c(1:6)
par(mfrow=c(2,3))
for (i in 1:6){
  tree.fit <- tree(as.factor(risk.label)~df.R[,i])
  mis[i]<-misclass.tree(tree.fit, detail = FALSE)
  plot(tree.fit)
  text(tree.fit,cex=.7)
  title(names(df.R)[i])
}
tree.fit$frame

mis

plot(tree.fit)
text(tree.fit,cex=.7)
title(names(df.R)[i])

# Part 3: LDA for the Golub data using PC1 & PC2                    #
pc.cr <- princomp((data.R),scores=TRUE,cor = TRUE)
loadings(pc.cr)
plot(pc.cr)
biplot(pc.cr)
summary(pc.cr)
Ux1<-as.vector(pc.cr$scores[,1])
Ux2<-as.vector(pc.cr$scores[,2])
x<-Ux1
y<-Ux2
plot(x,y,xlab="U1",ylab="U2")
points(x[risk.label==0],y[risk.label==0],col=2)

xmat<-cbind(x,y)
S<-cov(cbind(x,y))
S1<-solve(S)

x1<-x[risk.label==0]
y1<-y[risk.label==0]
x2<-x[risk.label==1]
y2<-y[risk.label==1]

mx1<-c(mean(x1),mean(y1))
mx2<-c(mean(x2),mean(y2))
plot(x,y)
points(x[risk.label==0],y[risk.label==0],col=2)
points(mx1[1],mx1[2],pch="+",lwd=5,cex=3)
points(mx2[1],mx2[2],pch="+",lwd=5,cex=3)


xvec<-mx2-mx1
w<-t(xvec)%*%S1
w
theta<-w%*%(mx1+mx2)/2
theta
z<-w[1]*x+w[2]*y
hist(z,nclass=50,main="histogram of z")
sort(z)
lines(c(theta,theta),c(0,10),col=2,lwd=2)

png("lda2.png")
plot(z,Ux1,pch=20,col=3,main="LDA with first PCA")
points(z[risk.label==0],Ux1[risk.label==0],col=2,pch=20)
#lines(c(theta,theta),c(-5,10),col=2,lwd=2)
dev.off()
png("lda3.png")
plot(z,Ux2,pch=20,col=3,main="LDA with second PCA")
points(z[risk.label==0],Ux2[risk.label==0],col=2,pc=20)
#lines(c(theta,theta),c(-5,10),col=2,lwd=2)
dev.off()
# Part 3: LDA for the Golub data using 6 significant   

xmat<-data.R
S.r<-cov(xmat)
S1.r<-solve(S.r)

x1<-x[risk.label==0]
y1<-y[risk.label==0]
x2<-x[risk.label==1]
y2<-y[risk.label==1]

mx1<-c(mean(xmat$var31[risk.label==0]),mean(xmat$var45[risk.label==0]),
       mean(xmat$var49[risk.label==0]),mean(xmat$var50[risk.label==0]),
       mean(xmat$var57[risk.label==0]),mean(xmat$var91[risk.label==0]))
mx2<-c(mean(xmat$var31[risk.label==1]),mean(xmat$var45[risk.label==1]),
       mean(xmat$var49[risk.label==1]),mean(xmat$var50[risk.label==1]),
       mean(xmat$var57[risk.label==1]),mean(xmat$var91[risk.label==1]))
plot(x,y)
points(x[risk.label==0],y[risk.label==0],col=2)
points(mx1[1],mx1[2],pch="+",lwd=5,cex=3)
points(mx2[1],mx2[2],pch="+",lwd=5,cex=3)


xvec<-mx2-mx1
w<-t(xvec)%*%S1.r
w
theta<-w%*%(mx1+mx2)/2
theta
z<-w[1]*xmat$var31+w[2]*xmat$var45+w[3]*xmat$var49+
  w[4]*xmat$var50+w[5]*xmat$var57+w[6]*xmat$var91
hist(z,nclass=50,main="histogram of z")
sort(z)
lines(c(theta,theta),c(0,10),col=2,lwd=2)


plot(z,xmat$var31)
points(z[risk.label==0],xmat$var31[risk.label==0],col=2)
lines(c(theta,theta),c(-5,10),col=2,lwd=2)


plot(z,Ux2)
points(z[risk.label==0],Ux2[risk.label==0],col=2)
lines(c(theta,theta),c(-5,10),col=2,lwd=2)
png("lda1.png")
par(mfrow=c(2,3))
for (i in 1:6){
  plot(z,xmat[,i],pch=20,col=3,ylab=names(xmat)[i])
  points(z[risk.label==0],xmat[,i][risk.label==0],col=2,pch=20)
  title(paste("LDA by",names(xmat)[i]))
}
dev.off()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#                     CLASSIFICATION CMA                        #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("CMA", version = "3.8")
library(CMA)
library(e1071)

dim(df)
df1<-df[,c(2,3:101)]
Y <- df1[,1]
X <- as.matrix(df1[,-1])
dim(X)

##generating learning sets by cross validation 
mccv <- GenerateLearningsets(y=Y, method = "MCCV", niter=10, ntrain=30)
loo <- GenerateLearningsets(y=Y, method="LOOCV")

### 1. gene selection: t test
selttest.mccv <- GeneSelection(X, Y, learningsets =mccv, method ="t.test")
selttest.loo <- GeneSelection(X, Y, learningsets =loo, method = "t.test")


show(selttest.mccv)
toplist(selttest.mccv, k = 10, iter = 1)

show(selttest.loo)
toplist(selttest.loo, k = 10, iter = 1)

png("iter_sel.png")
par(mfrow=c(3,2),cex=.5)
plot(selttest.mccv, iter = 1,xlab="metabolite",ylab="importance",
     main="Variable importance Iter 1 - MCCV")
plot(selttest.loo, iter = 1,xlab="metabolite",ylab="importance",
     main="Variable importance Iter 1 - LOOCV")
plot(selttest.mccv, iter = 2,xlab="metabolite",ylab="importance",
     main="Variable importance Iter 2 - MCCV")
plot(selttest.loo, iter = 2,xlab="metabolite",ylab="importance",
     main="Variable importance Iter 2 - LOOCV")
plot(selttest.mccv, iter = 3,xlab="metabolite",ylab="importance",
     main="Variable importance Iter 3 - MCCV")
plot(selttest.loo, iter = 3,xlab="metabolite",ylab="importance",
     main="Variable importance Iter 3 - LOOCV")
dev.off()
names(df)[57]

######################################
########### DLDA (LOOCV) #############
######################################


mccv <- GenerateLearningsets(y=Y, method = "MCCV", niter=100, ntrain=30)
loo <- GenerateLearningsets(y=Y, method="LOOCV")

dlda <- classification(X, Y, learningsets = loo,
                       genesellist = list(method  = "t.test"), classifier = dldaCMA,
                       nbgene = 20)

MCE <- evaluation(dlda,measure = c("misclassification"))
SPE <- evaluation(dlda,measure = c("specificity"))
SEN <- evaluation(dlda,measure = c("sensitivity"))
dldajoin <- join(dlda)
MCE
SPE
SEN

show(MCE)

attributes(dldajoin)

y1<-attributes(dldajoin)$y
length(y1)

y2<-attributes(dldajoin)$yhat
length(y2)

plot(y1)
points(y2,col=2)

table(y1,y2)

plot(attributes(dldajoin)$prob)

######################################
########### DLDA (MCCV)  #############
######################################


mccv <- GenerateLearningsets(y=Y, method = "MCCV", niter=100, ntrain=30)
loo <- GenerateLearningsets(y=Y, method="LOOCV")

dlda <- classification(X, Y, learningsets = mccv,
                       genesellist = list(method  = "t.test"), classifier = dldaCMA,
                       nbgene = 20)


MCE <- evaluation(dlda,measure = c("misclassification"))
eva.dlda <- evaluation(dlda,measure = c("misclassification"))
SPE <- evaluation(dlda,measure = c("specificity"))
SEN <- evaluation(dlda,measure = c("sensitivity"))
dldajoin <- join(dlda)
MCE
SPE
SEN

show(eva.dlda)
par(mfrow=c(2,2))
boxplot(MCE)
boxplot(SPE)
boxplot(SEN)

boxplot(eva.dlda)
eva.dlda
attributes(dldajoin)

y1<-attributes(dldajoin)$y
length(y1)

y2<-attributes(dldajoin)$yhat
length(y2)

table(y1,y2)

attributes(dldajoin)$prob


######################################
########### LDA  #####################
######################################


mccv <- GenerateLearningsets(y=Y, method = "MCCV", niter=100, ntrain=30)
loo <- GenerateLearningsets(y=Y, method="LOOCV")
?GenerateLearningsets

lda <- classification(X, Y, learningsets = mccv,
                      genesellist = list(method  = "t.test"), classifier = ldaCMA,
                      nbgene = 20)

lda <- classification(X, Y, learningsets = loo,
                      genesellist = list(method  = "t.test"), classifier = ldaCMA,
                      nbgene = 20)

eva.lda <- evaluation(lda,measure = "misclassification")
ldajoin <- join(lda)
show(eva.lda)
boxplot(eva.lda)

attributes(ldajoin)

y1<-attributes(ldajoin)$y
length(y1)

y2<-attributes(ldajoin)$yhat
length(y2)
table(y1,y2)

attributes(dldajoin)$prob


library(e1071)
svm <- classification(X, Y, learningsets = mccv,
                      genesellist = list(method  = "t.test"), classifier = svmCMA,
                      nbgene = 20)

eva.svm <- evaluation(svm,measure = "misclassification")
svmjoin <- join(svm)

boxplot(eva.svm)
attributes(eva.svm)
y<-attributes(eva.svm)$yhat
z<-attributes(eva.svm)$y

png("perf.png")
par(mfrow=c(1,3))
boxplot(eva.lda,main="LDA",col="light blue")
boxplot(eva.dlda,main="DLDA",col="light blue")
boxplot(eva.svm,main="SVM",col="light blue")
dev.off()

#######################################
########### loop ######################
#######################################

mccv <- GenerateLearningsets(y=Y, method = "MCCV", niter=100, ntrain=30)
loo <- GenerateLearningsets(y=Y, method="LOOCV")
genesel <- c("t.test","wilcox.test","f.test")
meta <- c(1,3,6,8,10)
MB<-length(meta)
classmethod <- c("knnCMA","dldaCMA","plsCMA","svmCMA")
eva.knn <- eva.dlda <- eva.lda <- eva.pls <- eva.svm <-matrix(0,1,MB)

for (i in 1:1){
  for(j in 1:MB) {
    #result.knn <- classification(golubX, golubY, learningsets = loo,
    #                       genesellist = list(method = genesel[i] ), classifier = knnCMA,
    #                        nbgene = ngene[j])
    #a <- evaluation(result.knn,measure = "misclassification")
    #eva.knn[i,j] <- mean(attributes(a)$score)
    
    
    result.lda <- classification(X, Y, learningsets = loo,
                                 genesellist = list(method = genesel[i] ), classifier = ldaCMA,
                                 nbgene = ngene[j])
    a <- evaluation(result.lda,measure = "misclassification")
    
    eva.lda[i,j] <- mean(attributes(a)$score)
    
    
    
    result.dlda <- classification(X, Y, learningsets = loo,
                                  genesellist = list(method = genesel[i] ), classifier = dldaCMA,
                                  nbgene = ngene[j])
    a <- evaluation(result.dlda,measure = "misclassification")
    
    eva.dlda[i,j] <- mean(attributes(a)$score)
    
    result.svm <- classification(X, Y, learningsets = loo,
                                 genesellist = list(method = genesel[i] ), classifier = svmCMA,
                                 nbgene = ngene[j])
    a <- evaluation(result.svm,measure = "misclassification")
    eva.svm[i,j] <- mean(attributes(a)$score)
  }
}

help(evaluation)

show(eva.dlda)
show(eva.lda)

par(mfrow=c(1,1))
plot(meta,eva.dlda,type="l",ylab="MCA",xlab="Number of metabolites")#,ylim=c(0,0.15))
lines(meta,eva.lda,col=2)
legend(5,0.48,c("DLDA","LDA"),lty=c(1,1),col=c(1,2))
help(legend)

plot(meta,eva.svm)

help(classification)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#                      UNSUPERVISED PCA                         #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
install.packages("mpm")
library(mpm)

############################################
############ spectral map with PCA #########
############################################
df1<-df[,c(2:101)]
par(mfrow=c(1,1))
r.pca <- mpm(data.frame(t(df1)), center = "column", normal = "column")
r <- plot(r.pca, label.tol = 20, scale = "uvc", col.group = (Golub.grp)[1:38], zoom = c(1,1.2), col.size = 5)
?mpm

############################################
############ spectral map###################
############################################

r.sma <- mpm(data.frame(t(df1)), row.weight = "mean", col.weight = "mean")
r <- plot(r.sma, label.tol = 20, scale = "uvc",
          col.group = (Golub.grp)[1:38], zoom = c(1,1.2), col.size = 5)
Golub[r$Rows$Select==1, 1]


############################################
#                                          #
#                                          #
# Part 1: clustering                       #
#                                          #
#                                          #
############################################

a <- dist(t(df1[,-1]))
a1 <- dist(t(df1[,-1]), method="euclidean",diag=TRUE, upper=TRUE)
b <- as.dist(cor(df1[,-1]))

metabolites<-dim(df)[2]-3
## t test ##
names(df)
for(i in 1:metabolites){	
  xi<-df[,i+2]
  x1<-xi[risk.label==1]
  x2<-xi[risk.label==0]
  t.i<-t.test(x1,x2, alternative="two.sided")
  t.t[i]<-t.i$statistic
  e.i[i]<-t.i$estimate[1]-t.i$estimate[2]
  p.val[i]<-t.i$p.value
}
hist(p.val,nclass=50)

which(p.val<0.05)
par(mar=c(4,4,4,4))
plot(t.t,-log(p.val))
plot(e.i,-log(p.val),pch=20,col=3)
abline(h=2,lty=2)
abline(v=c(-1,1),lty=2)
#sum(p.val<0.1)
mask <- with(t.t, abs(se.i) < .2 & p.value < .01)
png("cl1.png")
plot(t.t,p.val,pch=20,col=3,xlab="t statistic",ylab="p-value",
     main="t statistic and p-value by metabolites")
abline(0.05,0,col=2,lty=4,lwd=2)
dev.off()

index<-c(1:metabolites)
index1<-index[order(p.val)]
p.val.sort<-sort(p.val)

############################################
## reduced matrix                        ###
############################################

index2<-index1#[1:100]


data.g<-df1[-1,]
data.g<-df1[,index2]
dim(data.g) 


image(as.matrix(data.g),c(1:dim(data.g)[1]),c(1:dim(data.g)[2]),xlab="conditions",ylab="features",yaxt="n")
?image()

############################################
## clustering                            ###
############################################
a1 <- dist(df[,c(3:101)], method="euclidean",diag=TRUE, upper=TRUE)
hc1 <- hclust(a1, method="ave")
par(mfrow=c(1,1))
png("cl2.png")
plot(hc1,cex=0.35,xlab="metabolite",lty=1,lwd=1,col="dark blue",
     labels=risk.label)
dev.off()
#golub[,1]
hclusters <- cutree(hc1, h=70)
table(true=risk.label, cluster=hclusters)


hc2 <- hclust(b, method="ave")

par(mfrow=c(1,1))
plot(hc1)

install.packages("gplots")
library(gplots)
#x11()
par(mar=c(2,3,3,3))
dend1 <- as.dendrogram(hc1) 
png("heat.png")
heatmap.2(as.matrix(a1), 
          col=colorRampPalette(c("blue4", "white"))(10),
          keysize=0.85,scale="none", margins=c(5,3), 
          density.info="none",
          Colv=dend1, Rowv=dend1,
          dendrogram="none",
          trace="none", 
          labRow=" ",
          cexRow=0.25, cexCol=0.45,  
          ylab='',
          xlab="features",
          symkey=FALSE
)
dev.off()

install.packages("rafalib")
library(rafalib)
mypar()
d <- dist(df[,which(p.val<.05)+2])
hc <- hclust(d)
plot(hc,labels=risk.label,cex=0.5)
png("cl3.png")
myplclust(hc, labels=risk.label,main="Cluster by 6 significant metabolites",
          lab.col=as.fumeric(c("advanced","early")), cex=0.5)
dev.off()
hclusters <- cutree(hc, h=20)
table(true=risk.label, cluster=hclusters)

d1 <- dist(df[,c(3:101)])
hc1 <- hclust(d1)
plot(hc1,labels=risk.label,cex=0.5)
myplclust(hc1, labels=risk.label, lab.col=as.fumeric(c("advanced","early")), cex=0.5)
hclusters <- cutree(hc1, h=160)
table(true=risk.label, cluster=hclusters)

d2 <- dist(cbind(Ux1,Ux2))
hc2 <- hclust(d2)
plot(hc2,labels=risk.label,cex=0.5,lwd=1,col="dark blue")
png("cl4.png")
myplclust(hc2, labels=risk.label,lwd=1,main="Cluster by 2 PCAs",
          lab.col=as.fumeric(c("advanced","early")), cex=0.5)
dev.off()
hclusters <- cutree(hc2, h=7)
table(true=risk.label, cluster=hclusters)

# kmeans
set.seed(1)
km1 <- kmeans(df[,which(p.val<.05)+2], centers=5)
names(km)
png("kmean1.png")
plot(km1$cluster,pch=16,col=as.fumeric(c("advanced","early")),
     xlab="patient",ylab="cluster",
     main="Classification with 6 metabolites")
dev.off()
km1$withinss
set.seed(2)

km2 <- kmeans(cbind(Ux1,Ux2), centers=2)
png("kmean2a.png")
plot(Ux1,Ux2,col=as.fumeric(c("advanced","early")),pch=16,
     main="Original stage")
dev.off()
png("kmean2b.png")
plot(Ux1,Ux2,col=km2$cluster,pch=16,
     main="Classification by apply kmeans")
dev.off()

km3 <- kmeans(df[,c(3:101)], centers=2)
names(km)
png("kmean3.png")
plot(km3$cluster,pch=16,col=as.fumeric(c("advanced","early")),
     xlab="patient",ylab="cluster",
     main="Classification with all metabolites")
dev.off()


cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(rownames(df[,c(3:101)]))]
library(genefilter)
library(gplots)
library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
head(cbind(colnames(t(df[,c(3:101)]))),cols)

rv <- rowVars(t(df[,c(3:101)]))
idx <- order(-rv)[1:99]
heatmap.2(t(df[,c(3:101)]), labCol=metabolites,
          trace="none")#,ColSideColors=cols)#,col=hmcol)



km <- kmeans(df[,which(p.val<.05)+2][,1], centers=2)
names(km)
png("kmean1.png")
plot(df[,which(p.val<.05)+2][,1],km$cluster,pch=16,col=as.fumeric(c("advanced","early")),
     xlab="patient",ylab="cluster",
     main="Classification with 6 metabolites")
dev.off()
#########################################################
#       STANDARDIZE VALUE OF METABOLIC EXPRESSION       #
#########################################################
meta31<-(meta31-mean(meta31))/sd(meta31)
meta45<-(meta45-mean(meta45))/sd(meta45)
meta49<-(meta49-mean(meta49))/sd(meta49)
meta50<-(meta50-mean(meta50))/sd(meta50)
meta57<-(meta57-mean(meta57))/sd(meta57)
meta91<-(meta91-mean(meta91))/sd(meta91)

df.s<-cbind(meta31,meta45,meta49,meta50,meta57,meta91)
fit.s<-glm(risk.label~meta31+meta45+meta49+meta50+meta57+meta91,
           family = "binomial"(link="logit"))
summary(fit.s)
km.s <- kmeans(df.s, centers=2)
plot(km$cluster,pch=16,col=as.fumeric(c("advanced","early")),
     xlab="patient",ylab="cluster",
     main="Classification with 6 metabolites")
t.test(meta31[risk.label==1],meta31[risk.label==0])
t.test(meta91[risk.label==1],meta91[risk.label==0])
t.test(meta45[risk.label==1],meta45[risk.label==0])

ds <- dist(df.s)
hcs <- hclust(ds)
plot(hcs,labels=risk.label,cex=0.5)
myplclust(hcs, labels=risk.label, lab.col=as.fumeric(c("advanced","early")), cex=0.5)
hclusters <- cutree(hcs, h=9)
table(true=risk.label, cluster=hclusters)

pca1<-princomp(df.s)
plot(pca1)
s1<-pca1$scores[,1]
s2<-pca1$scores[,2]
kms1 <- kmeans(cbind(s1,s2), centers=2)
plot(s1,s2,col=as.fumeric(c("advanced","early")),pch=16,
     main="Original stage")
plot(s1,s2,col=kms1$cluster,pch=16,
     main="Classification by apply kmeans")





########### LASSO
library(glmnet)
y<-risk.label
x <- as.matrix(df[,c(3:101)])
x<-df.s
lasso.fit1 <- glmnet(x,y,alpha = 1,family = "binomial")
cv.lasso.fit1 <- cv.glmnet(x,y,alpha=1)
plot(cv.lasso.fit1)
plot(lasso.fit1, xvar="lambda")
abline(v=log(cv.lasso.fit1$lambda.min),col="blue",lwd=3)
cv.lasso.fit1$lambda.min
log(cv.lasso.fit1$lambda.min)

plot(lasso.fit1,xvar="lambda") 
abline(v=log(cv.lasso.fit1$lambda.min),col="blue",lwd=3)

set.seed(2)
#y<-df.R$SUV_max
x <- as.matrix(df[,which(p.val<.05)+2])
lasso.fit2 <- glmnet(x,y,alpha = 1,family = "binomial")
cv.lasso.fit2 <- cv.glmnet(x,y,alpha=1)

png("lamda3.png")
plot(cv.lasso.fit2)
dev.off()

plot(lasso.fit2,xvar="lambda") 
abline(v=log(cv.lasso.fit2$lambda.min),col="blue",lwd=3)

######## Elastic Net with alpha = 0.5 ######
set.seed(2)
#y<-df$SUV_max
x <- as.matrix(df[,c(3:101)])
Elnetfit1 <- glmnet(x,y,alpha=.7, family ="binomial")
cv.Elnetfit1<-cv.glmnet(x,y,alpha=.7,nfolds = 10)
plot(Elnetfit1,xvar="lambda")
abline(v=log(cv.Elnetfit1$lambda.min),col="blue",lwd=3)

cv.Elnetfit1.bin <- cv.glmnet(x,y,alpha=0.5, family ="binomial",nfolds = 6)
plot(cv.Elnetfit1.bin)


rpart.model <- rpart(risk.label ~ df.s, method = "anova")
rpart.model
plot(rpart.model)

### prediction with lasso,should be new X !!  ###
###if s = ... not specified then predicted for the entire grid of lambda used to build the model ##

### prediction error ###

y.cv.predict <- predict.cv.glmnet(cv.lasso.fit1, x, s="lambda.min")
plot(y,y.cv.predict,col=risk.label)
png("overfit_lasso.png")
plot(y,xlab="Observation",ylab="Scaled SUV max",
     main="Overfit checking with LASSO prediction")
lines(y.cv.predict,col=2)
legend(0,5,col=c(1,2),pch = c(1,1),legend=c("Observed","Predicted"))
dev.off()
abline(0,1)

data1<-t(df[,c(3:101)])
iqr<-c(1:metabolites)
for(i in 1:metabolites)
{
  iqr[i]<-quantile(t(df[,c(3:101)])[i,],prob=c(0.75))-
    quantile(t(df[,c(3:101)])[i,],prob=c(0.25))
}
iqr1<-quantile(iqr,prob=seq(from=0,to=1,by=0.05))

png("filter.png")
plot(c(1:metabolites),iqr,pch=16,col="light green",
     xlab="metabolites",ylab="IQR", main="IQR by metabolite")
lines(c(1:metabolites),iqr,col="light green")
abline(10,0,lty=2,col=4)
text(10,12,"Selected",col=2)
dev.off()

f2 <- function(x) (IQR(x) > 10)
ff <- filterfun(f2)

selected <- genefilter(data1,ff)
sum(selected)
index<-c(1:metabolites)
index1<-index[selected]
length(index1)
data2<-data1[index1,]

a2 <- dist(t(data2), method="euclidean",diag=TRUE, upper=TRUE)
hc.t <- hclust(a2, method="ave")
par(mfrow=c(1,1))
png("cl2.png")
plot(hc.t,cex=0.35,xlab="metabolite",lty=1,lwd=1,col="dark blue",
     labels=risk.label)
dev.off()
