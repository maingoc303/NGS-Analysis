df <- read.csv("Data2.csv",header = T)
df <- read.csv(choose.files(),header = T)
head(df)
# SUV_max is response
# 2 groups of variables:
## Z: clinical variables: BMI, Age, Gender, Diabetes, Smoking
## X: Metabolic variables

library(ggplot2)
library(plotly)
# - # - # - # - # - # - # - # - # - # - # - # - # - # - #
#           DATA EXPLORATION ANALYSIS                   #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - #
png("eda1.png")
qplot(Smoking,SUV_max,data=df, geom=c("boxplot", "jitter")
      ,fill = Smoking, xlab = "Diabetes",
      main="SUV_max by smoking and diabetes status")
dev.off()
png("eda2.png")
qplot(Age,SUV_max,data=df, geom=c("point", "smooth")
      ,method="lm", formula=y~x,color = Smoking,
      main="SUV max value by Age and smoking status")
dev.off()
png("eda3.png")
qplot(Age,SUV_max,data=df, geom=c("point", "smooth")
      ,method="lm", formula=y~x,color = Diabetes,
      main="SUV max value by Age and Diabetes status")
dev.off()
png("eda4.png")
qplot(BMI,SUV_max,data=df, geom=c("point", "smooth")
      ,method="lm", formula=y~x,color = Smoking,
      main="SUV max value by BMI and Smoking status")
dev.off()
qplot(BMI,SUV_max,data=df, geom=c("point", "smooth")
      ,method="lm", formula=y~x,color = Diabetes,
      main="SUV max value by BMI and Diabetes status")
dev.off()

# - # - # - # - # - # - # - # - # - # - # - # - # - # - #
#           CLINICAL VARIABLE ANALYSIS                  #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - #
m1 <- lm(df$SUV_max~df$BMI+df$Age+as.factor(df$Gender)+
           as.factor(df$Diabetes)+as.factor(df$Smoking))

# - # - # - # - # - # - # - # - # - # - # - # - # - # - #
#           METABOLIC VARIABLE ANALYSIS                 #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - #

png("eda6.png")
plot(cor(df[,c(2,8:117)])[1,c(2:111)],pch=20,col=4,
     main="Correlation between SUV_max and metabolic variables",
     xlab="Metabolic variables",ylab="Correlation value",
     ylim=c(-1,1))
dev.off()

m2 <- lm(df$SUV_max~as.matrix(df[,c(8:117)])+as.factor(df$Smoking))
summary(m2)

# predicting with simple linear regression
m3 <- lm(df$SUV_max~as.matrix(df[,c(8:117)])+as.factor(df$Smoking))
pca <- princomp(df[,c(2,8:117)],scores=TRUE,cor = TRUE)
plot(pca, type = "l", lwd = 2, col="blue",
     main="Principle Component Analysis")

length(m3$coefficients)

png("eda7.png")
par(mfrow=c(1,2))
plot((df$SUV_max),pch=20,col=3,
     xlab="Patient",ylab="SUV max",
     main="SUV max by patient")
boxplot(log2(df$SUV_max),col="light blue",
        ylab="Log2 of SUVmax", main="log2 of SUVmax value")
dev.off()

# splitting data to training set and vadidation (test) set
# using permutation instead of bootstrap

#######
# PCA to see what component contribute for variance explanation
########



n <- dim(df)[1]

df_new <- df[sample(n),]
train.index <- 1:round(0.7*n)
train.df <- df[train.index,]
test.index <- round(.7*n):n
test.df <- df[test.index,]


# Q1: Develop a predictive model using SPCA - Supervised Principle Component Analysis

### original data

summary(df[,c(2:4)])
data.frame(table(df[,c(5,6,7)]))

clinical_fit <- lm(df$SUV_max~df$BMI+df$Age+df$Gender+df$Diabetes+df$Smoking)

summary(clinical_fit)

# ------ Model Based Analysis ------ #
df$risk<-ifelse(df$SUV_max<=4,0,1)
df$risk<-as.factor(df$risk)
data1<-t(df[,c(8:118)])
t.t<-beta.i<-p.val<-rep(0,110)
# GLM #
for(i in 1:110){	
  xi<-as.numeric(data1[i,])
  modi<-glm(as.factor(data1[111,])~xi,binomial(link = "logit"))
  t.t[i]<-summary(modi)$coefficients[2,3]
  beta.i[i]<-summary(modi)$coefficients[2,1]
  p.val[i]<-summary(modi)$coefficients[2,4]
}

png("pval.png")
hist(p.val,nclass=50,col="light green", 
     main = "P-value of metabolic parameter")
dev.off()
plot(p.val)#,nclass=50)

library(Biobase)

library(multtest, verbose = FALSE)

rawp<-p.val
holm<-mt.rawp2adjp(rawp, proc=c("Holm"))
bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
bh<-mt.rawp2adjp(rawp, proc=c("BH"))
by<-mt.rawp2adjp(rawp, proc=c("BY"))
allp<-cbind(rawp, bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2])

png("multi.png")
mt.plot(allp,plottype="pvsr", 
        proc=c("rawp","Bonferroni","Holm","BH","BY"),
        leg=c(2,0.95),lty=1,col=1:5,lwd=2,
        main="Multiplicity Adjustment for Metabolic Signatures")
abline(h=.1,lty=2)
abline(h=.05,lty=2)
text(80,.12,"p-val = 10%",cex=.7,col=2)
text(80,.04,"p-val = 5%",cex=.7,col=2)
dev.off()
mt.reject(cbind(rawp,bonf$adjp[order(bonf$index),2],holm$adjp[order(holm$index),2],bh$adjp[order(bh$index),2],by$adjp[order(by$index),2]),seq(0,1,0.05))$r

# limma
meta.exp<-as.matrix(t(df[,c(8:117)]))
nn<-dim(meta.exp)[1]
nn
library(limma)

design <- ((df$risk))#cbind(rep(1,38),abs(classlabel-1))
design

fit <- lmFit(meta.exp,design)
summary(fit)

fit1 <- eBayes(fit)
head(fit1$t)
head(fit1$p.value)

topTable <- topTable(fit1,number=nrow(meta.exp))#,coef=2)

head(topTable)

plot(t.t,-log(p.val))
points(topTable[,2],-log(topTable[,4]),col=2)
plot(topTable[,2],-log(topTable[,4]))

############################################################################################
#### metabolite by metabolite regression                                                ####
############################################################################################

data1<-df
meta<-dim(data1[,c(8:117)])[2]
sd.i<-e.i<-t.t<-p.val<-c(1:meta)
t.tw<-p.valw<-c(1:meta)

## t test ##

for(i in 1:meta){	
  #cat(i)
  x1<-data1[i,]
  x2<-data1[i,]
  t.i<-t.test(x1,x2, alternative="two.sided")
  t.w<-wilcox.test(x1,x2)
  t.t[i]<-t.i$statistic
  p.val[i]<-t.i$p.value
  t.tw[i]<-t.w$statistic
  p.valw[i]<-t.w$p.value
  e.i[i]<-t.i$estimate[1]-t.i$estimate[2]
  sd.i[i]<-e.i[i]/t.t[i]
}

p.val.beta1<-t.val.beta1<-c(1:110)


for(i in 1:110){	
  #cat(i)
  fit.lm<-lm(df$SUV_max~as.numeric(df[,i+7]))#+as.factor(gr))
  p.val.beta1[i]<-summary(fit.lm)$coeff[2,4]
  t.val.beta1[i]<-summary(fit.lm)$coeff[2,3]
}

risk <- df$SUV_max
meta.sig<-df[,c(8:117)]
fit <- lmFit(t(meta.sig),risk)
summary(fit)
fit1 <- eBayes(fit)
topTable <- topTable(fit1,number=nrow(t(meta.sig)),coef=1)
head(topTable)
plot(t.val.beta1,-log(p.val.beta1))
points(topTable[,2],-log(topTable[,1]),col=2)
head(topTable)

png("eda5.png")
plot(p.val.beta1,pch=20,col="orange",xlab="Metabolic variable",
     ylab="p-value",
     main="Fitting linear regression by each metabolic variable")
dev.off()


# ------ SPCA top metabolites ------ #
K<-11
cutoff<-sort(p.val)[K]
sum(p.val < cutoff)
min(p.val)

index<-c(1:110)
index1<-index[p.val < cutoff]
index1

# ------ reduced data set with top metabolites ------ #

df.R<-df[,c(2,index1+7)]
dim(df.R)

# ------ PCA for full data set ------ #
pc.cr.f <- princomp(df[,c(2,8:117)],scores=TRUE,cor = TRUE)
names(pc.cr.f)
x<-prcomp(df[,c(8:117)])
names(x)

loadings(pc.cr.f)
png("pca1.png")
plot(pc.cr.f,main="PCA with all metabolic signature",col="light blue")
plot(pc.cr.f$sdev^2/sum(pc.cr.f$sdev^2),type="p",pch=20,
       col=4,cex=1,ylab="proportions of variance",xlab="Component",
     main="Proportions of variance by component")
lines(pc.cr.f$sdev^2/sum(pc.cr.f$sdev^2),lty=2,col="dark blue")
points(pc.cr.f$sdev^2/sum(pc.cr.f$sdev^2),ltw=2,
       col=3,pch=20)
dev.off()
biplot(pc.cr.f)
x<-sort((pc.cr.f$sdev)^2)

var.exp <- sort(pc.cr.f$loadings,decreasing = T)

# ------ 1st PCA ------ #
df$risk <- ifelse(df$SUV_max<=9,0,1)

library(ggplot2)
first.pca<-as.vector(pca$scores[,1])
sec.pca<-as.vector(pca$scores[,2])
qplot(x=first.pca,y=df$SUV_max,color=as.factor(df$risk))
qplot(x=sec.pca,y=df$SUV_max,color=as.factor(df$risk))
plot(first.pca,df$SUV_max)
summary(lm(df$SUV_max~first.pca))

# ------ PCA for reduced data set ------ #
png("pca2.png")
pca.r <- princomp(df.R,scores=TRUE,cor = TRUE)
plot(pca.r$sdev^2/sum(pca.r$sdev^2),type="p",pch=20,
     col=4,cex=1,main="Proportions of variance by component",
     ylab="proportions of variance")
lines(pca.r$sdev^2/sum(pca.r$sdev^2),lty=20,col="dark blue")
dev.off()

loadings(pca.r)
plot(pca.r,main="PCA with top 10 metabolic signatures")
biplot(pca.r)
summary(pc.cr.r)

install.packages("mdatools")
library(mdatools)

model = pca(df.R[,c(2:11)], scale = TRUE, cv = 10, info = 'Simple PCA model')
model = selectCompNum(model, 1)

#
res = model$calres
summary(res)

png("cv.pca.png")
plot(res)
dev.off()
png("scorepca.png")
plotScores(res, cgroup = df[, index1[3]], show.labels = TRUE)
dev.off()
#

summary(model)
png("cvpca.png")
plot(model, show.labels = TRUE)
dev.off()
## 3. Show scores and loadings plots for the model
par(mfrow = c(2, 2))
plotScores(model, comp = c(1, 3), show.labels = TRUE)
plotScores(model, comp = 2, type = 'h', show.labels = TRUE)
plotLoadings(model, comp = c(1, 3), show.labels = TRUE)
plotLoadings(model, comp = c(1, 2), type = 'h', show.labels = TRUE)
par(mfrow = c(1, 1))

## 4. Show residuals and variance plots for the model
par(mfrow = c(2, 2))
plotVariance(model, type = 'h')
plotCumVariance(model, show.labels = TRUE, legend.position = 'bottomright')
plotResiduals(model, show.labels = TRUE)
plotResiduals(model, ncomp = 2, show.labels = TRUE)
par(mfrow = c(1, 1))

# }

# model


# ------ LASSO ------ #
# install.packages("glmnet")
library(glmnet)

m1 <- lm(df$SUV_max~df$BMI+df$Age+as.factor(df$Gender)
         +as.factor(df$Diabetes)+as.factor(df$Smoking))

summary(m1)

set.seed(1)
y<-df$SUV_max
x <- as.matrix(df[,c(8:117)])
lasso.fit1 <- glmnet(x,y,alpha = 1,family = "gaussian")
cv.lasso.fit1 <- cv.glmnet(x,y,alpha=1)

png("lamda1.png")
plot(cv.lasso.fit1)
dev.off()
plot(lasso.fit1, xvar="lambda")
abline(v=cv.lasso.fit1$lambda.min,col="blue",lwd=3)
cv.lasso.fit1$lambda.min
log(cv.lasso.fit1$lambda.min)

png("lamda2.png")
plot(lasso.fit1,xvar="lambda") 
abline(v=log(cv.lasso.fit1$lambda.min),col="blue",lwd=3)
dev.off()

tmp_coeffs <- coef(cv.lasso.fit1, s = "lambda.min")
coeff <- as.data.frame(as.matrix(tmp_coeffs))
sum(coeff!=0)

summary(df)

# LASSO with top signatures
set.seed(2)
y<-df.R$SUV_max
x <- as.matrix(df.R[,c(2:11)])
lasso.fit2 <- glmnet(x,y,alpha = 1,family = "gaussian")
cv.lasso.fit2 <- cv.glmnet(x,y,alpha=1)

png("lamda3.png")
plot(cv.lasso.fit2)
dev.off()

plot(lasso.fit1, xvar="lambda")
abline(v=cv.lasso.fit1$lambda.min,col="blue",lwd=3)
cv.lasso.fit1$lambda.min
log(cv.lasso.fit1$lambda.min)

png("lamda4.png")
plot(lasso.fit2,xvar="lambda") 
abline(v=log(cv.lasso.fit2$lambda.min),col="blue",lwd=3)
dev.off()
cv.lasso.fit2$lambda.min

df.R$Smoke<-df$Smoking
n<-dim(df.R)[1]
id<-sample(n,n*.7)
df.R.t<-df.R[id,]
df.R.v<-df.R[-id,]
y<-df.R.t$SUV_max
x <- as.matrix(df.R.t[,c(2:11)])
lasso1 <- glmnet(x,y,alpha = 0.1986317,family = "gaussian")

pred1<-predict.glmnet(lasso1,df.R.v)

######## Elastic Net with alpha = 0.5 ######
set.seed(2)
y<-df$SUV_max
x <- as.matrix(df[,c(8:117)])
Elnetfit1 <- glmnet(x,y,alpha=.7, family ="gaussian")
cv.Elnetfit1<-cv.glmnet(x,y,alpha=.7,nfolds = 10)
plot(Elnetfit1,xvar="lambda")
abline(v=log(cv.Elnetfit1$lambda.min),col="blue",lwd=3)

png("elas1.png")
plot(cv.Elnetfit1)
dev.off()
png("elas2.png")
plot(Elnetfit1,xvar="lambda")
abline(v=log(cv.Elnetfit1$lambda.min),col="blue",lwd=3)
dev.off()

cv.Elnetfit1.bin <- cv.glmnet(x,y,alpha=0.5, family ="binomial",nfolds = 6)
plot(cv.Elnetfit1.bin)

plot(Elnetfit1,xvar="lambda")
abline(v=log(cv.Elnetfit1$lambda.min),col="blue",lwd=3)

# reduced data
set.seed(4)
y<-df.R$SUV_max
x <- as.matrix(df.R[,c(2:11)])
Elnetfit2 <- glmnet(x,y,alpha=.5, family ="gaussian")
cv.Elnetfit2<-cv.glmnet(x,y,alpha=.5,nfolds = 10)
plot(Elnetfit2,xvar="lambda")
abline(v=log(cv.Elnetfit2$lambda.min),col="blue",lwd=3)

png("elas3.png")
plot(cv.Elnetfit2)
dev.off()
png("elas4.png")
plot(Elnetfit1,xvar="lambda")
abline(v=log(cv.Elnetfit1$lambda.min),col="blue",lwd=3)
dev.off()
cv.Elnetfit1$lambda.min


# Apply LIMMA
install.packages("limma")
library(limma)
df1<-df[,c(2,8:117)]
fit1 <- lmFit(df1$SUV_max~as.matrix(df1[,c(2:111)]))
parm <- summary(fit1)
p.val <- sort(parm$coefficients[,4])
names(p.val[1:10])


########################################################################
# dummy variable for Smoking
names(df)

n <- dim(df)[1]

df_new <- df[sample(n),]
train.index <- 1:round(0.7*n)
train.df <- df[train.index,]
test.index <- round(.7*n):n
test.df <- df[test.index,]

##########################################################################################

##########################################################################################

##########################################################################################
# dummy variables for Smoking
smoking <- as.numeric(df$Smoking=='Yes')
former <- as.numeric(df$Smoking=='Former')


# correlation matrix for top 20 significant metabolites from reduced data
corMatrix <- df.R[,c(2:21)]
# Create the correlation matrix
M <- round(cor(corMatrix), 2)
# Create corrplot
install.packages("corrplot")
library(corrplot)
corrplot(M, diag = FALSE, method="color", order="FPC", tl.srt = 90)

#divide the new data
pca.train <- df[1:n*.7,c(2,8:117)]
pca.test <- df[-(1:n*.7),c(2,8:117)]

pca.train.r <- df.R[1:n*.7,]
pca.test.r <- df.R[-(1:n*.7),]

#principal component analysis
prin_comp <- prcomp(pca.train[,c(2:111)], scale. = T, scores=T)
summary(prin_comp)
names(prin_comp)

# reduced to significant metabolites
prin_comp.r <- prcomp(pca.train.r[,c(2:21)], scale. = T)
summary(prin_comp.r)

#outputs the mean of variables
prin_comp$center
#outputs the standard deviation of variables
prin_comp$scale
#pca loading
prin_comp$rotation
# 4 principle components of first 5 variables (metabolites)
prin_comp$rotation[1:5,1:4]
# the matrix x for the principal component score: 110 dimensions (metabolites)
dim(prin_comp$x)
# plot the resultant principal components
biplot(prin_comp, scale = 0)
biplot(prin_comp.r, scale = 0)

#compute standard deviation of each principal component
std_dev <- prin_comp$sdev
std_dev.r <- prin_comp.r$sdev
#compute variance
pr_var <- std_dev^2
pr_var.r <- std_dev^2
#check variance of first 10 components
pr_var[1:10]
#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:10]
prop_varex.r <- pr_var.r/sum(pr_var.r)
#scree plot
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

plot(prop_varex.r, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")
#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")

plot(cumsum(prop_varex.r), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")


#############################################################


################################################################
# PCA using correlation matrix:

pr.comp <- prcomp(pca.train[,c(2:110)], scale = TRUE, center = TRUE)
summary(pr.comp)

plot(pr.comp$x[, c(1, 2)],col = (SUV + 1),
     xlab = "PC1", ylab = "PC2")

#add a training set with principal components
train.data <- data.frame(SUV_max=pca.train$SUV_max, prin_comp$x)

# interest in first 20 PCAs
train.data <- train.data[,1:21]
train.data<-as.data.frame(train.data)

library(MASS)
# run model for new train data set
model1 <- lm(log(SUV_max)~., data = train.data)

model1 <- lda(SUV_max~pca1, data = train.t)
a<-predict(model1,newdata = test.n)

#transform test into PCA
test.data <- predict(prin_comp, newdata = pca.test)
test.data <- as.data.frame(test.data)

#select the first 20 components
test.data <- test.data[,1:20]
#make prediction on test data
prediction <- predict(model1, newdata = test.data)
prediction

install.packages("ROCR")
library(ROCR)

PRMSS <- mean(pred.fit - pca.test$SUV_max)^2

summary(prediction)

summary(model1)

#run a decision tree
install.packages("rpart")
library(rpart)
rpart.model <- rpart(SUV_max ~ .,data = train.data, method = "anova")
rpart.model

########################################################################
# PCA for important (significant metabolites)

# performance evaluation for prediction model based on PCA
PCA_Estimate=predict(Step_PCA_Reg,type='response',
                     newdata=cbind(Clean_Data[test,c(Target,categoric)],
                                   pca1$ind$coord[test,]))
format(cor(prediction, pca.test)^2, digits=4)


############## modeling
Ux<-as.vector(pca.r$scores[,1])
y<-df$risk+1
risk<-ifelse(y<=4,0,1)
summary(glm(risk~Ux+as.factor(df$Smoking),family = binomial))

#################

library(samr)

d=list(
  x=t(df[,c(8:117)]),y=y,
  metaid=as.character(1:dim(df)[2]),
  metanames=names(df[,c(8:117)])#paste("gene",
                  #as.character(1:dim(data)[1]))
  ,logged2=FALSE)

samr.obj <- samr(d,resp.type="Two class unpaired")
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
samr.plot(samr.obj,1.5)
sigmeta.table<-
  samr.compute.siggenes.table(samr.obj,del=1.5,
                              data=d, delta.table)

# loop to choose top significant
n <- dim(df)[1]



x<-c(1:10)
t.t<-p.val<-rep(0,110)
index<-c(1:110)
K<-11
set.seed(111)
for (i in 1:1000){
  df.n <- df[sample(n),]
  y=df.n$SUV_max
  for(j in 1:110){	
    xi<-as.numeric(df[,j+7])
    modi<-lm(y~xi)#,binomial(link = "logit"))
    t.t[j]<-summary(modi)$coefficients[2,3]
    #beta.j[j]<-summary(modi)$coefficients[2,1]
    p.val[j]<-summary(modi)$coefficients[2,4]
  }
  
  rawp<-p.val
  #holm<-mt.rawp2adjp(rawp, proc=c("Holm"))
  #bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
  bh<-mt.rawp2adjp(rawp, proc=c("BH"))
  cutoff<-sort(p.val)[K]
  #name <- paste("r",i)
  #top.meta$'i'=index[p.val < cutoff]
  x<-cbind(x,index[p.val < cutoff])
}

freq.meta<-c(1:110)
for(i in 1:110){
  freq.meta[i]=sum(x[,c(2:1001)]==i)
}

which(freq.meta>115)+7

png("topsig.png")
barplot(freq.meta[freq.meta>115],
        main="Top 10 metabolic signatures",
        cex.names = which(freq.meta>115))
dev.off()

########### top 8 signatures #############
top.sig<-which(freq.meta>120)+7
data1<-df[,c(2,7,top.sig)]
set.seed(123)  
n<-dim(data1)[1]
size <- floor(0.5*n)
MSE.test <- c(1:1000)
MSE.train <- c(1:1000)
Rsquare <- c(1:1000)
overfit<-c(1:1000)
first.var<-c(1:1000)
test.var<-c(1:1000)
for (i in 1:1000){
  # splitting data
  id<-sample(n,size)
  train <- data1[id,]
  test <- data1[-id,]
  pca<-princomp(train[,c(3:13)])
  first.pca<-as.vector(pca$scores[,1])
  first.var[i]<-(pca$sdev[1])^2/sum(pca$sdev^2)
  #sec.pca<-as.vector(pca$scores[,2])
  train.t<-train[,c(1:2)]
  train.t$pca1<-first.pca
  fit<-lm(train.t$SUV_max~train.t$Smoking+train.t$pca1)#train.t$Smoking+pca1)
  MSE.train[i]<-sum(fit$residuals^2)/dim(train.t)[1]
  pca.t<-princomp(test[,c(3:13)])
  fst.com<-as.vector(pca.t$scores[,1])
  #coef<-fit$coefficients
  test.n<-test[,c(1:2)]
  test.n$pca1<-fst.com
  pred<-predict.lm(fit,newdata = test.n)#as.factor(test$Smoking)*
  test$pre<-pred
  MSE.test[i]<-sum((test.n$SUV_max-pred)^2)/dim(test.n)[1]
  #Rsquare[i]<-cor(pred,test.n$SUV_max)^2
  overfit[i]<-MSE.test[i]-MSE.train[i]
  test.var[i]<-(pca.t$sdev[1])^2/sum(pca.t$sdev^2)
}

plot(MSE.test)
points(MSE.train,col=2)
hist(overfit)
hist(Rsquare)

hist(MSE.test,breaks = 25)
hist(MSE.train,breaks = 25)

plot(test$SUV_max,col=2)
points(test$pre,col=test$Smoking)

plot(first.var)
plot(test.var)
# glm
set.seed(456)
df.n<-df
df.n$risk<-as.factor(ifelse(df.n$SUV_max<=4,0,1))
n<-dim(df.n)[1]
MSE.test.n <- c(1:1000)
MSE.train.n <- c(1:1000)
Rsquare.n <- c(1:1000)
overfit.n<-c(1:1000)
first.var.n<-c(1:1000)
for (i in 1:1000){
  # splitting data
  train <- df.n[sample(n,n*.5,replace = T),c(2,7,top.sig,118)]
  test <- df.n[-sample(n,n*.5,replace = T),c(2,7,top.sig,118)]
  pca<-princomp(train[,c(3:13)])
  first.pca.n<-as.vector(pca$scores[,1])
  first.var.n[i]<-(pca$sdev[1])^2/sum(pca$sdev^2)
  #sec.pca<-as.vector(pca$scores[,2])
  fit<-glm(train$risk~as.factor(train$Smoking)+first.pca.n,family = "binomial")
  MSE.train.n[i]<-sum(fit$residuals^2)/dim(train)[1]
  pred<-predict.glm(fit,test)
  MSE.test.n[i]<-sum((test$SUV_max-pred)^2)/dim(test)[1]
  #Rsquare.n[i]<-cor(pred,test$SUV_max)^2
  overfit.n[i]<-MSE.test[i]-MSE.train[i]
}
plot(exp(pred))
points(test$SUV_max,col=2)



###################### top 8
top.sig<-which(freq.meta>120)+7
top3<-which(freq.meta>125)+7

data8<-df[,c(2,7,top.sig)]
data3<-df[,c(2,7,top3)]

model8 = pca(df[,c(3:10)], scale = TRUE, cv = 5, info = 'Simple PCA model')
model = selectCompNum(model, 1)

#
res = model8$calres
summary(res)
plot(res)
plotScores(res, cgroup = df[, index1[3]], show.labels = TRUE)
plot(model8, show.labels = TRUE)

set.seed(123)  
n<-dim(data)[1]
size <- floor(0.7*n)
MSE.test <- c(1:1000)
MSE.train <- c(1:1000)
Rsquare <- c(1:1000)
overfit<-c(1:1000)
first.var<-c(1:1000)
for (i in 1:1000){
  # splitting data
  id<-sample(n,size)
  train <- data[id,]
  test <- data[-id,]
  pca<-princomp(train[,c(3:13)])
  first.pca<-as.vector(pca$scores[,1])
  first.var[i]<-(pca$sdev[1])^2/sum(pca$sdev^2)
  #sec.pca<-as.vector(pca$scores[,2])
  train.t<-train[,c(1:2)]
  train.t$pca1<-first.pca
  fit<-lm(SUV_max~.,data=train.t)#train.t$Smoking+pca1)
  MSE.train[i]<-sum(fit$residuals^2)/dim(train.t)[1]
  pca.t<-princomp(test[,c(3:13)])
  fst.com<-as.vector(pca.t$scores[,1])
  #coef<-fit$coefficients
  test.n<-test[,c(1:2)]
  test.n$pca1<-fst.com
  pred<-predict.lm(fit,newdata = test.n)#as.factor(test$Smoking)*
  MSE.test[i]<-sum((test.n$SUV_max-pred)^2)/dim(test.n)[1]
  #Rsquare[i]<-cor(pred,test.n$SUV_max)^2
  overfit[i]<-MSE.test[i]-MSE.train[i]
}



# run model for new train data set

model1 <- lda(SUV_max~pca1, data = train.t)
a<-predict(model1,newdata = test.n)

head(a$posterior)
plot(a$x)

#transform test into PCA
test.data <- predict(prin_comp, newdata = pca.test)
test.data <- as.data.frame(test.data)

#select the first 20 components
test.data <- test.data[,1:20]
#make prediction on test data
top10<-which(freq.meta>117)+7
data1<-df[,c(2,7,top10)]
id<-sample(n,n*.5)
train <- df.R[id,]
test <- df.R[-id,]
pca1<-princomp(train[,c(2:11)])
pca.t<-princomp(test[,c(2:11)])
com1<-as.vector(pca.t$scores[,1])
test$comp<-com1
test.set<-test[,c(1,12:13)]
model1<-lm(train$SUV_max~train$Smoke+as.vector(pca1$scores[,1]))
pred1<-predict(model1,test.set)
sum((pred1-test.set$SUV_max)^2)
prediction <- predict(model1, newdata = test.set)
prediction

png("overfit.png")
plot(test.set$SUV_max,xlab="Observation",ylab="SUV max",
     main="Overfitting Checking with SPCA")
points(prediction,col=2)
legend(0,40,col=c(1,2),pch=c(1,1),legend = c("Observed","Predict"))
plot(test.set$SUV_max,prediction)
cor(test.set$SUV_max,prediction)
dev.off()

dim(df.R)
names(df.R)

install.packages("rggobi")
library(rggobi)
g <- ggobi(df.s)


fit.n <- lda(risk.factor ~ meta31+meta45+meta49+meta50
           +meta57+meta91, data=data.frame(df.s), CV=TRUE)
ct <- table(risk.label, fit.n$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))
fit.q <- qda(risk.factor ~ meta31+meta45+meta49+meta50
             +meta57+meta91, data=data.frame(df.s))#, prior=c(1,1,1,1)/3)
plot(fit.n)
plot(fit.n, dimen=1, type="both")

install.packages("klaR")
library(klaR)

partimat(risk.factor ~ meta31+meta45+meta49+meta50
         +meta57+meta91, data=data.frame(df.s),method="lda")

pairs(df.s[,c(1:6)], main="My Title ", pch=20, 
      bg=c("red", "blue","yellow")[unclass(risk.factor)])

install.packages("FactoMineR")
library(FactoMineR)
result <- PCA(df.s)

d <- dist(df.s) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS")
text(x, y, labels = row.names(df.s), cex=.7)

d <- dist(df.s) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Nonmetric MDS",col=risk.label)

text(x, y, labels = row.names(df.s), cex=.7)


# Model Based Clustering
install.packages("mclust")
library(mclust)
fit <- Mclust(df.s)
plot(fit) # plot results 
summary(fit) # display the best model

mydata <- scale(df[,c(3:101)])

# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 3) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram
groups <- cutree(fit, k=2) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=2, border="red")

# Ward Hierarchical Clustering with Bootstrapped p values
install.packages("pvclust")
library(pvclust)

fit <- pvclust(mydata, method.hclust="ward.D2",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 2)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
install.packages("fpc")
library(fpc)
plotcluster(mydata, fit$cluster)
cluster.stats(d, fit1$cluster, fit2$cluster)



Y
s <- svd(Y)
What <- t(s$v[,c(1:2)])
colnames(What)<-colnames(Y)
round(What,2)

e<-scale(df[,c(3:101)])

ind <- which(df$Stage != "early")
y <- df$Stage[ind]
X <- df[ind,c(3:101)] 
dim(X)
length(y)
library(class)
pred <- knn(train =  X, test = X, cl=y, k=2)
mean(y != pred)

library(caret)
set.seed(1)
idx <- createFolds(y, k=10)
sapply(idx, length)
y[idx[[1]]]

head( X[idx[[1]], which(p.val<.05)] )
sapply(idx, function(i) table(y[i]))

library(rafalib)
Xsmall <- cmdscale(dist(X))

plot(Xsmall,col=as.fumeric(c("advanced","early")))
legend("topleft",levels(factor(c("advanced","early"))),fill=seq_along(levels(factor(y))))

pred <- knn(train=Xsmall[ -idx[[1]] , ], test=Xsmall[ idx[[1]], ], cl=y[ -idx[[1]] ], k=5)
table(true=y[ idx[[1]] ], pred)
mean(y[ idx[[1]] ] != pred)
for (i in 1:10) {
  pred <- knn(train=Xsmall[ -idx[[i]] , ], test=Xsmall[ idx[[i]], ], cl=y[ -idx[[i]] ], k=5)
  print(paste0(i,") error rate: ", round(mean(y[ idx[[i]] ] != pred),3)))
}

