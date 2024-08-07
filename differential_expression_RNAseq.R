##Load all the objects and libraries

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
load("Robjects/preprocessing.Rdata")

##Create design matrix
sampleinfo <- read.delim("data/SampleInfo_Corrected.txt")
design

##estimating the dispersion

dgeObj <- estimateCommonDisp(dgeObj)

##estimate genewise dispersion estimates
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)

#Plot dispersions
plotBCV(dgeObj)

##Testing for differential expression


# Fit the linear model
fit <- glmFit(dgeObj, design)
names(fit)
head(coef(fit))

#Conduct likelihood ratio tests for luminal vs basal and show the top genes:
lrt.BvsL <- glmLRT(fit, coef=2) topTags(lrt.BvsL)

