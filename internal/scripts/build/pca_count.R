# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2021/3/17 0017

rm(list = ls())
library(gmodels)
library(ggplot2)
library(ggthemes)
library(ggpubr)
data = read.table("/data/output/expression_result/gene_count.csv",sep = "\t",header = T,row.names = 1)
data = as.data.frame(lapply(data,as.numeric))
data[is.na(data)] = 0
head(data)
pca.info <- fast.prcomp(data)
pca.data <- data.frame(sample = rownames(pca.info$rotation),Type = colnames(data),pca.info$rotation)
png("/data/output/report_result/pca_count.png",width=800,height=800)
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type") + theme_bw()
dev.off()
pdf("/data/output/report_result/pca_count.pdf")
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type") + theme_bw()
dev.off()