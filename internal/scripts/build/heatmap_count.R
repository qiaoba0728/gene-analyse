# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2021/3/17 0017

rm(list = ls())
library(pheatmap)
data = read.table("/data/output/expression_result/gene_count.csv",sep = "\t",header = T,row.names = 1)
data = as.data.frame(lapply(data,as.numeric))
data[is.na(data)] = 0
head(data)
cor_data = cor(data, method = "pearson")
png("/data/output/report_result/heatmap_count.png",width=800,height=800)
#write.table(cor_matr, file="cor_matr.xls",row.names=F, col.names=T,quote=FALSE,sep="\t")
pheatmap(cor_data,display_numbers = T,clustering_method = "average")
dev.off()

pdf("/data/output/report_result/heatmap_count.pdf")
#write.table(cor_matr, file="cor_matr.xls",row.names=F, col.names=T,quote=FALSE,sep="\t")
pheatmap(cor_data,display_numbers = T,clustering_method = "average")
dev.off()
