# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2021/3/16 0016

rm(list = ls())
fpkm = read.table("/data/output/expression_result/gene_fpkm.csv",sep = "\t",header = T)
#读取基因表达值数据
#fpkm <- read.delim('gene_fpkm.csv', row.names = 1, sep = '\t')

#过滤标准:保留基因在大于floor(ncol(exprSet)/10)个样本中表达量大于1
#fpkm_filter=fpkm[apply(fpkm,1,function(x) sum(x>1) > floor(ncol(fpkm)/10)),]
dim(fpkm)
result = NULL
for (i in 2:ncol(fpkm)) {
    temp = data.frame(value = fpkm[,i],sample = colnames(fpkm)[i])
    result = rbind(result,temp)
}
result = as.data.frame(lapply(result,as.numeric))
result[is.na(result)] = 0
head(result)

result$value = log10(result$value)
library(ggplot2)
#result = data.matrix(log10(result+0.01))
pdf("/data/output/report_result/fpkm.pdf")
ggplot(result,aes(value,color=sample)) +
    xlab("log10(fpkm)") +
    geom_density(alpha = 0.3) +
    geom_rug() + theme_bw()
dev.off()

png("/data/output/report_result/fpkm.png",width=1200,height=600)
ggplot(result,aes(value, color=sample)) +
    xlab("log10(fpkm)") +
    geom_density(alpha = 0.3) +
    geom_rug() + theme_bw()
dev.off()
# 添加分面
library(ggpubr)

pdf("/data/output/report_result/fpkm_all.pdf")
ggdensity(result, x = "value",
facet.by = "sample", linetype = "sample",
rug = TRUE, xlab = "log10(fpkm)",
color = "sample", fill = "sample")
dev.off()
png("/data/output/report_result/fpkm_all.png",width=1200,height=600)
ggdensity(result, x = "value",
facet.by = "sample", linetype = "sample",
rug = TRUE, xlab = "log10(fpkm)",
color = "sample", fill = "sample")
dev.off()

png("/data/output/report_result/fpkm_violin.png")
ggplot(result, aes(sample,value))+
    geom_violin(aes(fill = sample),trim = FALSE)+
    geom_boxplot(width = 0.2)+
#scale_fill_manual(values=c(brewer.pal(5,"Set2")[c(1,3,2,5)]))+
    theme_classic()+
    labs(x='Stage',y='The expression level',title='Gene name')+
    theme(panel.background=element_rect(fill="white",colour="black",size=0.25),axis.line=element_line(colour="black",size=0.25),axis.title=element_text(size=5,face="plain",color="black"),axis.text = element_text(size=12,face="plain",color="black"),legend.position="none")
dev.off()

pdf("/data/output/report_result/fpkm_violin.pdf")
ggplot(result, aes(sample,value))+
    geom_violin(aes(fill = sample),trim = FALSE)+
    geom_boxplot(width = 0.2)+
#scale_fill_manual(values=c(brewer.pal(5,"Set2")[c(1,3,2,5)]))+
    theme_classic()+
    labs(x='Stage',y='The expression level',title='Gene name')+
    theme(panel.background=element_rect(fill="white",colour="black",size=0.25),axis.line=element_line(colour="black",size=0.25),axis.title=element_text(size=5,face="plain",color="black"),axis.text = element_text(size=12,face="plain",color="black"),legend.position="none")
dev.off()