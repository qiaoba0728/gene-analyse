# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2021/3/17 0017

rm(list = ls())
data = read.table("/data/output/expression_result/gene_count.csv",sep = "\t",header = T,row.names = 1)
for (i in 1:ncol(data)) {
    temp = data.frame(value = data[,i],sample = colnames(data)[i])
    temp$value.type[temp$value < 5] = "0-5"
    temp$value.type[temp$value >= 5 & temp$value < 100] = "5-100"
    temp$value.type[temp$value >= 100 & temp$value < 500] = "100-500"
    temp$value.type[temp$value >= 500 & temp$value < 1000] = "500-1000"
    temp$value.type[temp$value >= 1000] = ">=1000"
    a = data.frame(table(temp$value.type))
    x = factor(a$Var1,levels = c("0-5","5-100","100-500","500-1000",">=1000"))
    filename = paste0("/data/output/report_result/",colnames(data)[i],"_freq_count.png")
    png(filename,width=800,height=800)
    barplot(Freq~x, data=a,main="gene count",xlab="gene expression level",ylab="count")
    dev.off()
    filename = paste0("/data/output/report_result/",colnames(data)[i],"_freq_count.pdf")
    pdf(filename)
    barplot(Freq~x, data=a,main="gene count",xlab="gene expression level",ylab="count")
    dev.off()
}