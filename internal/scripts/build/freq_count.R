# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2021/3/17 0017

rm(list = ls())
library(ggplot2)
data <- read.table("/data/output/expression_result/gene_count_number.csv",sep = "\t",header = T,row.names = 1,stringsAsFactors = F)
#data <- data.frame(lapply(data,as.numeric))
ncol(data)

for (i in 1:ncol(data)) {
    print(i)
    a <- table(cut(data[,i],breaks=c(0,10,100,1000,10000,1000000)))
    d <- as.data.frame(a)
    group <- c("0-10","10-100","100-1000","1000-10000",">=10000")
    d$group  <- group
    print(d)
    p <- ggplot(data=d,mapping=aes(x=group,y=Freq,fill=group,group=factor(1)))+geom_bar(stat="identity")
    filename <- paste0("/data/output/report_result/",colnames(data)[i],"_freq_count.png")
    print(filename)
    png(filename,width=800,height=800)
    print(p)
    dev.off()

    filename <- paste0("/data/output/report_result/",colnames(data)[i],"_freq_count.pdf")
    print(filename)
    pdf(filename)
    print(p)
    dev.off()
    print("finished")
}