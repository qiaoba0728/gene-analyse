# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2022/2/4 0004

### 火山图
data <-read.table(file="/data/output/diff/volcano_plot.txt",header = TRUE, row.names =1,sep = "\t")
pdf("/data/output/diff/volcano_plot.pdf")
volcano <-ggplot(data = volcano_data,aes(x=log2FoldChange,y= -1*log10(padj)))+geom_point(aes(color=significant))+scale_color_manual(values = c("red","grey","blue")) + labs(title="Volcano_Plot",x=expression((log[2](FC)), y=expression(-log[10](padj)) ))+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)
volcano
dev.off()
png("/data/output/diff/volcano_plot.png",width=800,height=800)
volcano
dev.off()