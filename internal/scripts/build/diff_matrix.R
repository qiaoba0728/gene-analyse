rm(list = ls())
library(DESeq2)
args<-commandArgs(trailingOnly=TRUE)
geneFile<-args[1]
start<-args[2]
end<-args[3]
repeated <- args[4]
output<-args[5]
name<-args[6]
#geneFile<-"feature.all.counts.txt"
#start<-"6"
#end<-"9"
#repeated <- "3"
#output<-"./out/"
#name<-"test"

print(geneFile)
print(start)
print(end)
print(repeated)
print(output)
print(name)
start<- as.numeric(start)
end<- as.numeric(end)
repeated<- as.numeric(repeated)
#输入必须是count格式
cf <- c(rep("sample1",repeated),rep("sample2",repeated))
counts <- read.table(geneFile, check.names = F, sep = "\t", row.names = 1, header = T)
m1 <- counts[c(start:(start+repeated - 1))]
m2 <- counts[c(end:(end+repeated - 1))]
count <- cbind(m1,m2)
#count <- data.matrix(count)


count = as.data.frame(lapply(count,as.numeric))
count[is.na(count)] = 0
#count <- data.matrix(count)
rownames(count) = rownames(counts)

head(count)
countCondition <- factor(cf)
countCondition
levels(countCondition)
colnames(count)
coldata <- data.frame(row.names=colnames(count), countCondition)
coldata
dds <- DESeqDataSetFromMatrix(countData=count, colData=coldata, design=~countCondition)
# 过滤掉没表达或者几乎不表达的
# 所有counts 加和 不足1 的 过滤掉
dds <- dds[rowSums(counts(dds)) > 1,]
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$padj), ]
resdata <- merge(as.data.frame(resOrdered), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
resdata_filter <- subset(resdata,resdata$padj<0.05 & abs(resdata$log2FoldChange) >1)
dim(resdata_filter)
head(resdata_filter)
filterName <- paste0(output,"/diffexpr-",name,"-0.05.txt")
filterName
#过滤后的文件
write.table(resdata_filter, file=filterName,sep="\t",row.names=F)
resName <- paste0(output,"/diffexpr-",name,".all.txt")
resName
write.table(resdata, file=paste(resName,sep=""),quote=F,sep="\t",row.names=F)

#MAplot
pdf(paste0(output,"/",name,"_MAplot.pdf"))
plotMA(res, main=paste0(output,"/",name,("_MAplot"), ylim=c(-10,10)))
dev.off()
#相关性
ddCor <- cor(count)
library(pheatmap)
heatmapName <- paste0(name,"_pheatmap.pdf")
pheatmap(file = paste0(output,"/",heatmapName), ddCor, clustering_method = "average", display_numbers = T)
colData(dds)
#PCA
vsd <- varianceStabilizingTransformation(dds)
library(ggplot2)
data <- plotPCA(vsd, intgroup=c("countCondition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pcaName <- paste0(output,"/",name,"_PCA.pdf")
pdf(pcaName)
ggplot(data, aes(PC1, PC2, color = colnames(dds))) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + labs(title = "Sample Vs counts PCA") + theme_bw()
dev.off()
#火山图
g <- read.table(resName, header = T, row.names = 1)
g <- na.omit(g)
g <- transform(g, padj = -1*log10(g$padj))
down <- g[g$log2FoldChange <= -1,]
up <- g[g$log2FoldChange >=1,]
no <- g[g$log2FoldChange > -1 & g$log2FoldChange < 1,]

pdf(paste0(output,"/",name,"_volcano.pdf"))
plot(no$log2FoldChange, no$padj, xlim = c(-10,10), ylim = c(0,100), col = "blue", pch = 16, cex = 0.8, main = "Gene Expression", xlab = "log2FoldChange", ylab = "-log10(Qvalue)")
points(up$log2FoldChange, up$padj, col = "red", pch = 16, cex = 0.8)
points(down$log2FoldChange, down$padj, col = "green", pch = 16, cex = 0.8)
dev.off()