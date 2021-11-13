package scripts

const (
	DIFF_MATRIX = `rm(list = ls())
library(DESeq2)
args<-commandArgs(trailingOnly=TRUE)
geneFile<-args[1]
start<-args[2]
end<-args[3]
startRepeated <- args[4]
endRepeated <- args[5]
output<-args[6]
name<-args[7]
#geneFile<-"feature.all.counts.txt"
#start<-"6"
#end<-"9"
#repeated <- "3"
#output<-"./out/"
#name<-"test"

print(geneFile)
print(start)
print(end)
print(startRepeated)
print(endRepeated)
print(output)
print(name)
start<- as.numeric(start)
end<- as.numeric(end)
startRepeated<- as.numeric(startRepeated)
endRepeated<- as.numeric(endRepeated)
#输入必须是count格式
cf <- c(rep("sample1",startRepeated),rep("sample2",endRepeated))
counts <- read.table(geneFile, check.names = F, sep = "\t", row.names = 1, header = T)
m1 <- counts[c(start:(start+startRepeated - 1))]
head(m1)
m2 <- counts[c(end:(end+endRepeated - 1))]
head(m2)
count <- cbind(m1,m2)


count = as.data.frame(lapply(count,as.numeric))
count[is.na(count)] = 0
rownames(count) = rownames(counts)

head(count)
countCondition <- factor(cf)
countCondition
levels(countCondition)
colnames(count)
coldata <- data.frame(row.names=colnames(count), countCondition)
coldata
input <- data.matrix(count)
dds <- DESeqDataSetFromMatrix(countData=input, colData=coldata, design=~countCondition)
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

filterNameXls <- paste0(output,"/diffexpr-",name,"-0.05.csv")
#过滤后的文件
write.table(resdata_filter, file=filterName,sep="\t",row.names=F)
write.csv(resdata_filter, filterNameXls,row.names = T)


resName <- paste0(output,"/diffexpr-",name,".all.txt")
resName

resNameXls <- paste0(output,"/diffexpr-",name,".all.xls")
write.table(resdata, file=paste(resName,sep=""),quote=F,sep="\t",row.names=F) 
write.csv(resdata, resNameXls,row.names = T)
#MAplot
pdf(paste0(output,"/",name,"_MAplot.pdf"))
plotMA(res, main=paste0(output,"/",name,("_MAplot"), ylim=c(-10,10)))
dev.off()

png(paste0(output,"/",name,"_MAplot.png"),width=1200,height=600)
plotMA(res, main=paste0(output,"/",name,("_MAplot"), ylim=c(-10,10)))
dev.off()
#相关性
ddCor <- cor(count)

library(pheatmap)
heatmapName <- paste0(name,"_pheatmap.pdf")
pheatmap(file = paste0(output,"/",heatmapName), ddCor, clustering_method = "average", display_numbers = T)

heatmapName <- paste0(name,"_pheatmap.png")
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
pcaName <- paste0(output,"/",name,"_PCA.png")
png(pcaName,width=1200,height=600)
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

png(paste0(output,"/",name,"_volcano.png"),width=1200,height=600)
plot(no$log2FoldChange, no$padj, xlim = c(-10,10), ylim = c(0,100), col = "blue", pch = 16, cex = 0.8, main = "Gene Expression", xlab = "log2FoldChange", ylab = "-log10(Qvalue)")
points(up$log2FoldChange, up$padj, col = "red", pch = 16, cex = 0.8)
points(down$log2FoldChange, down$padj, col = "green", pch = 16, cex = 0.8)
dev.off()

s <- read.table(filterName, header = T, row.names = 1)
s <- na.omit(s)
s <- transform(s, padj = -1*log10(s$padj))
down <- s[s$log2FoldChange <= -1,]
up <- s[s$log2FoldChange >=1,]

data = data.frame(type=c("down","up","total"),value = c(dim(down)[1],dim(up)[1]),dim(down)[1] + dim(up)[1])
write.table (data,file =paste0(output,"/",name,"_down_up.csv"), row.names = FALSE, col.names =FALSE)

pdf(paste0(output,"/",name,"_down_up.pdf"))
ggplot(data=data,mapping=aes(x=type,y=value,fill=type))+
  geom_bar(stat="identity",width = 0.4) + geom_text(aes(label=(data$value)),vjust=-0.25) + 
  xlab(name)+ylab("Number of Genes")
dev.off()

png(paste0(output,"/",name,"_down_up.png"),width=1200,height=600)
ggplot(data=data,mapping=aes(x=type,y=value,fill=type))+
  geom_bar(stat="identity",width = 0.4) + geom_text(aes(label=(data$value)),vjust=-0.25) + 
  xlab(name)+ylab("Number of Genes")
dev.off()
`
)
