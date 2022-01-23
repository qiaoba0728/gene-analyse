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
plot(no$log2FoldChange, no$padj, xlim = c(-10,10), ylim = c(0,100), col = "grey", pch = 16, cex = 0.8, main = "Gene Expression", xlab = "log2FoldChange", ylab = "-log10(Qvalue)")
points(up$log2FoldChange, up$padj, col = "green", pch = 16, cex = 0.8)
points(down$log2FoldChange, down$padj, col = "red", pch = 16, cex = 0.8)
dev.off()

png(paste0(output,"/",name,"_volcano.png"),width=1200,height=600)
plot(no$log2FoldChange, no$padj, xlim = c(-10,10), ylim = c(0,100), col = "grey", pch = 16, cex = 0.8, main = "Gene Expression", xlab = "log2FoldChange", ylab = "-log10(Qvalue)")
points(up$log2FoldChange, up$padj, col = "green", pch = 16, cex = 0.8)
points(down$log2FoldChange, down$padj, col = "red", pch = 16, cex = 0.8)
dev.off()

s <- read.table(filterName, header = T, row.names = 1)
s <- na.omit(s)
s <- transform(s, padj = -1*log10(s$padj))
down <- s[s$log2FoldChange <= -1,]
up <- s[s$log2FoldChange >=1,]
dim(down)[1]
dim(up)[1]
data = data.frame(type=c("up","down"),value = c(dim(up)[1],dim(down)[1]))
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
	DIFF_VEEN = `args=commandArgs(T)
print(args[1])
install.packages("VennDiagram")
library (VennDiagram)
files = strsplit(args[1], ",")[[1]]
print(files)
print(length(files))

if (length(files) != 3)
{
        return(message("must input 3 param!"))
}
#for (f in files)
#{
#        #data <- fun.read(f)
#       #data
#       path <- paste("/data/output/diff","/",f,".txt",sep="")
#        print(path)
#        res <- read.table(path, header = TRUE, stringsAsFactors = FALSE)$Gene
#       head(res)
#        datas <- list(datas,res)
#}
path <- paste("/data/output/diff/diffexpr","-",files[1],"-0.05.txt",sep="")
A = read.table(path, header = TRUE, stringsAsFactors = FALSE)$Gene
path <- paste("/data/output/diff/diffexpr","-",files[2],"-0.05.txt",sep="")
B = read.table(path, header = TRUE, stringsAsFactors = FALSE)$Gene
path <- paste("/data/output/diff/diffexpr","-",files[3],"-0.05.txt",sep="")
C = read.table(path, header = TRUE, stringsAsFactors = FALSE)$Gene

result = list(A = A,B = B,C = C)
names(result)[1] = files[1]
names(result)[2] = files[2]
names(result)[3] = files[3]
#head(datas[1])

venn.diagram(x=result, "/data/output/report_result/veen.png", height = 450, width = 450,  resolution =300, imagetype="png",  col="transparent",fill=c("cornflowerblue","green","yellow"), alpha = 0.50, cex=0.45, cat.cex=0.45)
`
	MFUZZ = `
#安装 Mfuzz包；
BiocManager::install("Mfuzz")
#载入Mfuzz包；
library(Mfuzz)
# 2:4:xxx,5:7:yyy
args<-commandArgs(trailingOnly=TRUE)
groups<-args[1]
#设置工作目录；
#读入测试数据；
#输入数据一般为差异基因或蛋白的表达矩阵；
#推荐使用归一化（Normalisation）后的数值，比如FPKM值；
df <- read.csv("/data/output/expression_result/gene_count_number.csv",header = T,sep='\t')
head(df)

df1=data.frame(gene=df$gene)
#names = strsplit("xxx,yyy", ",")[[1]]
groups = strsplit(groups, ",")[[1]]
print(groups)
for(group in groups)
{
        indexs=strsplit(group, ":")[[1]]
        start=as.integer(indexs[1])
        end=as.integer(indexs[2])
        print(indexs[1])
        print(indexs[2])
        print(indexs[3])
        num=apply(df[,start:end], 1, mean)
        df1[indexs[3]]=num
}
head(df1)
row.names(df1) <- df1$gene
df1=df1[,-1]
head(df1)
#转成矩阵；
mat <- as.matrix(df1)
#查看前6行；
head(mat)
#创建用于Mfuzz的对象；
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
#Features Samples
#9345 4
#过滤缺失值超过25%的基因;
dt.r <- filter.NA(dt, thres=0.25)
#以均值的方式填充缺失值；
dt.f <- fill.NA(dt.r,mode="mean")
#过滤低表达或表达量变化不大的基因；
#由于是差异基因，这里不做过滤；
tmp <- filter.std(dt.f,min.std=0)

#对数据进行标准化；
dt.s <- standardise(tmp)
#查看标准化后的数据；
df.s <- dt.s@assayData$exprs
head(df.s)
#使用mestimate函数估计m值；
m1 <- mestimate(dt.s)
#执行模糊聚类；
set.seed(007)
cl <- mfuzz(dt.s,c=8,m=m1)




#对聚类结果进行可视化；
png("/data/output/report_result/mfuzz.png",,width=1200,height=600)
mfuzz.plot(dt.s,cl,mfrow=c(2,4),
new.window= FALSE,
time.labels=colnames(dt.s))
dev.off()

pdf("/data/output/report_result/mfuzz.pdf",,width=1200,height=600)
mfuzz.plot(dt.s,cl,mfrow=c(2,4),
new.window= FALSE,
time.labels=colnames(dt.s))
dev.off()`
)
