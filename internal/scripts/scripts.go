package scripts

const (
	FEATURE_SCRIPT = `library(Rsubread)
library(limma)
library(edgeR)
args<-commandArgs(trailingOnly=TRUE)
bamFile<-args[1]
gtfFile<-args[2]
nthreads<-args[3]
outFilePref<-args[4]
outStatsFilePath<-paste(outFilePref,'.stat',sep='');
outCountsFilePath<-paste(outFilePref,'.count',sep='');
fCountsList=featureCounts(bamFile,annot.ext=gtfFile,isGTFAnnotationFile=TRUE,nthreads=nthreads,isPairedEnd=TRUE)
dgeList=DGEList(counts=fCountsList$counts,genes=fCountsList$annotation)
fpkm=rpkm(dgeList,dgeList$gene$Lenght)
tpm=exp(log(fpkm)-log(sum(fpkm))+log(1e6))
write.table(fCountsList$stat,outStatsFilePath,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
featureCounts=cbind(fCountsList$annotation[,1],fCountsList$counts,fpkm,tpm)
head(featureCounts)
colnames(featureCounts)=c('gene_id','counts','fpkm','tpm')
write.table(featureCounts,outCountsFilePath,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
`
	MATRIX_COUNT = `args<-commandArgs(trailingOnly=TRUE)
inFile<-args[1]
outFile<-args[2]
a=read.csv(inFile,header=F,sep="\t")
colnames(a)=c('sample','gene','counts')
head(a)
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
counts = counts[-1,]
head(counts)
write.table(counts,file=outFile,sep="\t",quote=FALSE,row.names=FALSE)`
	MATRIX_TPM = `args<-commandArgs(trailingOnly=TRUE)
inFile<-args[1]
outFile<-args[2]
a=read.csv(inFile,header=F,sep="\t")
colnames(a)=c('sample','gene','tpm')
head(a)
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
counts = counts[-1,]
head(counts)
write.table(counts,file=outFile,sep="\t",quote=FALSE,row.names=FALSE)
`
	MATRIX_FPKM = `args<-commandArgs(trailingOnly=TRUE)
inFile<-args[1]
outFile<-args[2]
a=read.csv(inFile,header=F,sep="\t")
colnames(a)=c('sample','gene','fpkm')
head(a)
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
counts = counts[-1,]
head(counts)
write.table(counts,file=outFile,sep="\t",quote=FALSE,row.names=FALSE)
`

	FAST_REPORT = `rm(list = ls())
library(ggplot2)
library(grid)
args <- commandArgs(trailingOnly=TRUE)
name = args[1]
dir = "/data/output/clean/"
filename1 = paste0(dir,name,"-read1BeforeContentCurves.csv")
filename2 = paste0(dir,name,"-read2BeforeContentCurves.csv")

content1 = read.table(filename1,sep = ",",header = T)
content2 = read.table(filename2,sep = ",",header = T)
result1 = NULL
for (i in 2:ncol(content1)) {
    if (i >=6 ) {
        next
    }
    temp = data.frame(index = content1[,1],value = content1[,i],type = colnames(content1)[i],id = "read1")
    result1 = rbind(result1,temp)
}
result2 = NULL
for (i in 2:ncol(content2)) {
    if (i >=6 ) {
        next
    }
    temp = data.frame(index = content2[,1],value = content2[,i],type = colnames(content2)[i],id = "read2")
    result2 = rbind(result2,temp)
}

result = rbind(result1,result2)
filename = paste0("/data/output/report_result/",name,"_content_before.png")
png(filename,width=1200,height=600)
ggplot(result, aes(x=index, y=value, colour=type,group=type)) +
    geom_line(size=0.5)+ facet_wrap(~id,scales="free")+ xlab("position in read(bp)") +
    ylab("value") + theme(panel.background = element_rect(fill = 'grey', colour = 'grey'))
dev.off()

filename1 = paste0(dir,name,"-read1AfterContentCurves.csv")
filename2 = paste0(dir,name,"-read2AfterContentCurves.csv")

content1 = read.table(filename1,sep = ",",header = T)
content2 = read.table(filename2,sep = ",",header = T)
result1 = NULL
for (i in 2:ncol(content1)) {
    if (i >=6 ) {
        next
    }
    temp = data.frame(index = content1[,1],value = content1[,i],type = colnames(content1)[i],id = "read1")
    result1 = rbind(result1,temp)
}
result2 = NULL
for (i in 2:ncol(content2)) {
    if (i >=6 ) {
        next
    }
    temp = data.frame(index = content2[,1],value = content2[,i],type = colnames(content2)[i],id = "read2")
    result2 = rbind(result2,temp)
}

result = rbind(result1,result2)
filename = paste0("/data/output/report_result/",name,"_content_after.png")
png(filename,width=1200,height=600)
ggplot(result, aes(x=index, y=value, colour=type,group=type)) +
    geom_line(size=0.5)+ facet_wrap(~id,scales="free")+ xlab("position in read(bp)") +
    ylab("value") + theme(panel.background = element_rect(fill = 'grey', colour = 'grey'))
dev.off()
################################################################


filename1 = paste0(dir,name,"-read1BeforeQualityCurves.csv")
filename2 = paste0(dir,name,"-read2BeforeQualityCurves.csv")
quality1 = read.table(filename1,sep = ",",header = T)
quality2 = read.table(filename2,sep = ",",header = T)

qualityResult1 = NULL
for (i in 2:ncol(quality1)) {
    temp = data.frame(index = quality1[,1],value = quality1[,i],type = colnames(quality1)[i],id = "read1")
    qualityResult1 = rbind(qualityResult1,temp)
}

qualityResult2 = NULL
for (i in 2:ncol(quality2)) {
    temp = data.frame(index = quality2[,1],value = quality2[,i],type = colnames(quality2)[i],id = "read2")
    qualityResult2 = rbind(qualityResult2,temp)
}

result = rbind(qualityResult1,qualityResult2)
filename = paste0("/data/output/report_result/",name,"_quality_before.png")
png(filename,width=1200,height=600)
ggplot(result, aes(x=index, y=value, colour=type,group=type)) +
    geom_line(size=0.5)+ facet_wrap(~id,scales="free") + xlab("position in read(bp)") +
    ylab("quality") + theme(panel.background = element_rect(fill = 'grey', colour = 'grey'))
dev.off()


filename1 = paste0(dir,name,"-read1AfterQualityCurves.csv")
filename2 = paste0(dir,name,"-read2AfterQualityCurves.csv")
quality1 = read.table(filename1,sep = ",",header = T)
quality2 = read.table(filename2,sep = ",",header = T)

qualityResult1 = NULL
for (i in 2:ncol(quality1)) {
    temp = data.frame(index = quality1[,1],value = quality1[,i],type = colnames(quality1)[i],id = "read1")
    qualityResult1 = rbind(qualityResult1,temp)
}

qualityResult2 = NULL
for (i in 2:ncol(quality2)) {
    temp = data.frame(index = quality2[,1],value = quality2[,i],type = colnames(quality2)[i],id = "read2")
    qualityResult2 = rbind(qualityResult2,temp)
}

result = rbind(qualityResult1,qualityResult2)
filename = paste0("/data/output/report_result/",name,"_quality_after.png")
png(filename,width=1200,height=600)
ggplot(result, aes(x=index, y=value, colour=type,group=type)) +
    geom_line(size=0.5)+ facet_wrap(~id,scales="free") + xlab("position in read(bp)") +
    ylab("quality") + theme(panel.background = element_rect(fill = 'grey', colour = 'grey'))
dev.off()
`
	FPKM_REPORT = `rm(list = ls())
fpkm = read.table("/data/output/expression_result/gene_fpkm.csv",sep = "\t",header = T)
#读取基因表达值数据
#fpkm <- read.delim('gene_fpkm.csv', row.names = 1, sep = '\t')

#过滤标准:保留基因在大于floor(ncol(exprSet)/10)个样本中表达量大于1
#fpkm_filter=fpkm[apply(fpkm,1,function(x) sum(x>1) > floor(ncol(fpkm)/10)),]
dim(fpkm)
fpkm = apply(fpkm,2,as.numeric)
fpkm[is.na(fpkm)] = 0
result = NULL
colnames(fpkm)
for (i in 2:ncol(fpkm)) {
    temp = data.frame(value = fpkm[,i],sample = colnames(fpkm)[i])
    result = rbind(result,temp)
}
#result = as.data.frame(lapply(result,as.numeric))
#result[is.na(result)] = 0
head(result)

result$value = log10(result$value)
library(ggplot2)

pdf("/data/output/report_result/fpkm.pdf")
ggplot(result,aes(value,fill=sample, color=sample)) +
    xlab("log10(fpkm)") +
    geom_density(alpha = 0.3) +
    geom_rug() + theme_bw()
dev.off()

png("/data/output/report_result/fpkm.png",width=1200,height=600)
ggplot(result,aes(value,fill=sample, color=sample)) +
    xlab("log10(fpkm)") +
    geom_density(alpha = 0.3) +
    geom_rug() + theme_bw()
dev.off()
# 添加分面
library(ggpubr)

pdf("/data/output/report_result/fpkm_all.pdf")
ggdensity(result, x = "value",
facet.by = "sample",
rug = TRUE, xlab = "log10(fpkm)",
color = "sample", fill = "sample")
dev.off()
png("/data/output/report_result/fpkm_all.png",width=1200,height=600)
ggdensity(result, x = "value",
facet.by = "sample",
rug = TRUE, xlab = "log10(fpkm)",
color = "sample", fill = "sample")
dev.off()

png("/data/output/report_result/fpkm_violin.png")
ggplot(result, aes(sample,value))+
  geom_violin(aes(fill = sample),trim = FALSE)+
  geom_boxplot(width = 0.2)+
  #scale_fill_manual(values=c(brewer.pal(5,"Set2")[c(1,3,2,5)]))+
  theme_classic()+
  labs(x='Stage',y='The expression level')+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),axis.line=element_line(colour="black",size=0.25),axis.title=element_text(size=5,face="plain",color="black"),axis.text = element_text(size=12,face="plain",color="black"),legend.position="none")
dev.off()

pdf("/data/output/report_result/fpkm_violin.pdf")
ggplot(result, aes(sample,value))+
  geom_violin(aes(fill = sample),trim = FALSE)+
  geom_boxplot(width = 0.2)+
  #scale_fill_manual(values=c(brewer.pal(5,"Set2")[c(1,3,2,5)]))+
  theme_classic()+
  labs(x='Stage',y='The expression level')+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),axis.line=element_line(colour="black",size=0.25),axis.title=element_text(size=5,face="plain",color="black"),axis.text = element_text(size=12,face="plain",color="black"),legend.position="none")
dev.off()
`
	PCA_REPORT = `rm(list = ls())
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
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type",size = 5) + theme_bw()
dev.off()
pdf("/data/output/report_result/pca_count.pdf")
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type",size = 5) + theme_bw()
dev.off()`
	FREQ_REPORT = `rm(list = ls())
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
}`
	HEATMAP_REPORT = `rm(list = ls())
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
dev.off()`
	NOMODE_KEGG = `rm(list = ls())
library("clusterProfiler")
# 导入基因列表
args<-commandArgs(trailingOnly=TRUE)
geneFile<-args[1]
name<-args[2]
geneFile
name
gene <- read.csv(geneFile,header = T,sep="\t")
input <- as.factor(gene$Gene)

# 导入注释文件
term2gene <- read.csv("/data/input/kegg/term2gene.csv",header=T,sep=",")
colnames(term2gene) = c("ko_term","gene")

#input = as.factor(term2gene[100:1000,2])
#term2gene = data.frame(gene = term2gene$Gene,term = term2gene$Term)
term2name <- read.csv("/data/input/kegg/term2name.csv",header=T,sep=",")
colnames(term2name) = c("ko_term","name")


#test = merge(term2gene,term2name,by = "ko_term")
# 富集分析
data = enricher(input,TERM2GENE=term2gene,TERM2NAME = term2name,pvalueCutoff = 0.05, pAdjustMethod = "BH",qvalueCutoff = 1)
#x = enricher(gene,TERM2GENE=term2gene,TERM2NAME = term2name)
head(data)
filterName <- paste0("/data/output/kegg_result/",name,".kegg.data.csv")
write.csv(data,file = filterName)
# 绘制条形图
filterName <- paste0("/data/output/kegg_result/",name,".kegg.bar.pdf")
pdf(filterName)
barplot(data,showCategory = 10)
# 绘制气泡图
filterName <- paste0("/data/output/kegg_result/",name,".kegg.dot.pdf")
pdf(filterName)
dotplot(data)`
)
