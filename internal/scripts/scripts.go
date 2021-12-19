package scripts

const (
	FEATURE_SCRIPT = `
if(! require(xlsx))
        install.packages("xlsx")
library(Rsubread)
library(limma)
library(edgeR)
args<-commandArgs(trailingOnly=TRUE)
bamFile<-args[1]
gtfFile<-args[2]
nthreads<-args[3]
outFilePref<-args[4]
outStatsFilePath<-paste(outFilePref,'.stat',sep='');
outCountsFilePath<-paste(outFilePref,'.count',sep='');
outCountsFilePathXls<-paste(outFilePref,'.csv',sep='');
fCountsList=featureCounts(bamFile,annot.ext=gtfFile,isGTFAnnotationFile=TRUE,nthreads=nthreads,isPairedEnd=TRUE)
dgeList=DGEList(counts=fCountsList$counts,genes=fCountsList$annotation)
fpkm=rpkm(dgeList,dgeList$gene$Lenght)
tpm=exp(log(fpkm)-log(sum(fpkm))+log(1e6))
write.table(fCountsList$stat,outStatsFilePath,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
featureCounts=cbind(fCountsList$annotation[,1],fCountsList$counts,fpkm,tpm)
head(featureCounts)
colnames(featureCounts)=c('gene_id','counts','fpkm','tpm')
write.table(featureCounts,outCountsFilePath,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.csv(featureCounts, outCountsFilePathXls,row.names = T)
`
	MATRIX_COUNT = `
if(! require(xlsx))
        install.packages("xlsx")
args<-commandArgs(trailingOnly=TRUE)
inFile<-args[1]
outFile<-args[2]
outFileNumber<-args[3]
outFileXls<-args[4]
a=read.csv(inFile,header=F,sep="\t",stringsAsFactors = F)
colnames(a)=c('sample','gene','counts')
head(a)
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
#counts = counts[-1,]
head(counts)


write.table(counts,file=outFile,sep="\t",quote=FALSE,row.names=FALSE)

names(counts)
names = counts$gene
head(names)
data = as.data.frame(lapply(counts[,-1],as.numeric))
data = cbind(gene = counts[,1],data)
rownames(data) = names
#data$gene = counts$gene
head(data)
data = na.omit(data)
summary(data)
write.table(data,file=outFileNumber,sep="\t",quote=FALSE,row.names=FALSE)

#write.table(counts,file=outFile,sep="\t",quote=FALSE,row.names=FALSE)
write.csv(counts, outFileXls,row.names = T)
`
	MATRIX_TPM = `
if(! require(xlsx))
        install.packages("xlsx")
args<-commandArgs(trailingOnly=TRUE)
inFile<-args[1]
outFile<-args[2]
outFileXls<-args[3]
a=read.csv(inFile,header=F,sep="\t",stringsAsFactors = F)
colnames(a)=c('sample','gene','tpm')
head(a)
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
#counts = counts[-1,]
head(counts)
write.table(counts,file=outFile,sep="\t",quote=FALSE,row.names=FALSE)
write.csv(counts, outFileXls,row.names = T)
`
	MATRIX_FPKM = `
if(! require(xlsx))
        install.packages("xlsx")
args<-commandArgs(trailingOnly=TRUE)
inFile<-args[1]
outFile<-args[2]
outFileXls<-args[3]
a=read.csv(inFile,header=F,sep="\t",stringsAsFactors = F)
colnames(a)=c('sample','gene','fpkm')
head(a)
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
#counts = counts[-1,]
head(counts)
write.table(counts,file=outFile,sep="\t",quote=FALSE,row.names=FALSE)
write.csv(counts, outFileXls,row.names = T)
`

	FAST_REPORT = `
if(! require(xlsx))
        install.packages("xlsx")
rm(list = ls())
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
	READ_REPORT = `args=commandArgs(T)
print(args[1])
print(args[2])
data = read.table(args[1], header = TRUE, stringsAsFactors = FALSE)
library(ggplot2)
data = data[c("Group","Tag_count")]
png(args[2],width=800,height=600)
ggplot(data=data,mapping=aes(x=Group,y=Tag_count,fill=Group,group=factor(1)))+
  geom_bar(stat="identity")
dev.off()
`
	FPKM_REPORT = `
if(! require(xlsx))
        install.packages("xlsx")
rm(list = ls())
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
}`
	HEATMAP_REPORT = `install.packages("pheatmap")
library(pheatmap)
data = read.table("/data/output/expression_result/gene_count.csv",header = T,sep ="\t",row.names = 1)
data = as.data.frame(lapply(data,as.numeric))
data[data==0] <- NA
data[is.na(data)] <- min(data,na.rm = T)*0.01
head(data)
png("/data/output/report_result/cluster.png")
pheatmap(log10(data))
dev.off()

pdf("/data/output/report_result/cluster.pdf")
pheatmap(log10(data))
dev.off()`

	INSERTSECT = `
args=commandArgs(T)
print(args[1])
print(args[2])
sets = strsplit(args[1], ",")[[1]]
res = NULL
for(i in sets)
{
	if (is.null(res)) {
	res = read.table(i, header = TRUE, stringsAsFactors = FALSE)$Gene
}else {
	res = intersect(res,read.table(i, header = TRUE, stringsAsFactors = FALSE)$Gene)
}
}
res
data = data.frame(res)
colnames(data) = c("Gene")
write.table(data,file=paste(args[2],"merge.txt",sep = "_"))
write.csv(data,file=paste(args[2],"merge.csv",sep = "_"),row.names = F)
`
	NOMODE_GO = `
args=commandArgs(T)
print(args[1])
print(args[2])
display_number = c(10, 10, 10)
## GO enrichment with clusterProfiler
if(! require(clusterProfiler))
        install.packages("clusterProfiler")
if(! require("GOplot"))
        install.packages("GOplot")
install.packages("%s", repos = NULL, type = "source")
library(%s)
library(clusterProfiler)
library(GOplot)
glist = read.table(args[1], header = TRUE, stringsAsFactors = FALSE)$Gene
ego_MF <- enrichGO(OrgDb=%s,
             gene = glist,
 	           keyType="GID",
             #pvalueCutoff = %.2f,
			 pvalueCutoff = 1,
			 qvalueCutoff = 1,
             ont = "MF")
#png(paste(args[2],"go_enrich_q_mf_0.05_all.png",sep = "_"))
#goplot(ego_MF)
#pdf(paste(args[2],"go_enrich_q_mf_0.05_all.pdf",sep = "_"))
#goplot(ego_MF)

ego_MF=as.data.frame(ego_MF)
ego_MF=ego_MF[order(ego_MF$Count,decreasing = T),]
ego_MF
ego_result_MF <-  na.omit(ego_MF[1:display_number[1], ])
# ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]
nrow(ego_result_MF)
ego_CC <- enrichGO(OrgDb=%s,
                   gene = glist,
                   #pvalueCutoff = %.2f,
                   pvalueCutoff = 1,
			       qvalueCutoff = 1,
		               keyType="GID",
                   ont = "CC")
                  # readable=TRUE)
#png(paste(args[2],"go_enrich_q_cc_0.05_all.png",sep = "_"))
#goplot(ego_CC)
#pdf(paste(args[2],"go_enrich_q_cc_0.05_all.pdf",sep = "_"))
#goplot(ego_CC)

ego_CC=as.data.frame(ego_CC)
ego_CC=ego_CC[order(ego_CC$Count,decreasing = T),]
ego_result_CC <-  na.omit(ego_CC[1:display_number[2], ])
# ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]
nrow(ego_result_CC)
ego_BP <- enrichGO(OrgDb=%s,
                   gene = glist,
		               keyType="GID",
                   #pvalueCutoff = %.2f,
                   pvalueCutoff = 1,
			       qvalueCutoff = 1,
                   ont = "BP")
                   #readable=TRUE)
#png(paste(args[2],"go_enrich_q_bp_0.05_all.png",sep = "_"))
#goplot(ego_BP)
#pdf(paste(args[2],"go_enrich_q_bp_0.05_all.pdf",sep = "_"))
#goplot(ego_BP)

ego_BP=as.data.frame(ego_BP)
ego_BP=ego_BP[order(ego_BP$Count,decreasing = T),]
ego_result_BP <- na.omit(ego_BP[1:display_number[3], ])
# ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]
nrow(ego_result_BP)
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                                   Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                                   GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                                   type=factor(c(rep("biological process", nrow(ego_result_BP)), rep("cellular component", nrow(ego_result_CC)),
                                          rep("molecular function", nrow(ego_result_MF))), levels=c("biological process", "cellular component","molecular function")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                       collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
nrow(go_enrich_df)
labels=(sapply(
  levels(factor(go_enrich_df$Description)),
  shorten_names))

names(labels)=factor(rev(1:nrow(go_enrich_df)))
#names(labels) = rev(1:nrow(go_enrich_df))
# nrow(go_enrich_df)
# rev(1:nrow(go_enrich_df))

## colors for bar // green, blue, orange
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
pdf(paste(args[2],"go_enrich_q_0.05_all.pdf",sep = "_"))
p
dev.off()
png(paste(args[2],"go_enrich_q_0.05_all.png",sep = "_"))
p
dev.off()
`
	NOMODE_GO_EX = `
args=commandArgs(T)
print(args[1])
print(args[2])
display_number = c(10, 10, 10)
## GO enrichment with clusterProfiler
if(! require(clusterProfiler))
        install.packages("clusterProfiler")
if(! require("GOplot"))
        install.packages("GOplot")
install.packages("%s", repos = NULL, type = "source")
library(%s)
library(clusterProfiler)
library(GOplot)
glist = read.table(args[1], header = TRUE, stringsAsFactors = FALSE)$Gene
ego_MF <- enrichGO(OrgDb=%s,
             gene = glist,
 	           keyType="GID",
             pvalueCutoff = %.2f,
             ont = "MF")
#png(paste(args[2],"go_enrich_q_mf_0.05.png",sep = "_"))
#goplot(ego_MF)
#pdf(paste(args[2],"go_enrich_q_mf_0.05.pdf",sep = "_"))
#goplot(ego_MF)

ego_MF=as.data.frame(ego_MF)
ego_MF=ego_MF[order(ego_MF$Count,decreasing = T),]
ego_MF
ego_result_MF <-  na.omit(ego_MF[1:display_number[1], ])
# ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]
nrow(ego_result_MF)
ego_CC <- enrichGO(OrgDb=%s,
                   gene = glist,
                   pvalueCutoff = %.2f,
		               keyType="GID",
                   ont = "CC")
                  # readable=TRUE)
#png(paste(args[2],"go_enrich_q_cc_0.05.png",sep = "_"))
#goplot(ego_CC)
#pdf(paste(args[2],"go_enrich_q_cc_0.05.pdf",sep = "_"))
#goplot(ego_CC)

ego_CC=as.data.frame(ego_CC)
ego_CC=ego_CC[order(ego_CC$Count,decreasing = T),]
ego_result_CC <-  na.omit(ego_CC[1:display_number[2], ])
# ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]
nrow(ego_result_CC)
ego_BP <- enrichGO(OrgDb=%s,
                   gene = glist,
		               keyType="GID",
                   pvalueCutoff = %.2f,
                   ont = "BP")
                   #readable=TRUE)
#png(paste(args[2],"go_enrich_q_bp_0.05.png",sep = "_"))
#goplot(ego_BP)
#pdf(paste(args[2],"go_enrich_q_bp_0.05.pdf",sep = "_"))
#goplot(ego_BP)

ego_BP=as.data.frame(ego_BP)
ego_BP=ego_BP[order(ego_BP$Count,decreasing = T),]
ego_result_BP <- na.omit(ego_BP[1:display_number[3], ])
# ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]
nrow(ego_result_BP)
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                                   Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                                   GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                                   type=factor(c(rep("biological process", nrow(ego_result_BP)), rep("cellular component", nrow(ego_result_CC)),
                                          rep("molecular function", nrow(ego_result_MF))), levels=c("biological process", "cellular component","molecular function")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                       collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
nrow(go_enrich_df)
labels=(sapply(
  levels(factor(go_enrich_df$Description)),
  shorten_names))

names(labels)=factor(rev(1:nrow(go_enrich_df)))
#names(labels) = rev(1:nrow(go_enrich_df))
# nrow(go_enrich_df)
# rev(1:nrow(go_enrich_df))

## colors for bar // green, blue, orange
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
pdf(paste(args[2],"go_enrich_q_0.05.pdf",sep = "_"))
p
dev.off()
png(paste(args[2],"go_enrich_q_0.05.png",sep = "_"))
p
dev.off()
`
	NOMODO_KEGG = `
library(purrr)
if(! require(clusterProfiler))
        install.packages("clusterProfiler")
if(! require(tidyverse))
        install.packages("tidyverse")

########################################################################################
# clusterProfiler 可能是目前最优秀的富集分析软件，参考网站：
# https://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
# 这里富集分析使用的是通用富集分析方法，参考：http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#universal-enrichment-analysis
########################################################################################

################################################
# 导入自己构建的 OrgDb
################################################
library(%s)
columns(%s)

########################################################################################
# 导入需要进行富集分析的基因列表，并转换为向量
#########################################################################################
#gene_list <- read_csv(Args[6], col_names = FALSE)
args=commandArgs(T)
print(args[1])
print(args[2])


glist = read.table(args[1], header = TRUE, stringsAsFactors = FALSE)$Gene

################################################
# 从 OrgDB 提取 Pathway 和基因的对应关系
################################################

pathway2gene <- AnnotationDbi::select(%s, 
                                      keys = keys(%s), 
                                      columns = c("Pathway")) %%>%%
  na.omit() %%>%%
  dplyr::select(Pathway, GID)

################################################
# 导入 Pathway 与名称对应关系
################################################
# 对json文件操作
if(!file.exists('kegg_info.RData')){
  library(jsonlite)
  library(purrr)
  library(RCurl)

  update_kegg <- function(json = "ko00001.json",file=NULL) {
    pathway2name <- tibble(Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())

    kegg <- fromJSON(json)

    for (a in seq_along(kegg[["children"]][["children"]])) {
      A <- kegg[["children"]][["name"]][[a]]

      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]

        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]

          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %%>%% str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))

          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]

          kos <- str_match(kos_info, "K[0-9]*")[,1]

          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
      }
    }

    save(pathway2name, ko2pathway, file = file)
  }

  update_kegg(json = "ko00001.json",file="kegg_info.RData")

}
load("kegg_info.RData")

#KEGG pathway 富集
ekp <- enricher(glist, 
                TERM2GENE = pathway2gene, 
                TERM2NAME = pathway2name, 
                pvalueCutoff = 1, 
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 5)
ekp_results <- as.data.frame(ekp)
if(dim(ekp_results)[1] == 0){
	return(message("data frame len is zero!"))
}
nrow(ekp_results)
write.table(ekp_results, file = paste(args[2],"DEG.enrichKEGG.txt",sep = "_"))
write.csv(ekp_results, paste(args[2],"DEG.enrichKEGG.csv",sep = "_"),row.names = T)

pdf(paste(args[2],"KEGG_enrichment.barplot_1.pdf",sep = "_"))
barplot(ekp,showCategory=50,drop=T,font.size = 8)
dev.off()

png(paste(args[2],"KEGG_enrichment.barplot_1.png",sep = "_"))
barplot(ekp,showCategory=50,drop=T,font.size = 8)
dev.off()


pdf(paste(args[2],"KEGG_enrichment.dotplot_1.pdf",sep = "_"))
dotplot(ekp,showCategory=50, font.size = 8)
dev.off()

png(paste(args[2],"KEGG_enrichment.dotplot_1.png",sep = "_"))
dotplot(ekp,showCategory=50, font.size = 8)
dev.off()

ekp <- enricher(glist,
                TERM2GENE = pathway2gene,
                TERM2NAME = pathway2name,
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 5)

ekp_results <- as.data.frame(ekp)
if(dim(ekp_results)[1] == 0){
	return(message("data frame len is zero!"))
}
nrow(ekp_results)
write.table(ekp_results, file = paste(args[2],"DEG.enrichKEGG.p.txt",sep = "_"))
write.csv(ekp_results, paste(args[2],"DEG.enrichKEGG.p.csv",sep = "_"),row.names = T)

pdf(paste(args[2],"KEGG_enrichment.barplot_p.pdf",sep = "_"))
barplot(ekp,showCategory=10,drop=T,font.size = 8)
dev.off()

png(paste(args[2],"KEGG_enrichment.barplot_p.png",sep = "_"))
barplot(ekp,showCategory=10,drop=T,font.size = 8)
dev.off()


pdf(paste(args[2],"KEGG_enrichment.dotplot_p.pdf",sep = "_"))
dotplot(ekp,showCategory=10, font.size = 8)
dev.off()

png(paste(args[2],"KEGG_enrichment.dotplot_p.png",sep = "_"))
dotplot(ekp,showCategory=10, font.size = 8)
dev.off()

ekp <- enricher(glist,
                TERM2GENE = pathway2gene,
                TERM2NAME = pathway2name,
                pvalueCutoff = 1,
                qvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                minGSSize = 5)

ekp_results <- as.data.frame(ekp)
if(dim(ekp_results)[1] == 0){
	return(message("data frame len is zero!"))
}
nrow(ekp_results)
write.table(ekp_results, file = paste(args[2],"DEG.enrichKEGG.q.txt",sep = "_"))
write.csv(ekp_results, paste(args[2],"DEG.enrichKEGG.q.csv",sep = "_"),row.names = T)


pdf(paste(args[2],"KEGG_enrichment.barplot_q.pdf",sep = "_"))
barplot(ekp,showCategory=10,drop=T,font.size = 8)
dev.off()

png(paste(args[2],"KEGG_enrichment.barplot_q.png",sep = "_"))
barplot(ekp,showCategory=10,drop=T,font.size = 8)
dev.off()


pdf(paste(args[2],"KEGG_enrichment.dotplot_q.pdf",sep = "_"))
dotplot(ekp,showCategory=10, font.size = 8)
dev.off()

png(paste(args[2],"KEGG_enrichment.dotplot_q.png",sep = "_"))
dotplot(ekp,showCategory=10, font.size = 8)
dev.off()`
	BUILD_STAT = `
library(ggplot2)
data =  read.csv("/data/output/report_result/stat.csv",header = T)
data
png("/data/output/report_result/build_stat.png")
bp<- ggplot(data, aes(x="", y=Percentage, fill=Group))+
geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y")+ theme(panel.grid = element_blank(),axis.text = element_blank(), axis.title = element_blank(),axis.ticks = element_blank()) +guides(fill = guide_legend(title = "Class_code"))
newlegend <- paste(data$Group, " (", round(data$Percentage/sum(data$Percentage)* 100, 2), "%)", sep = "")
bp = bp + scale_fill_discrete(breaks = data$Group, labels = newlegend)
bp
dev.off()
pdf("/data/output/report_result/build_stat.pdf")
bp
dev.off()
`
)
