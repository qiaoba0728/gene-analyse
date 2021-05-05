library("clusterProfiler",ask = F,update = F)
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
dotplot(data)