a=read.csv('gene.count',header=T,sep="\t",stringsAsFactors = F)
head(a)
colnames(a)=c('sample','gene','counts')
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
#counts = counts[-1,]
head(counts)
summary(counts)
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
write.table(data,file="gene_count_number.csv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(counts,file="gene_count.csv",sep="\t",quote=FALSE,row.names=FALSE)