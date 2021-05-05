a=read.csv('gene.count',header=F,sep="\t")
head(a)
colnames(a)=c('sample','gene','counts')
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
counts = counts[-1,]
head(counts)
write.table(counts,file="gene_count.csv",sep="\t",quote=FALSE,row.names=FALSE)