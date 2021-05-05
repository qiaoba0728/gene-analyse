a=read.csv('gene.fpkm',header=F,sep="\t")
colnames(a)=c('sample','gene','fpkm')
head(a)
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
counts = counts[-1,]
head(counts)
write.table(counts,file="gene_fpkm.csv",sep="\t",quote=FALSE,row.names=FALSE)
