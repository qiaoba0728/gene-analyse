a=read.csv('gene.tpm',header=F,sep="\t")
colnames(a)=c('sample','gene','tpm')
head(a)
library(reshape2)
counts=dcast(a,formula=gene~sample)
head(counts)
counts = counts[-1,]
head(counts)
write.table(counts,file="gene_tpm.csv",sep="\t",quote=FALSE,row.names=FALSE)
