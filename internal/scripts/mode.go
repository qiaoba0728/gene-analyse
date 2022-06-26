package scripts

const (
	MODE_GO_MOUSE = `
args=commandArgs(T)
print(args[1])
print(args[2])
display_number = c(10, 10, 10)
## GO enrichment with clusterProfiler
library(%s)
library(clusterProfiler)
library(GOplot)
library(topGO)
library(Rgraphviz)
glist = read.table(args[1], header = TRUE, stringsAsFactors = FALSE)$Gene
ego_MF <- enrichGO(OrgDb=%s,
             gene = glist,
                   keyType="SYMBOL",
             pvalueCutoff = %f,
             ont = "MF")

pdf(paste(args[2],"go_enrich_mf_top10.pdf",sep = "_"),width=800,height=800)
plotGOgraph(ego_MF)
dev.off()
png(paste(args[2],"go_enrich_mf_top10.png",sep = "_"),width=800,height=800)
plotGOgraph(ego_MF)
dev.off()

ego_MF=as.data.frame(ego_MF)
ego_MF=ego_MF[order(ego_MF$Count,decreasing = T),]
ego_MF
name=paste(args[2],"MF_GO.csv",sep = "_")
write.table(ego_MF,file=name,sep = "\t")
ego_result_MF <-  na.omit(ego_MF[1:display_number[1], ])
nrow(ego_result_MF)
ego_result_MF
ego_CC <- enrichGO(OrgDb=%s,
                   gene = glist,
                   pvalueCutoff = %f,
                               keyType="SYMBOL",
                   ont = "CC")
                  # readable=TRUE)


pdf(paste(args[2],"go_enrich_cc_top10.pdf",sep = "_"),width=800,height=800)
plotGOgraph(ego_CC)
dev.off()
png(paste(args[2],"go_enrich_cc_top10.png",sep = "_"),width=800,height=800)
plotGOgraph(ego_CC)
dev.off()


ego_CC=as.data.frame(ego_CC)
ego_CC=ego_CC[order(ego_CC$Count,decreasing = T),]
name=paste(args[2],"CC_GO.csv",sep = "_")
write.table(ego_CC,file=name,sep = "\t")

ego_result_CC <-  na.omit(ego_CC[1:display_number[2], ])
# ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]
nrow(ego_result_CC)
ego_result_CC
ego_BP <- enrichGO(OrgDb=%s,
                   gene = glist,
                               keyType="SYMBOL",
                   pvalueCutoff = %f,
                   ont = "BP")
                   #readable=TRUE)

pdf(paste(args[2],"go_enrich_bp_top10.pdf",sep = "_"),width=800,height=800)
plotGOgraph(ego_BP)
dev.off()
png(paste(args[2],"go_enrich_bp_top10.png",sep = "_"),width=800,height=800)
plotGOgraph(ego_BP)
dev.off()

ego_BP=as.data.frame(ego_BP)
ego_BP=ego_BP[order(ego_BP$Count,decreasing = T),]
name=paste(args[2],"BP_GO.csv",sep = "_")
write.table(ego_BP,file=name,sep = "\t")
ego_result_BP <- na.omit(ego_BP[1:display_number[3], ])
nrow(ego_result_BP)
ego_result_BP
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                                   Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                               GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                                   type=factor(c(rep("molecular function", nrow(ego_result_MF)), rep("cellular component", nrow(ego_result_CC)),rep("biological process", nrow(ego_result_BP))), levels=c( "molecular function","cellular component","biological process")))
name=paste(args[2],"ALL_GO.csv",sep = "_")
write.table(go_enrich_df,file=name,sep = "\t")
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

#labels=(sapply(
#  levels(factor(go_enrich_df$Description)),
#  shorten_names))
nrow(go_enrich_df)
go_enrich_df$Description = lapply(go_enrich_df$Description,shorten_names)
head(go_enrich_df)

#names(labels)=factor(rev(1:nrow(go_enrich_df)))


## colors for bar // green, blue, orange
CPCOLS <- c("#66C3A5","#8DA1CB","#FD8D62")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = CPCOLS) + theme_bw() +
  scale_x_discrete(labels=go_enrich_df$Description) +
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
pdf(paste(args[2],"go_enrich_q_0.05.pdf",sep = "_"),width=800,height=800)
p
dev.off()
png(paste(args[2],"go_enrich_q_0.05.png",sep = "_"),width=800,height=800)
p
dev.off()
`

	MODE_KEGG_MOUSE = `
args=commandArgs(T)
print(args[1])
print(args[2])
library(%s)
library(ggplot2)
library(clusterProfiler)
glist = read.table(args[1], header = TRUE, stringsAsFactors = FALSE)$Gene
mapdt <- bitr(glist, fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb = %s)
data=enrichKEGG(mapdt$ENTREZID,organism = %s,keyType="kegg",,pvalueCutoff = %f,pAdjustMethod = "BH",qvalueCutoff = 0.1)

keggData=as.data.frame(data)
head(keggData)
write.table(keggData, file = paste(args[2],"DEG.enrichKEGG-0.05.txt",sep = "_"))
write.csv(keggData, paste(args[2],"DEG.enrichKEGG-0.05.csv",sep = "_"),row.names = T)

pdf(paste(args[2],"KEGG_enrichment.dotplot-0.05.pdf",sep = "_"))
dotplot(data, showCategory=20,font.size = 8) #气泡图
dev.off()
png(paste(args[2],"KEGG_enrichment.dotplot-0.05.png",sep = "_"))
dotplot(data, showCategory=20,font.size = 8) #气泡图
dev.off()

pdf(paste(args[2],"KEGG_enrichment.barplot-0.05.pdf",sep = "_"))
barplot(data, showCategory=20) #气泡图
dev.off()
png(paste(args[2],"KEGG_enrichment.barplot-0.05.png",sep = "_"))
barplot(data,showCategory=20,drop=T) #柱状图
dev.off()

`
)
