# Title     : TODO
# Objective : TODO
# Created by: Administrator
# Created on: 2021/3/16 0016

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

