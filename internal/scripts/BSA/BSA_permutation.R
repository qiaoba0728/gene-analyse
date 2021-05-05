options(stringsAsFactors=F)
Args <- commandArgs()
dat <- read.table(Args[6])
dat <- dat[!(dat[,4]>0.9 & dat[,5] > 0.9),]

sample_sizeA <- 20
sample_sizeB <- 20

resample_delta <- function(sizeA,sizeB,covA,covB) {
	res <- rep(0,1000)
	for (i in 1:1000) {
		A <- sample(c(0,1),sizeA,replace=T)
		B <- sample(c(0,1),sizeB,replace=T)
		AT <- sample(A,size=covA,replace=T)
		BT <- sample(B,size=covB,replace=T)
		res[i] <- abs((sum(AT) / covA) - (sum(BT) / covB))
	}
	res <- sort(res)
	myres <- c(res[c(950,990)]) 
}

delta_permutation <- matrix(0,nrow=nrow(dat),ncol=2)
## set progress bar ##
pb <- txtProgressBar(style=3)
for (i in 1:nrow(dat)) {
	delta_permutation[i,] <- resample_delta(sample_sizeA,sample_sizeB,as.numeric(dat[i,7]),as.numeric(dat[i,12]))
	setTxtProgressBar(pb, i/nrow(dat))
}
close(pb)

chr <- as.character(dat[,1])
chr <- as.numeric(unlist(strsplit(chr,split="[A-Z,a-z]+"))[c(F,T)])
delta_permutation <- cbind(chr,delta_permutation)

## plot ##
mywindow <- as.numeric(Args[7])
mystep <- as.numeric(Args[8])

#chromosome_len <- 1:9
delta_plot <- cbind(delta_permutation,dat[,c(2,6)])
delta_plot <- delta_plot[delta_plot[,1] > 0,]
delta_plot <- cbind(delta_plot,data[4,5])
chromosome_len <- 1:length(unique(delta_plot[,1]))

write.table(delta_plot,"/data/output/delta_permutation.txt",quote=F,row.names=F,col.names=F,sep="\t")

for (i in 1:length(chromosome_len)) {
	tmp_mat <- delta_plot[delta_plot[,1] == i,]
	chromosome_len[i] <- tmp_mat[nrow(tmp_mat),4]
}
chromosome_len <- c(0,chromosome_len)

format_dat <- NULL
text_pos <- 1:(length(chromosome_len)-1)
for (i in 1:(length(chromosome_len)-1) ) {
	tmp_mat <- rbind(NULL,delta_plot[delta_plot[,1] == i,])
	mystart <- seq(tmp_mat[1,4],tmp_mat[nrow(tmp_mat),4],mystep)
	myend <- seq(tmp_mat[1,4] + mywindow,tmp_mat[nrow(tmp_mat),4],mystep)
	mymin <- min(length(mystart),length(myend))
	mystart <- mystart[1:mymin]
	myend <- myend[1:mymin]
	seq_mat <- cbind(mystart,myend)
	tmp_res <- matrix(0,nrow=nrow(seq_mat),ncol=4)
	for (j in 1:nrow(seq_mat)) {
		tmp <- tmp_mat[tmp_mat[,4] >= seq_mat[j,1] & tmp_mat[,4] <= seq_mat[j,2],]
		if (nrow(tmp) == 0) {
			mymean <- 0
                	mymean05 <- 0.4
                	mymean001 <- 0.5
		}else{
			mymean <- mean(tmp[,5])
			mymean05 <- mean(tmp[,2])
			mymean001 <- mean(tmp[,3])
		}
		mypos <- mean(seq_mat[j,1],seq_mat[j,2])
		tmp_res[j,] <- c(mypos,mymean,mymean05,mymean001)
	}
	tmp_res[,1] <- tmp_res[,1] + sum(chromosome_len[1:i])
	format_dat <- rbind(format_dat,tmp_res)
	text_pos[i] <- ( sum(chromosome_len[1:i]) + sum(chromosome_len[1:i+1]) ) / 2
}

if (length(chromosome_len)-1 > 12) {
	c_n = 0.5
}else{
	c_n = 1
}

pdf("Average_Delta_with_confidence_interval.pdf",width=14,height=3.5)
plot(format_dat[,c(1,2)],pch=20,bty="l",type="n",ylab="Average delta value",xlab="Chromosome",xaxt="n",xaxs="i",yaxs="i",ylim=c(0,max(format_dat[,2],na.rm=T)+0.05))
lines(format_dat[,1],format_dat[,3],col="#EBD3E8")
lines(format_dat[,1],format_dat[,4],col="#D1E9E9")
points(format_dat[,c(1,2)],pch=20)
axis(side=1,labels=paste("Chr",1:(length(chromosome_len)-1),sep=""),at=text_pos,tick=0,cex.axis=c_n)
for (i in 2:(length(chromosome_len)-1)) {
	tmp <- sum(chromosome_len[1:i])
	abline(v=tmp,col="gray",lty=2,lwd=2)
}
dev.off()

write.table(format_dat,file="format_final_dat.txt",sep="\t",quote=F,row.names=F,col.names=F)

pdf("Average_Delta_with_confidence_interval_0.05.pdf",width=14,height=3.5)
plot(format_dat[,c(1,2)],pch=20,bty="l",ylab="Average delta value",xlab="Chromosome",xaxt="n",xaxs="i",yaxs="i",ylim=c(0,max(format_dat[,2],na.rm=T)+0.05))
axis(side=1,labels=paste("Chr",1:(length(chromosome_len)-1),sep=""),at=text_pos,tick=0,cex.axis=c_n)
for (i in 2:(length(chromosome_len)-1)) {
        tmp <- sum(chromosome_len[1:i])
        abline(v=tmp,col="gray",lty=2,lwd=2)
}
lines(format_dat[,1],format_dat[,3],col="#EBD3E8")

dev.off()

pdf("Average_Delta_with_confidence_interval_base0.1.pdf",width=14,height=3.5)
plot(format_dat[,c(1,2)],pch=20,bty="l",ylab="Average delta value",xlab="Chromosome",xaxt="n",xaxs="i",yaxs="i",ylim=c(0.1,max(format_dat[,2],na.rm=T)+0.05))
abline(h=0.1)
axis(side=1,labels=paste("Chr",1:(length(chromosome_len)-1),sep=""),at=text_pos,tick=0,cex.axis=c_n)
for (i in 2:(length(chromosome_len)-1)) {
        tmp <- sum(chromosome_len[1:i])
        abline(v=tmp,col="gray",lty=2,lwd=2)
}
lines(format_dat[,1],format_dat[,3],col="#EBD3E8")
lines(format_dat[,1],format_dat[,4],col="#D1E9E9")

dev.off()


pdf("Average_Delta_with_confidence_interval_base0.1_0.05.pdf",width=14,height=3.5)
plot(format_dat[,c(1,2)],pch=20,bty="l",ylab="Average delta value",xlab="Chromosome",xaxt="n",xaxs="i",yaxs="i",ylim=c(0.1,max(format_dat[,2],na.rm=T)+0.05))
abline(h=0.1)
axis(side=1,labels=paste("Chr",1:(length(chromosome_len)-1),sep=""),at=text_pos,tick=0,cex.axis=c_n)
for (i in 2:(length(chromosome_len)-1)) {
        tmp <- sum(chromosome_len[1:i])
        abline(v=tmp,col="gray",lty=2,lwd=2)
}
lines(format_dat[,1],format_dat[,3],col="#EBD3E8")

dev.off()




