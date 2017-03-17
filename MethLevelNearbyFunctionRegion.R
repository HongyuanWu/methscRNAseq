setwd("/media/Home_Raid1/zhl002/NAS1/projects/mouse_WGBS/bedfiles/BED")
files=Sys.glob("*enhancer*txt")
data<-read.table(files[1]);
main="enhancer"
xlab=paste("Distance to center of"," main(bp)",sep="")
ylab="Methylation level(Beta)"
pdf("enhancer.dist.density.pdf")
plot(data[,2],data[,3],type = "n",ylim=c(0,1),main="enhancer",xlab = xlab,ylab=ylab)
for(i in 1:length(files)){
  data<-read.table(files[i]);
  fit<-smooth.spline(data[,2], data[,3],spar=0.1)
  lines(fit,col=i,lwd=2)
}
legend=unlist(lapply(files,function(x) unlist(strsplit(x,"[.]"))[3]))
legend("bottomright",legend =legend,col=1:length(files),lty=1,cex=0.7,bty = "n")
dev.off()

# bivalentdomain
files=Sys.glob("*bivalentdomain*txt")
files
data<-read.table(files[1]);
main="Bivalentdomain"
xlab=paste("Distance to center of",main," (bp)",sep="")
ylab="Methylation level(Beta)"
pdf("bivalentdomain.nearby.density.pdf")
plot(data[,2],data[,3],type = "n",ylim=c(0,1),main=main,xlab = xlab,ylab=ylab)
for(i in 1:length(files)){
  data<-read.table(files[i]);
  fit<-smooth.spline(data[,2], data[,3])
  lines(fit,col=i,lwd=2)
}

legend=unlist(lapply(files,function(x) unlist(strsplit(x,"[.]"))[3]))
legend("bottomright",legend =legend,col=1:length(files),lty=1,cex=0.7,bty = "n")
dev.off()


files=Sys.glob("*Promoter*txt")
files
data<-read.table(files[1]);
main="Promoter"
xlab=paste("Distance to center of",main," (bp)",sep="")
ylab="Methylation level(Beta)"
pdf("Promoter.nearby.density.pdf")
plot(data[,2],data[,3],type = "n",ylim=c(0,1),main=main,xlab = xlab,ylab=ylab)
for(i in 1:length(files)){
  data<-read.table(files[i]);
  fit<-smooth.spline(data[,2], data[,3])
  lines(fit,col=i,lwd=2)
}
legend=unlist(lapply(files,function(x) unlist(strsplit(x,"[.]"))[4]))
legend("bottomright",legend =legend,col=1:length(files),lty=1,cex=0.7,bty = "n")
dev.off()


files=Sys.glob("*Exon*txt")
files
data<-read.table(files[1]);
main="Exon"
xlab=paste("Distance to center of",main," (bp)",sep="")
ylab="Methylation level(Beta)"
pdf("Exon.nearby.density.pdf")
plot(data[,2],data[,3],type = "n",ylim=c(0,1),main=main,xlab = xlab,ylab=ylab)
for(i in 1:length(files)){
  data<-read.table(files[i]);
  fit<-smooth.spline(data[,2], data[,3])
  lines(fit,col=i,lwd=2)
}
legend=unlist(lapply(files,function(x) unlist(strsplit(x,"[.]"))[4]))
legend("bottomright",legend =legend,col=1:length(files),lty=1,cex=0.7,bty = "n")
dev.off()



files=Sys.glob("*Intron*txt")
files
data<-read.table(files[1]);
main="Intron"
xlab=paste("Distance to center of",main," (bp)",sep="")
ylab="Methylation level(Beta)"
pdf("Intron.nearby.density.pdf")
plot(data[,2],data[,3],type = "n",ylim=c(0,1),main=main,xlab = xlab,ylab=ylab)
for(i in 1:length(files)){
  data<-read.table(files[i]);
  fit<-smooth.spline(data[,2], data[,3])
  lines(fit,col=i,lwd=2)
}
legend=unlist(lapply(files,function(x) unlist(strsplit(x,"[.]"))[4]))
legend("bottomright",legend =legend,col=1:length(files),lty=1,cex=0.7,bty = "n")
dev.off()


files=Sys.glob("*CpGI*txt")
files
data<-read.table(files[1]);
main="CpGI"
xlab=paste("Distance to center of",main," (bp)",sep="")
ylab="Methylation level(Beta)"
pdf("CpGI.nearby.density.pdf")
plot(data[,2],data[,3],type = "n",ylim=c(0,1),main=main,xlab = xlab,ylab=ylab)
for(i in 1:length(files)){
  data<-read.table(files[i]);
  fit<-smooth.spline(data[,2], data[,3])
  lines(fit,col=i,lwd=2)
}
legend=unlist(lapply(files,function(x) unlist(strsplit(x,"[.]"))[4]))
legend("bottomright",legend =legend,col=1:length(files),lty=1,cex=0.7,bty = "n")
dev.off()

