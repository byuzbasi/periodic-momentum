# Render saved pilot CSVs only. No optimizer, QR solve, eigenanalysis or fitting.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L, capabilities("cairo"))
analysis <- normalizePath(args[1],mustWork=TRUE)
out <- args[2]
if(dir.exists(out))stop("Refusing to overwrite manuscript figures")
dir.create(out,recursive=TRUE)
t <- read.csv(file.path(analysis,"trajectories.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(t)==24L*513L)
methods <- c("CHB","HB","GD","NAG-Q")
colors <- c("#176B91","#B45322","#60656F","#157F68")
styles <- c(1,3,4,5)
device <- function(name,height) {
  path <- file.path(out,name)
  if(file.exists(path))stop("Refusing to overwrite a figure")
  cairo_pdf(path,width=8,height=height,family="Arial",pointsize=14,onefile=TRUE)
}
device("Fig2.pdf",5.7)
layout(matrix(c(1,1,2,3,4,5),nrow=3,byrow=TRUE),heights=c(.5,3,3))
par(cex=1,mar=c(0,0,0,0));plot.new()
legend("center",methods,col=colors,lty=styles,lwd=1.7,ncol=4,bty="n",cex=1)
par(cex=1,mar=c(3.4,3.9,1.9,.6),las=1,mgp=c(2.4,.7,0))
panel <- 0L
for(profile in c("exact","scalar_half"))for(case in c("synthetic_R_half","synthetic_R_eighth")) {
  panel <- panel+1L
  z <- t[t$case==case & t$profile==profile,]
  gap <- if(case=="synthetic_R_half")"1/2" else "1/8"
  title <- paste0("(",letters[panel],") R = ",gap,if(profile=="exact")": exact" else ": scalar half")
  plot(NA,xlim=c(0,512),ylim=range(z$log10_D),xlab="Gradient queries k",ylab="",main=title,cex.main=1)
  mtext(expression(log[10](D[k])),side=2,line=2.8,las=0)
  abline(h=0,col="#BFC6CC",lty=3)
  for(i in seq_along(methods)) {
    v <- z[z$method==methods[i],]
    lines(v$k,v$log10_D,col=colors[i],lty=styles[i],lwd=1.6)
  }
}
dev.off()

device("Fig3.pdf",6.1)
layout(matrix(c(1,1,2,3,4,5),nrow=3,byrow=TRUE),heights=c(.8,3,3))
par(cex=1,mar=c(0,0,0,0));plot.new()
legend("center",c("CHB even","CHB odd","HB","GD","NAG-Q"),
       col=c(colors[1],colors[1],colors[2:4]),lty=c(1,2,styles[2:4]),
       pch=c(1,1,NA,NA,NA),lwd=1.6,ncol=3,bty="n",cex=1)
par(cex=1,mar=c(3.4,3.9,1.9,.6),las=1,mgp=c(2.4,.7,0))
panel <- 0L
for(metric in c("log10_D","log10_F"))for(profile in c("exact","scalar_half")) {
  panel <- panel+1L
  z <- t[t$case=="yacht" & t$profile==profile,]
  title <- paste0("(",letters[panel],") ",if(profile=="exact")"Exact gradients" else "Scalar half")
  plot(NA,xlim=c(0,512),ylim=range(z[[metric]],finite=TRUE),xlab="Gradient queries k",ylab="",main=title,cex.main=1)
  mtext(if(metric=="log10_D")expression(log[10](D[k])) else expression(log[10](F[k])),side=2,line=2.8,las=0)
  abline(h=0,col="#94A3B8",lty=3)
  for(parity in 0:1) {
    v <- z[z$method=="CHB" & z$k%%2==parity,]
    lines(v$k,v[[metric]],col=colors[1],lty=parity+1L,lwd=1.6)
    mark <- v[v$k%%64==parity | v$k%in%c(511,512),]
    points(mark$k,mark[[metric]],col=colors[1],pch=1,cex=.55)
  }
  for(i in 2:4) {
    v <- z[z$method==methods[i],]
    lines(v$k,v[[metric]],col=colors[i],lty=styles[i],lwd=1.6)
  }
}
dev.off()
