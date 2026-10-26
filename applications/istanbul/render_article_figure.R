# Render only saved, validated trajectories. No optimizer or numerical solver.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L,capabilities("cairo"))
d <- read.csv(file.path(args[1],"trajectories.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(d)==4104L)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
out <- file.path(args[2],"Fig3.pdf")
stopifnot(!file.exists(out))
methods <- c("CHB","HB","GD","NAG-Q")
colors <- c("#176B91","#B45322","#60656F","#157F68")
styles <- c(1,3,4,5)
cairo_pdf(out,width=8,height=3.45,family="Arial",pointsize=14,onefile=TRUE)
layout(matrix(c(1,1,2,3),nrow=2,byrow=TRUE),heights=c(.5,3))
par(mar=rep(0,4));plot.new()
legend("center",methods,col=colors,lty=styles,lwd=1.7,ncol=4,bty="n",cex=1)
par(mar=c(3.4,4.8,1.9,.6),las=1,mgp=c(2.4,.7,0))
z <- d[d$profile=="exact",]
for(metric in c("log10_D","log10_F")) {
  label <- if(metric=="log10_D") "(a) Iterate distance" else "(b) Objective gap"
  plot(NA,xlim=c(0,512),ylim=range(z[[metric]]),xlab="Gradient queries k",ylab="",main=label,cex.main=1)
  mtext(if(metric=="log10_D")expression(log[10](D[k])) else expression(log[10](F[k])),side=2,line=2.8,las=0)
  abline(h=0,col="#BFC6CC",lty=3)
  for(i in seq_along(methods)) {
    v <- z[z$method==methods[i],];v <- v[order(v$k),]
    lines(v$k,v[[metric]],col=colors[i],lty=styles[i],lwd=1.6)
  }
}
invisible(dev.off())
