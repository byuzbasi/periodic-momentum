args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
input <- args[1]; output <- args[2]
stopifnot(!dir.exists(output))
d <- read.csv(file.path(input,"trajectories.csv"),check.names=FALSE)
stopifnot(nrow(d)==4104L)
methods <- c("CHB","HB","GD","NAG-Q")
profiles <- c("exact","certified_mixing")
colors <- c("#0072B2","#D55E00","#666666","#009E73")
styles <- c(1,2,3,4)
draw <- function() {
  layout(matrix(c(1,1,2,3,4,5),nrow=3,byrow=TRUE),heights=c(0.07,0.465,0.465))
  par(mar=rep(0,4)); plot.new()
  legend("center",legend=methods,col=colors,lty=styles,lwd=2,horiz=TRUE,bty="n",cex=1.05)
  par(mar=c(4.1,4.6,2.3,0.8),mgp=c(2.6,0.7,0),tcl=-0.3)
  panel <- 0
  for (metric in c("log10_D","log10_F")) for (profile in profiles) {
    panel <- panel+1
    sub <- d[d$profile==profile,]
    ylabel <- if(metric=="log10_D") expression(log[10](D[k])) else expression(log[10](F[k]))
    title <- paste0("(",letters[panel],") ",if(profile=="exact") "Exact gradients" else "Certified mixing errors")
    plot(range(sub$k),range(sub[[metric]]),type="n",xlab="Gradient queries k",ylab=ylabel,main=title,las=1)
    abline(h=pretty(range(sub[[metric]])),col="#E5E7EB",lwd=0.7)
    for(i in seq_along(methods)) {
      z <- sub[sub$method==methods[i],]; z <- z[order(z$k),]
      lines(z$k,z[[metric]],col=colors[i],lty=styles[i],lwd=1.8)
    }
    box()
  }
}
dir.create(output)
png(file.path(output,"istanbul_four_methods.png"),width=1500,height=1050,res=140,type="cairo",pointsize=12)
draw(); invisible(dev.off())
svg(file.path(output,"istanbul_four_methods.svg"),width=10.7,height=7.5,pointsize=12)
draw(); invisible(dev.off())
writeLines(capture.output(sessionInfo()),file.path(output,"render_environment.txt"))
