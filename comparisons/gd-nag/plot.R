# Saved CSVs only: base R/Cairo SVG and PNG; no optimization or data fitting.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L, capabilities("cairo"))
analysis <- normalizePath(args[1],mustWork=TRUE)
out <- args[2]
if (dir.exists(out)) stop("Refusing to overwrite figures")
dir.create(out,recursive=TRUE)
t <- read.csv(file.path(analysis,"trajectories.csv"),stringsAsFactors=FALSE)
p <- read.csv(file.path(analysis,"paired_maxima.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(t)==24L*513L,nrow(p)==24L*256L)
methods <- c("CHB","HB","GD","NAG-Q")
colors <- c("#176B91","#B45322","#60656F","#157F68")
styles <- c(1,2,3,4)
cases <- c("synthetic_R_half","synthetic_R_eighth","yacht")
titles <- c("Synthetic R = 1/2","Synthetic R = 1/8","Yacht")
profiles <- c("exact","scalar_half")
draw <- function(data,metric,xcol,ylabel,yacht_only=FALSE) {
  chosen <- if (yacht_only) 3L else 1:3
  layout(matrix(c(rep(1L,2L),seq_len(2L*length(chosen))+1L),ncol=2L,byrow=TRUE),
         heights=c(.48,rep(3,length(chosen))))
  par(cex=1,mar=c(0,0,0,0));plot.new()
  legend("center",methods,col=colors,lty=styles,lwd=2,ncol=4,bty="n",cex=.95)
  par(cex=1,mar=c(3.8,4.4,2.1,.8),las=1,mgp=c(2.5,.7,0))
  for (ci in chosen) for (profile in profiles) {
    d <- data[data$case==cases[ci] & data$profile==profile,]
    plot(NA,xlim=c(0,512),ylim=range(d[[metric]],finite=TRUE),
         xlab="Gradient queries",ylab="",main=paste(titles[ci],if(profile=="exact")"- exact" else "- scalar half"),cex.main=.95)
    mtext(ylabel,side=2,line=3.1,las=0)
    abline(h=0,lty=3,col="#BFC6CC")
    for (i in seq_along(methods)) {
      v <- d[d$method==methods[i],]
      lines(v[[xcol]],v[[metric]],col=colors[i],lty=styles[i],lwd=1.7)
    }
  }
}
for (name in c("distance","paired_distance","yacht_objective")) {
  is_yacht <- name=="yacht_objective"
  height <- if(is_yacht)4.3 else 9.4
  for (format in c("svg","png")) {
    path <- file.path(out,paste0(name,".",format))
    if (format=="svg") svg(path,width=9,height=height,family="Arial",pointsize=14) else
      png(path,width=1800,height=as.integer(200*height),res=200,type="cairo",family="Arial",pointsize=14)
    if(name=="distance")draw(t,"log10_D","k",expression(log[10](D[k])))
    if(name=="paired_distance")draw(p,"log10_max_D","even_k",expression(log[10](D[pair])))
    if(is_yacht)draw(t,"log10_F","k",expression(log[10](F[k])),TRUE)
    dev.off()
  }
}
writeLines(c(R.version.string,"Base R/Cairo; accepted CSVs only; no numerical optimizer execution."),file.path(out,"render_environment.txt"))
