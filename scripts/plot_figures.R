# Render the frozen reporting CSVs only; no optimizer, QR or spectrum calculation.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, capabilities("cairo"))
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
if (dir.exists(out)) stop("Refusing to overwrite an existing figure directory")
dir.create(out, recursive = TRUE)
synthetic <- file.path(root, "periodic_momentum_numerics_v06", "analysis")
s <- read.csv(file.path(synthetic, "scale.csv"), stringsAsFactors = FALSE)
t <- read.csv(file.path(synthetic, "trajectories.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(s) == 5L, nrow(t) == 16L * 513L)
blue <- "#176B91"; orange <- "#B45322"
device <- function(name, height) {
  grDevices::cairo_pdf(file.path(out, name), width = 8, height = height,
                       family = "Arial", pointsize = 14, onefile = TRUE)
}
device("Fig1.pdf", 3.8)
par(mfrow = c(1,2), cex = 1, mar = c(3.6,4.2,2.0,.6), las = 1, mgp = c(2.3,.7,0))
yl <- range(c(s$ratio,s$limit))
plot(s$R,s$ratio, log = "x", xaxt = "n", pch = 19, type = "b", col = blue,
     xlab = "Normalized gap R", ylab = "", main = "(a) Scaled upper bound",
     cex.main = 1, ylim = yl + c(-1,1) * diff(yl) * .12)
mtext(expression(epsilon[U]/R^4),side = 2,line = 2.9,las = 0)
axis(1,at = rev(s$R),labels = rev(s$R_label))
abline(h = s$limit[1],lty = 2,col = orange,lwd = 1.5)
legend("topleft",c("Scalar formula","Limit"),col = c(blue,orange),
       lty = c(1,2),pch = c(19,NA),bty = "n",cex = 1)
e <- s[!is.na(s$local_exponent),]
yl <- range(c(4,e$local_exponent))
plot(e$R,e$local_exponent,log = "x",xaxt = "n",pch = 19,type = "b",col = blue,
     xlab = "Smaller R in each pair",ylab = "",main = "(b) Local exponent",
     cex.main = 1,ylim = yl + c(-1,1) * diff(yl) * .15)
mtext("Local log ratio exponent",side = 2,line = 2.9,las = 0)
axis(1,at = rev(e$R),labels = rev(e$R_label))
abline(h = 4,col = orange,lty = 2,lwd = 1.5)
dev.off()

device("Fig2.pdf",6.4)
layout(matrix(c(1,1,2,3,4,5),nrow = 3,byrow = TRUE),heights = c(.45,3,3))
# layout() reduces cex for three rows; restore it explicitly for legibility.
par(cex = 1, mar = c(0,0,0,0));plot.new()
legend("center",c("Two-periodic CHB","Stationary HB"),col = c(blue,orange),
       lty = c(1,2),lwd = 1.5,bty = "n",ncol = 2,cex = 1)
par(cex = 1, mar = c(3.5,4.3,1.8,.5),las = 1,mgp = c(2.3,.7,0))
panel <- 0L
for (profile in c("exact","scalar_above")) for (gap in c("1/2","1/8")) {
  z <- t[t$profile == profile & t$R_label == gap,]
  panel <- panel + 1L
  title <- if (profile == "exact")
    bquote(.(paste0("(",letters[panel],")"))~R == .(gap)*": exact") else
    bquote(.(paste0("(",letters[panel],")"))~R == .(gap)*": "*1.5*epsilon[U])
  plot(NA,xlim = c(0,512),ylim = range(z$log10_x_over_qk),xlab = "Gradient queries k",
       ylab = "",main = title,cex.main = 1)
  mtext(expression(log[10](paste("||",x[k],"|| / (||",x[0],"|| ",q^k,")"))),
        side = 2,line = 2.9,las = 0)
  abline(h = 0,col = "#CBD5E1",lty = 3)
  for (method in c("CHB","HB")) {
    zz <- z[z$method == method,]
    lines(zz$k,zz$log10_x_over_qk,col = if (method == "CHB") blue else orange,
          lty = if (method == "CHB") 1 else 2,lwd = 1.5)
  }
}
dev.off()

writeLines(R.version.string,file.path(out,"render_environment.txt"))
