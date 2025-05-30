# -------------------------------------------------------------------------
# Figure 5
#
# Manuscript reference : Section 5.2
# Simulation code reference : simulations/DC_pointcloud.R
# -------------------------------------------------------------------------


## Import main functions
source("simulations/DC_pointcloud.R")

## Generate Figure 5
pdf("Figure5.pdf", width = 20, height = 13)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(2, 3))                 # 2 row, 3 columns


### First Row of Figure 5 : Simulated data clouds on spherical, hyperbolic, and flat manifold, respectively.
col = rainbow_hcl(5)

scatterplot3d::scatterplot3d(df.pc.pos,  pch=19, main="Point cloud data on spherical manifold",zlim=c(-0.4, 1.4), color=col[5],cex.main=2.5,cex.lab=2,cex.axis=2,box = FALSE,xlab="", ylab= "", zlab = "",x.ticklabs="",y.ticklabs="",z.ticklabs="")
scatterplot3d::scatterplot3d(df.pc.neg,  pch=19, main="Point cloud data on hyperbolic manifold", color=col[2],cex.main=2.5,cex.lab=2,cex.axis=2,box = FALSE,xlab="", ylab= "", zlab = "",x.ticklabs="",y.ticklabs="",z.ticklabs="")
scatterplot3d::scatterplot3d(df.pc.flat,  pch=19, main="Point cloud data on flat manifold",zlim=c(-2, 1), color=col[4],cex.main=2.5,cex.lab=2,cex.axis=2,box = FALSE,xlab="", ylab= "", zlab = "",x.ticklabs="",y.ticklabs="",z.ticklabs="")


### Second Row of Figure 5 : Intrinsic curvature test - Plot confidence region for the joint distribution of intrinsic metric and Frechet variance of data clouds on spherical, hyperbolic, and flat manifold, respectively.
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))

par(mar=c(5,6,3,2))
plot(ellipse(cov.est.pc.pos,center = c(Vm.est.pc.pos,Vf.est.pc.pos),level=0.99),type='l',col=col[1],xlim = c(1.6,2.6), ylim =c(1.6,2.6) ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Point cloud data on spherical manifold")
points(Vm.est.pc.pos,Vf.est.pc.pos,pch=16)
lines(ellipse(cov.est.pc.pos,center = c(Vm.est.pc.pos,Vf.est.pc.pos),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.pc.pos,center = c(Vm.est.pc.pos,Vf.est.pc.pos),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

par(mar=c(5,6,3,2))
plot(ellipse(cov.est.pc.neg,center = c(Vm.est.pc.neg,Vf.est.pc.neg),level=0.99),type='l',col=col[1],xlim = c(5,8), ylim =c(5,8),lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Point cloud data on hyperbolic manifold")
points(Vm.est.pc.neg,Vf.est.pc.neg,pch=16)
lines(ellipse(cov.est.pc.neg,center = c(Vm.est.pc.neg,Vf.est.pc.neg),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.pc.neg,center = c(Vm.est.pc.neg,Vf.est.pc.neg),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

par(mar=c(5,6,3,2))
plot(ellipse(cov.est.pc.flat,center = c(Vm.est.pc.flat,Vf.est.pc.flat),level=0.99),type='l',col=col[1],xlim = c(0.4,0.65), ylim =c(0.4,0.65),lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2.5,cex.lab=2,cex.axis=2,main="Point cloud data on flat manifold")
points(Vm.est.pc.flat,Vf.est.pc.flat,pch=16)
lines(ellipse(cov.est.pc.flat,center = c(Vm.est.pc.flat,Vf.est.pc.flat),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.pc.flat,center = c(Vm.est.pc.flat,Vf.est.pc.flat),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

par(old_par)                          # restore original settings
dev.off()