# -------------------------------------------------------------------------
# Figure S7
#
# Sensitivity analysis for the choice of input distance using Gait synchronization data in the space of symmetric positive definite matrices
# Data available at https://github.com/deepcharles/gait-data
# Reference: Truong et al. 2019, 'A data set for the study of human locomotion with inertial measurements units', Image Processing On Line 9, 381–390."
# Manuscript reference : Section S.6
# Data application code reference : Applications/DC_gait.R
# -------------------------------------------------------------------------

## Import main functions
source("Applications/DC_gait.R")


## Generate Figure S.7: Sensitivity analysis for different input distances (Bures-Wasserstein, Frobenius, and Cholesky metric)
pdf("FigureS7.pdf", width = 10, height = 17)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(3, 2),                   # 3 rows × 2 columns, filled row-wise
    mar   = c(5, 5, 7, 2),             # inner margins for every panel
    oma   = c(0, 1, 9, 0))             # 9 outer lines at the top

rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))


## First row of Figure S.7: Intrinsic curvature test with Bures-Wasserstein input distances
plot(ellipse(cov.est.healthy.BW,center = c(Vm.est.healthy.BW,Vf.est.healthy.BW),level=0.99), xlim=c(400, 1900), ylim = c(400,1900), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Healthy Group")
points(Vm.est.healthy.BW,Vf.est.healthy.BW,pch=16)
lines(ellipse(cov.est.healthy.BW,center = c(Vm.est.healthy.BW,Vf.est.healthy.BW),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.healthy.BW,center = c(Vm.est.healthy.BW,Vf.est.healthy.BW),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

plot(ellipse(cov.est.ortho.BW,center = c(Vm.est.ortho.BW,Vf.est.ortho.BW),level=0.99), xlim=c(200, 1700), ylim = c(200,1700), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Orthopedic Group")
points(Vm.est.ortho.BW,Vf.est.ortho.BW,pch=16)
lines(ellipse(cov.est.ortho.BW,center = c(Vm.est.ortho.BW,Vf.est.ortho.BW),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.ortho.BW,center = c(Vm.est.ortho.BW,Vf.est.ortho.BW),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

mtext("Input distance: Bures–Wasserstein metric",side = 3, outer = TRUE, line = 0, font = 2, cex = 1.6)


## Second row of Figure S.7: Intrinsic curvature test with Cholesky input distances
plot(ellipse(cov.est.healthy.Chol,center = c(Vm.est.healthy.Chol,Vf.est.healthy.Chol),level=0.99), xlim=c(800, 2800), ylim = c(800,2800), , type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Healthy Group")
points(Vm.est.healthy.Chol,Vf.est.healthy.Chol,pch=16)
lines(ellipse(cov.est.healthy.Chol,center = c(Vm.est.healthy.Chol,Vf.est.healthy.Chol),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.healthy.Chol,center = c(Vm.est.healthy.Chol,Vf.est.healthy.Chol),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

plot(ellipse(cov.est.ortho.Chol,center = c(Vm.est.ortho.Chol,Vf.est.ortho.Chol),level=0.99), xlim=c(200, 1900), ylim = c(200,1900), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Orthopedic Group")
points(Vm.est.ortho.Chol,Vf.est.ortho.Chol,pch=16)
lines(ellipse(cov.est.ortho.Chol,center = c(Vm.est.ortho.Chol,Vf.est.ortho.Chol),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.ortho.Chol,center = c(Vm.est.ortho.Chol,Vf.est.ortho.Chol),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

mtext("Input distance: Cholesky metric", side = 3, outer = TRUE, line = -41, font = 2, cex = 1.6)


## Third row of Figure S.7: Intrinsic curvature test with Frobenius input distances
plot(ellipse(cov.est.healthy.Fb,center = c(Vm.est.healthy.Fb,Vf.est.healthy.Fb),level=0.99), xlim=c(0, 50000000), ylim = c(0,50000000), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Healthy Group")
points(Vm.est.healthy.Fb,Vf.est.healthy.Fb,pch=16)
lines(ellipse(cov.est.healthy.Fb,center = c(Vm.est.healthy.Fb,Vf.est.healthy.Fb),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.healthy.Fb,center = c(Vm.est.healthy.Fb,Vf.est.healthy.Fb),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

plot(ellipse(cov.est.ortho.Fb,center = c(Vm.est.ortho.Fb,Vf.est.ortho.Fb),level=0.99), xlim=c(-5000000, 39000000), ylim = c(-5000000,39000000), type='l',col=col[1] ,lwd=2,xlab=expression(V[paste(italic(I),",",italic(M))]), ylab = expression(V[paste(italic(I),",",italic(F))]),cex.main=2,cex.lab=2,cex.axis=2,main="Orthopedic Group")
points(Vm.est.ortho.Fb,Vf.est.ortho.Fb,pch=16)
lines(ellipse(cov.est.ortho.Fb,center = c(Vm.est.ortho.Fb,Vf.est.ortho.Fb),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.ortho.Fb,center = c(Vm.est.ortho.Fb,Vf.est.ortho.Fb),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

mtext("Input distance: Frobenius metric", side = 3, outer = TRUE, line = -81, font = 2, cex = 1.6)

par(old_par)                          # restore original settings
dev.off()