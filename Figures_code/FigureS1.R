# -------------------------------------------------------------------------
# Figure S.1
#
# Manuscript reference : Section S.1.1
# Simulation code reference : simulations/DC_spd.R
# -------------------------------------------------------------------------


## Import main functions
source("simulations/DC_spd.R")


## Generate Figure S.1
pdf("FigureS1.pdf", width = 20, height = 7)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings
par(mfrow = c(1, 3))                 # 1 row, 3 columns


### Figure S.1(a) : Ratio between metric variance over Frechet variance
par(mar=c(5,5,3,2))
plot(Vf.spd/Vm.spd,col='black',xlab = "Metrics",ylab = "Ratio of Variances",xaxt='n',ylim=c(0.7,1.3),cex=2,cex.axis=2,cex.lab=2,pch=19)
axis(side = 1, at = 1:6, labels = c("Frob", "log-E", "p-Frob", "Chol","AIR", "BW"), cex.axis=2)
legend("topleft",legend = c(expression(V[F]~"/"~V[M])),fill=c("black"),cex=2)


### Figure S.1(b) : Curvature test results with Affine-Invariant Riemannian (AIR) metric
rgPal <- colorRampPalette(c('cadetblue','coral3'))
col = rgPal((3))

par(mar=c(5,5,3,2))
plot(ellipse(cov.est.AIRM,center = c(Vm.est.AIRM,Vf.est.AIRM),level=0.99),type='l',col=col[1],main="Space of SPD matrix with AIR metric",cex.lab=2,cex.axis=2,cex.main=2.5,lwd=2,xlab=expression(V[M]), ylab = expression(V[F]),xlim=c(15,45),ylim=c(15,45))
points(Vm.est.AIRM,Vf.est.AIRM,pch=16)
lines(ellipse(cov.est.AIRM,center = c(Vm.est.AIRM,Vf.est.AIRM),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.AIRM,center = c(Vm.est.AIRM,Vf.est.AIRM),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)


### Figure S.1(c) : Curvature test results with Procrustes size-and-shape (PSS) metric (a.k.a Bures-Wasserstein (BW) metric)
par(mar=c(5,5,3,2))
plot(ellipse(cov.est.PSS,center = c(Vm.est.PSS,Vf.est.PSS),level=0.99),type='l',col=col[1],main="Space of SPD matrix with BW metric",cex.lab=2,cex.axis=2,cex.main=2.5,lwd=2,xlab=expression(V[M]), ylab = expression(V[F]),xlim=c(4500,9500),ylim=c(4500,9500))
points(Vm.est.PSS,Vf.est.PSS,pch=16)
lines(ellipse(cov.est.PSS,center = c(Vm.est.PSS,Vf.est.PSS),level=0.95),col=col[2],lwd=2)
lines(ellipse(cov.est.PSS,center = c(Vm.est.PSS,Vf.est.PSS),level=0.90),col=col[3],lwd=2)
abline(coef = c(0,1),lty=2)
legend("topleft", legend = c(expression(alpha~"= 0.01"),expression(alpha~"= 0.05"),expression(alpha~"= 0.1")),lty=1,lwd=2, col = col,cex=2)

par(old_par)                          # restore original settings
dev.off()